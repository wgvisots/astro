/**********************************************************************
 * close_approach.cpp
 *
 * Copyright 2024 W. Greg Visotski
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 *         limitations under the License.
 **********************************************************************/

#include "astro/orbit/close_approach.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <optional>
#include <tuple>
#include <utility>

#include "../core/enum_io.h"
#include "astro/core/constants.h"
#include "astro/core/error.h"
#include "astro/core/rootfind.h"
#include "astro/orbit/two_body.h"

namespace astro {

std::ostream& operator<<(std::ostream& os, TrajectoryModel value) {
    return detail::writeEnum(os, "TrajectoryModel", value);
}

namespace {

/// Shortest sampling step (s), so that a trajectory through the centre (|r| = 0) still advances.
constexpr double kMinStep = 1e-3;

/// The SecularOrder of a secular model; throws astro::Error for any other value.
SecularOrder secularOrder(TrajectoryModel model) {
    switch (model) {
    case TrajectoryModel::SecularJ2:   return SecularOrder::J2;
    case TrajectoryModel::SecularJ2J4: return SecularOrder::J2SquaredJ4;
    case TrajectoryModel::TwoBody:     break;
    }
    throw Error("Trajectory: invalid TrajectoryModel");
}

/// Two-body propagation from s0: seconds after the epoch -> state.
std::function<StateVector(double)> twoBody(const StateVector& s0, double mu) {
    return [s0, mu](double dt) { return propagateKepler(s0, dt, mu); };
}

/// Secular propagation with velocity = d(position)/dt, which the range rate needs (see
/// SecularPropagatorCorrected). Shared, so copies of the function stay cheap.
std::function<StateVector(double)> secular(std::shared_ptr<const SecularPropagatorCorrected> p) {
    return [p = std::move(p)](double dt) { return p->propagate(dt); };
}

[[nodiscard]] bool isFinite(const StateVector& s) noexcept {
    return std::isfinite(s.r.x) && std::isfinite(s.r.y) && std::isfinite(s.r.z) && std::isfinite(s.v.x)
        && std::isfinite(s.v.y) && std::isfinite(s.v.z);
}

/// Both states at s seconds after the window start, and the range rate there.
struct Sample {
    double      s{0.0}; // seconds after the window start
    StateVector a;      // first object
    StateVector b;      // second object
    double      g{0.0}; // (rb - ra) . (vb - va): negative while closing (m^2/s)
};

/// The sampled range-rate search shared by closestApproach and closeApproaches.
class Search {
public:
    Search(const Trajectory& first, const Trajectory& second, const Time& start, const CloseApproachOptions& opt)
        : fa_(first.secondsAfter(start))
        , fb_(second.secondsAfter(start))
        , opt_(opt) {}

    [[nodiscard]] Sample at(double s) const {
        Sample x{s, fa_(s), fb_(s), 0.0};
        if (!isFinite(x.a) || !isFinite(x.b)) {
            throw Error("closestApproach: a trajectory returned a non-finite state");
        }
        x.g = dot(x.b.r - x.a.r, x.b.v - x.a.v);
        return x;
    }

    /// Samples [0, T] and calls onMinimum with each refined local minimum, in time order. Returns
    /// the samples at 0 and T.
    template <class OnMinimum>
    std::pair<Sample, Sample> run(double T, OnMinimum&& onMinimum) const {
        const Sample first = at(0.0);
        Sample       prev  = first;
        while (prev.s < T) {
            Sample cur = at(std::min(prev.s + step(prev), T));
            if (prev.g < 0.0 && cur.g >= 0.0) { // closing, then separating: a minimum in (prev, cur]
                const double tca = brentRoot([this](double s) { return at(s).g; }, prev.s, cur.s, opt_.timeTolerance);
                onMinimum(at(tca));
            }
            prev = std::move(cur);
        }
        return {first, prev};
    }

private:
    /// |r| / |v|: the time scale on which the object's position turns (infinite at rest).
    [[nodiscard]] static double timeScale(const StateVector& s) noexcept {
        const double v = norm(s.v);
        return v > 0.0 ? norm(s.r) / v : std::numeric_limits<double>::infinity();
    }

    [[nodiscard]] double step(const Sample& x) const noexcept {
        return std::max(opt_.stepFraction * std::min(timeScale(x.a), timeScale(x.b)), kMinStep);
    }

    std::function<StateVector(double)> fa_;
    std::function<StateVector(double)> fb_;
    CloseApproachOptions               opt_;
};

[[nodiscard]] CloseApproach approach(const Sample& x, const Time& start, bool atBoundary) {
    CloseApproach c;
    c.time          = start + x.s;
    c.distance      = norm(x.b.r - x.a.r);
    c.relativeSpeed = norm(x.b.v - x.a.v);
    c.first         = x.a;
    c.second        = x.b;
    c.atBoundary    = atBoundary;
    return c;
}

/// Window length (s), after checking the window and the options.
[[nodiscard]] double windowLength(const Time& start, const Time& end, const CloseApproachOptions& opt) {
    if (!(opt.stepFraction > 0.0) || !std::isfinite(opt.stepFraction)) {
        throw Error("closestApproach: stepFraction must be positive and finite");
    }
    if (!(opt.timeTolerance > 0.0) || !std::isfinite(opt.timeTolerance)) {
        throw Error("closestApproach: timeTolerance must be positive and finite");
    }
    const double T = end - start;
    if (!(T >= 0.0) || !std::isfinite(T)) {
        throw Error("closestApproach: end must not be before start");
    }
    return T;
}

} // namespace

Trajectory::Trajectory(const Time& epoch, const StateVector& state, TrajectoryModel model, const ZonalParameters& zp)
    : epoch_(epoch) {
    if (model == TrajectoryModel::TwoBody) {
        fromEpoch_ = twoBody(state, zp.mu);
    } else {
        fromEpoch_ = secular(std::make_shared<const SecularPropagatorCorrected>(state, zp, secularOrder(model)));
    }
}

Trajectory::Trajectory(const Time&              epoch,
                       const KeplerianElements& elements,
                       TrajectoryModel          model,
                       const ZonalParameters&   zp)
    : epoch_(epoch) {
    if (model == TrajectoryModel::TwoBody) {
        fromEpoch_ = twoBody(elementsToState(elements, zp.mu), zp.mu);
    } else {
        fromEpoch_ = secular(std::make_shared<const SecularPropagatorCorrected>(elements, zp, secularOrder(model)));
    }
}

Trajectory::Trajectory(const Time&                epoch,
                       const EquinoctialElements& elements,
                       TrajectoryModel            model,
                       const ZonalParameters&     zp)
    : Trajectory(epoch, equinoctialToKeplerian(elements), model, zp) {}

Trajectory::Trajectory(StateFunction stateAt)
    : stateAt_(std::move(stateAt)) {
    if (!stateAt_) {
        throw Error("Trajectory: empty state function");
    }
}

StateVector Trajectory::state(const Time& t) const {
    return fromEpoch_ ? fromEpoch_(t - epoch_) : stateAt_(t);
}

std::function<StateVector(double)> Trajectory::secondsAfter(const Time& origin) const {
    if (fromEpoch_) {
        return [f = fromEpoch_, offset = origin - epoch_](double s) { return f(offset + s); };
    }
    return [f = stateAt_, origin](double s) { return f(origin + s); };
}

CloseApproach closestApproach(const Trajectory&           first,
                              const Trajectory&           second,
                              const Time&                 start,
                              const Time&                 end,
                              const CloseApproachOptions& options) {
    const double                 T = windowLength(start, end, options);
    const Search                 search(first, second, start, options);
    std::optional<CloseApproach> best;
    const auto                   consider = [&best](CloseApproach c) {
        if (!best || c.distance < best->distance) {
            best = std::move(c);
        }
    };
    // Minima first, so a true closest approach wins a tie with a window end.
    const auto [s0, sT] = search.run(T, [&](const Sample& x) { consider(approach(x, start, false)); });
    consider(approach(s0, start, true));
    consider(approach(sT, start, true));
    return *best;
}

std::vector<CloseApproach> closeApproaches(const Trajectory&           first,
                                           const Trajectory&           second,
                                           const Time&                 start,
                                           const Time&                 end,
                                           double                      maxDistance,
                                           const CloseApproachOptions& options) {
    const double T = windowLength(start, end, options);
    if (std::isnan(maxDistance)) {
        throw Error("closeApproaches: maxDistance is NaN");
    }
    const Search               search(first, second, start, options);
    std::vector<CloseApproach> out;
    (void)search.run(T, [&](const Sample& x) {
        CloseApproach c = approach(x, start, false);
        if (c.distance <= maxDistance) {
            out.push_back(std::move(c));
        }
    });
    return out;
}

namespace {

/// 0.5 (erf b - erf a) for a <= b. Short intervals by 8-point Gauss-Legendre on the integrand
/// exp(-t^2) / sqrt(pi) (a difference of error functions would cancel); longer ones in
/// complementary form when both ends are on the same side, so that small probabilities keep their
/// relative precision.
[[nodiscard]] double erfInterval(double a, double b) noexcept {
    if (b - a < 0.5) {
        // clang-format off
        static constexpr double x[4] = {0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363};
        static constexpr double w[4] = {0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763};
        // clang-format on
        const double m = 0.5 * (a + b), h = 0.5 * (b - a);
        double       sum = 0.0;
        for (int k = 0; k < 4; ++k) {
            const double t1 = m - h * x[k], t2 = m + h * x[k];
            sum += w[k] * (std::exp(-t1 * t1) + std::exp(-t2 * t2));
        }
        return h * sum / std::sqrt(constants::pi);
    }
    if (a > 0.0) {
        return 0.5 * (std::erfc(a) - std::erfc(b));
    }
    if (b < 0.0) {
        return 0.5 * (std::erfc(-b) - std::erfc(-a));
    }
    return 0.5 * (std::erf(b) - std::erf(a));
}

/// One Gauss-Kronrod 7-15 step on [a, b] (QUADPACK's QK15): the Kronrod estimate, and its
/// difference from the embedded Gauss rule as the error estimate.
template <class F>
std::pair<double, double> kronrod15(const F& f, double a, double b) {
    // clang-format off
    static constexpr double xk[8] = {0.991455371120812639206854697526329, 0.949107912342758524526189684047851,
                                     0.864864423359769072789712788640926, 0.741531185599394439863864773280788,
                                     0.586087235467691130294144845693013, 0.405845151377397166906606412076961,
                                     0.207784955007898467600689403773245, 0.0};
    static constexpr double wk[8] = {0.022935322010529224963732008058970, 0.063092092629978553290700663189204,
                                     0.104790010322250183839876322541518, 0.140653259715525918745189590510238,
                                     0.169004726639267902826583426598550, 0.190350578064785409913256402421014,
                                     0.204432940075298892414161999234649, 0.209482141084727828012999174891714};
    static constexpr double wg[4] = {0.129484966168869693270611432679082, 0.279705391489276667901467771423780,
                                     0.381830050505118944950369775488975, 0.417959183673469387755102040816327};
    // clang-format on
    const double m = 0.5 * (a + b), h = 0.5 * (b - a);
    const double fc = f(m);
    double       k = wk[7] * fc, g = wg[3] * fc;
    for (int i = 0; i < 7; ++i) {
        const double pair  = f(m - h * xk[i]) + f(m + h * xk[i]);
        k                 += wk[i] * pair;
        if (i % 2 == 1) { // the Gauss nodes are the odd Kronrod ones
            g += wg[i / 2] * pair;
        }
    }
    return {k * h, std::abs((k - g) * h)};
}

/// Adaptive Gauss-Kronrod on [a, b] to absolute tolerance `tol` on |Kronrod - Gauss|, which
/// overstates the Kronrod estimate's error by orders of magnitude for smooth integrands.
template <class F>
double adaptiveKronrod(const F& f, double a, double b, double tol, int depth) {
    const auto [value, error] = kronrod15(f, a, b);
    if (depth <= 0 || error <= tol) {
        return value;
    }
    const double m = 0.5 * (a + b);
    return adaptiveKronrod(f, a, m, 0.5 * tol, depth - 1) + adaptiveKronrod(f, m, b, 0.5 * tol, depth - 1);
}

/// Integral of f over [a, b] to ~1e-12 relative: Gauss-Kronrod on 8 panels gives the scale (and
/// makes sure no feature narrower than the interval is missed), then each panel is refined.
template <class F>
double integrate(const F& f, double a, double b) {
    constexpr int kPanels = 8;
    const double  h       = (b - a) / kPanels;
    double        value[kPanels], error[kPanels], coarse = 0.0;
    for (int k = 0; k < kPanels; ++k) {
        std::tie(value[k], error[k])  = kronrod15(f, a + h * k, k == kPanels - 1 ? b : a + h * (k + 1));
        coarse                       += value[k];
    }
    const double tol = std::max(1e-11 * std::abs(coarse), 1e-300) / kPanels;
    double       sum = 0.0;
    for (int k = 0; k < kPanels; ++k) {
        const double lo = a + h * k, hi = k == kPanels - 1 ? b : a + h * (k + 1);
        sum += error[k] <= tol ? value[k]
                               : adaptiveKronrod(f, lo, 0.5 * (lo + hi), 0.5 * tol, 30)
                                     + adaptiveKronrod(f, 0.5 * (lo + hi), hi, 0.5 * tol, 30);
    }
    return sum;
}

/// Principal axes of a 2x2 covariance: variances (major, minor) and the major axis angle.
struct Principal {
    double major{0.0}, minor{0.0}, angle{0.0};
};

[[nodiscard]] Principal principalAxes(double xx, double xy, double yy) noexcept {
    const double mean = 0.5 * (xx + yy);
    const double d    = std::hypot(0.5 * (xx - yy), xy);
    return {mean + d, std::max(mean - d, 0.0), 0.5 * std::atan2(2.0 * xy, xx - yy)};
}

} // namespace

double probabilityInCircle(double mx, double my, double xx, double xy, double yy, double radius) {
    if (!(radius > 0.0) || !std::isfinite(radius)) {
        throw Error("probabilityInCircle: radius must be positive and finite");
    }
    if (!std::isfinite(mx) || !std::isfinite(my) || !std::isfinite(xx) || !std::isfinite(xy) || !std::isfinite(yy)) {
        throw Error("probabilityInCircle: non-finite mean or covariance");
    }
    const double scale = std::max(std::abs(xx), std::abs(yy));
    if (xx < 0.0 || yy < 0.0 || xx * yy - xy * xy < -1e-12 * scale * scale) {
        throw Error("probabilityInCircle: covariance is not positive semi-definite");
    }
    const Principal p = principalAxes(xx, xy, yy);
    const double    c = std::cos(p.angle), s = std::sin(p.angle);
    const double    u0 = c * mx + s * my;  // mean along the major axis
    const double    v0 = -s * mx + c * my; // mean along the minor axis
    const double    su = std::sqrt(p.major), sv = std::sqrt(p.minor);
    const double    R2 = radius * radius;

    if (su == 0.0) { // no uncertainty at all
        return mx * mx + my * my < R2 ? 1.0 : 0.0;
    }
    constexpr double kSqrt2 = 1.4142135623730951;
    if (sv == 0.0) { // all the uncertainty along one line, at offset v0 from the centre
        const double half = v0 * v0 < R2 ? std::sqrt(R2 - v0 * v0) : 0.0;
        return erfInterval((-half - u0) / (kSqrt2 * su), (half - u0) / (kSqrt2 * su));
    }

    // Major axis u = R sin t over the part of the circle within 15 sigma of the mean (outside it the
    // density is below exp(-112)). The chord at u has half-length h = R cos t.
    constexpr double kWidth = 15.0;
    const double     ua     = std::max(-radius, u0 - kWidth * su);
    const double     ub     = std::min(radius, u0 + kWidth * su);
    if (!(ua < ub)) {
        return 0.0;
    }
    const double ta     = std::asin(std::clamp(ua / radius, -1.0, 1.0));
    const double tb     = std::asin(std::clamp(ub / radius, -1.0, 1.0));
    const double norm_u = 1.0 / (std::sqrt(constants::twoPi) * su);
    const auto   f      = [&](double t) {
        const double u = radius * std::sin(t), h = radius * std::cos(t);
        const double z = (u - u0) / su;
        return h * norm_u * std::exp(-0.5 * z * z) * erfInterval((-h - v0) / (kSqrt2 * sv), (h - v0) / (kSqrt2 * sv));
    };
    return std::clamp(integrate(f, ta, tb), 0.0, 1.0);
}

CollisionProbability collisionProbability(const CloseApproach& approach,
                                          const Mat3&          covariance1,
                                          const Mat3&          covariance2,
                                          double               hardBodyRadius) {
    if (!(hardBodyRadius > 0.0) || !std::isfinite(hardBodyRadius)) {
        throw Error("collisionProbability: hardBodyRadius must be positive and finite");
    }
    const Vec3   dr = approach.relativePosition();
    const Vec3   dv = approach.relativeVelocity();
    const double sp = norm(dv);
    if (!(sp > 0.0)) {
        throw Error("collisionProbability: zero relative velocity, no encounter plane");
    }
    CollisionProbability out;
    out.zAxis        = dv / sp;
    const Vec3 m     = dr - dot(dr, out.zAxis) * out.zAxis; // miss vector in the encounter plane
    out.missDistance = norm(m);
    if (out.missDistance > 0.0) {
        out.xAxis = m / out.missDistance;
    } else { // any direction in the plane
        const Vec3 helper = std::abs(out.zAxis.x) < 0.9 ? Vec3{1.0, 0.0, 0.0} : Vec3{0.0, 1.0, 0.0};
        out.xAxis         = normalized(cross(helper, out.zAxis));
    }
    out.yAxis = cross(out.zAxis, out.xAxis);

    Mat3 C  = covariance1;
    C      += covariance2;
    out.xx  = dot(out.xAxis, C * out.xAxis);
    out.yy  = dot(out.yAxis, C * out.yAxis);
    out.xy  = 0.5 * (dot(out.xAxis, C * out.yAxis) + dot(out.yAxis, C * out.xAxis)); // symmetric part
    if (!std::isfinite(out.xx) || !std::isfinite(out.xy) || !std::isfinite(out.yy)) {
        throw Error("collisionProbability: non-finite covariance");
    }
    try {
        out.probability = probabilityInCircle(out.missDistance, 0.0, out.xx, out.xy, out.yy, hardBodyRadius);
    } catch (const Error&) {
        throw Error("collisionProbability: covariance is not positive semi-definite in the encounter plane");
    }
    const Principal p = principalAxes(out.xx, out.xy, out.yy);
    out.sigmaMajor    = std::sqrt(p.major);
    out.sigmaMinor    = std::sqrt(p.minor);
    const double det  = out.xx * out.yy - out.xy * out.xy;
    out.mahalanobis   = det > 0.0               ? out.missDistance * std::sqrt(out.yy / det)
                      : out.missDistance == 0.0 ? 0.0
                                                : std::numeric_limits<double>::infinity();
    return out;
}

CollisionProbability collisionProbability(const CloseApproach& approach,
                                          const Mat6&          covariance1,
                                          const Mat6&          covariance2,
                                          double               hardBodyRadius) {
    return collisionProbability(approach, block(covariance1, 0, 0), block(covariance2, 0, 0), hardBodyRadius);
}

namespace {

/// The approach with both covariances and its collision probability.
ConjunctionAssessment assess(const CloseApproach&       approach,
                             const UncertainTrajectory& first,
                             const UncertainTrajectory& second,
                             double                     hardBodyRadius) {
    ConjunctionAssessment out;
    out.approach         = approach;
    out.firstCovariance  = first.stateCovariance(approach.time).covariance;
    out.secondCovariance = second.stateCovariance(approach.time).covariance;
    out.collision        = collisionProbability(approach, out.firstCovariance, out.secondCovariance, hardBodyRadius);
    return out;
}

void checkHardBodyRadius(double hardBodyRadius) {
    if (!(hardBodyRadius > 0.0) || !std::isfinite(hardBodyRadius)) {
        throw Error("assessConjunction: the hard-body radius must be positive and finite");
    }
}

} // namespace

UncertainTrajectory::UncertainTrajectory(const Time&            epoch,
                                         const StateVector&     state,
                                         const Mat6&            covariance,
                                         TrajectoryModel        model,
                                         const ZonalParameters& zp)
    : trajectory_(epoch, state, model, zp) {
    for (const auto& row : covariance) {
        for (double x : row) {
            if (!std::isfinite(x)) {
                throw Error("UncertainTrajectory: covariance is not finite");
            }
        }
    }
    if (model == TrajectoryModel::TwoBody) {
        covarianceAt_ = [epoch, state, covariance, mu = zp.mu](const Time& t) {
            return propagateKeplerWithCovariance(state, covariance, t - epoch, mu);
        };
    } else {
        auto p        = std::make_shared<const SecularPropagatorCorrected>(state, zp, secularOrder(model));
        covarianceAt_ = [epoch, covariance, p = std::move(p)](const Time& t) {
            return p->propagateWithCovariance(t - epoch, covariance);
        };
    }
}

UncertainTrajectory::UncertainTrajectory(Trajectory trajectory, StateCovarianceFunction stateCovarianceAt)
    : trajectory_(std::move(trajectory))
    , covarianceAt_(std::move(stateCovarianceAt)) {
    if (!covarianceAt_) {
        throw Error("UncertainTrajectory: empty state-covariance function");
    }
}

StateCovariance UncertainTrajectory::stateCovariance(const Time& t) const {
    return covarianceAt_(t);
}

ConjunctionAssessment assessConjunction(const UncertainTrajectory&  first,
                                        const UncertainTrajectory&  second,
                                        const Time&                 start,
                                        const Time&                 end,
                                        double                      hardBodyRadius,
                                        const CloseApproachOptions& options) {
    checkHardBodyRadius(hardBodyRadius);
    return assess(closestApproach(first.trajectory(), second.trajectory(), start, end, options),
                  first,
                  second,
                  hardBodyRadius);
}

std::vector<ConjunctionAssessment> assessConjunctions(const UncertainTrajectory&  first,
                                                      const UncertainTrajectory&  second,
                                                      const Time&                 start,
                                                      const Time&                 end,
                                                      double                      hardBodyRadius,
                                                      double                      maxDistance,
                                                      const CloseApproachOptions& options) {
    checkHardBodyRadius(hardBodyRadius);
    std::vector<ConjunctionAssessment> out;
    for (const CloseApproach& c :
         closeApproaches(first.trajectory(), second.trajectory(), start, end, maxDistance, options)) {
        out.push_back(assess(c, first, second, hardBodyRadius));
    }
    return out;
}

} // namespace astro
