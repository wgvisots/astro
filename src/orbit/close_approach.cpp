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
#include <limits>
#include <memory>
#include <optional>
#include <tuple>
#include <utility>
#include <vector>

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
    double      s{0.0};     // seconds after the window start
    StateVector a;          // first object
    StateVector b;          // second object
    double      g{0.0};     // (rb - ra) . (vb - va): negative while closing (m^2/s)
    double      noise{0.0}; // |g| below this is rounding error of the two states (m^2/s)
};

/// Rounding units of g that count as noise. Objects that keep their distance (the same orbit, a
/// fixed phase apart) have g = 0 up to ~5 units, measured for two-body and secular trajectories
/// from LEO to GEO; at a real minimum the samples either side are many orders of magnitude above it.
constexpr double kRangeRateNoise = 64.0;

/// The sampled range-rate search shared by closestApproach and closeApproaches.
class Search {
public:
    Search(const Trajectory& first, const Trajectory& second, const Time& start, const CloseApproachOptions& opt)
        : fa_(first.secondsAfter(start))
        , fb_(second.secondsAfter(start))
        , opt_(opt) {}

    [[nodiscard]] Sample at(double s) const {
        Sample x{s, fa_(s), fb_(s), 0.0, 0.0};
        if (!isFinite(x.a) || !isFinite(x.b)) {
            throw Error("closestApproach: a trajectory returned a non-finite state");
        }
        const Vec3 dr = x.b.r - x.a.r, dv = x.b.v - x.a.v;
        x.g     = dot(dr, dv);
        x.noise = kRangeRateNoise * std::numeric_limits<double>::epsilon()
                * ((norm(x.a.r) + norm(x.b.r)) * norm(dv) + (norm(x.a.v) + norm(x.b.v)) * norm(dr));
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
            // Closing, then separating: a minimum in (prev, cur]. A change of sign within the
            // rounding error of both samples is not one: the distance is constant there.
            if (prev.g < 0.0 && cur.g >= 0.0 && (-prev.g > prev.noise || cur.g > cur.noise)) {
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

Trajectory::Trajectory(const Time& epoch, const StateVector& state, TrajectoryModel model, const CentralBody& body)
    : epoch_(epoch) {
    if (model == TrajectoryModel::TwoBody) {
        fromEpoch_ = twoBody(state, body.mu);
    } else {
        fromEpoch_ = secular(std::make_shared<const SecularPropagatorCorrected>(state, body, secularOrder(model)));
    }
}

Trajectory::Trajectory(const Time&              epoch,
                       const KeplerianElements& elements,
                       TrajectoryModel          model,
                       const CentralBody&       body)
    : epoch_(epoch) {
    if (model == TrajectoryModel::TwoBody) {
        fromEpoch_ = twoBody(elementsToState(elements, body.mu), body.mu);
    } else {
        fromEpoch_ = secular(std::make_shared<const SecularPropagatorCorrected>(elements, body, secularOrder(model)));
    }
}

Trajectory::Trajectory(const Time&                epoch,
                       const EquinoctialElements& elements,
                       TrajectoryModel            model,
                       const CentralBody&         body)
    : Trajectory(epoch, equinoctialToKeplerian(elements), model, body) {}

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

/// 0.5 (erf b - erf a) for a <= b, keeping its relative precision for small probabilities.
/// Where exp(-t^2) changes by less than a factor e over the interval, a difference of error
/// functions would cancel: 8-point Gauss-Legendre on exp(-t^2) / sqrt(pi) is exact to rounding
/// there. Elsewhere the complementary forms are used when both ends are on the same side.
[[nodiscard]] double erfInterval(double a, double b) noexcept {
    if (b - a < 0.5 && (b - a) * (std::abs(a) + std::abs(b)) < 1.0) {
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

/// One Gauss-Kronrod 7-15 step on [a, b] (QUADPACK's QK15, Piessens et al. 1983): the Kronrod
/// estimate and QUADPACK's error estimate, which sharpens |Kronrod - Gauss| with the integral of
/// |f - mean| and never goes below 50 rounding units of the integral of |f|.
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
    double       lo[7], hi[7];
    double       k = wk[7] * fc, g = wg[3] * fc, absSum = wk[7] * std::abs(fc);
    for (int i = 0; i < 7; ++i) {
        lo[i]   = f(m - h * xk[i]);
        hi[i]   = f(m + h * xk[i]);
        k      += wk[i] * (lo[i] + hi[i]);
        absSum += wk[i] * (std::abs(lo[i]) + std::abs(hi[i]));
        if (i % 2 == 1) { // the Gauss nodes are the odd Kronrod ones
            g += wg[i / 2] * (lo[i] + hi[i]);
        }
    }
    const double mean = 0.5 * k;
    double       asc  = wk[7] * std::abs(fc - mean);
    for (int i = 0; i < 7; ++i) {
        asc += wk[i] * (std::abs(lo[i] - mean) + std::abs(hi[i] - mean));
    }
    asc          *= std::abs(h);
    double error  = std::abs((k - g) * h);
    if (asc != 0.0 && error != 0.0) {
        error = asc * std::min(1.0, std::pow(200.0 * error / asc, 1.5));
    }
    constexpr double eps = std::numeric_limits<double>::epsilon();
    return {k * h, std::max(error, 50.0 * eps * absSum * std::abs(h))};
}

/// Adaptive Gauss-Kronrod on [a, b] for a non-negative integrand. A panel is accepted when its
/// error estimate is within `relTol` of its own value (which bounds the relative error of the sum)
/// or below `absTol`, which is halved at each bisection. It is also accepted when halving it no
/// longer reduces an already small error estimate: the integrand's own rounding noise has been
/// reached, and going on would cost 2^depth evaluations for nothing. `budget` (Kronrod steps left)
/// bounds the cost whatever the integrand.
template <class F>
double adaptiveKronrod(const F& f,
                       double   a,
                       double   b,
                       double   value,
                       double   error,
                       double   relTol,
                       double   absTol,
                       int      depth,
                       int&     budget) {
    if (depth <= 0 || budget <= 0 || error <= std::max(relTol * std::abs(value), absTol)) {
        return value;
    }
    const double m               = 0.5 * (a + b);
    const auto [leftV, leftE]    = kronrod15(f, a, m);
    const auto [rightV, rightE]  = kronrod15(f, m, b);
    budget                      -= 2;
    if (error <= 1e-6 * std::abs(value) && leftE + rightE >= 0.5 * error) {
        return leftV + rightV;
    }
    return adaptiveKronrod(f, a, m, leftV, leftE, relTol, 0.5 * absTol, depth - 1, budget)
         + adaptiveKronrod(f, m, b, rightV, rightE, relTol, 0.5 * absTol, depth - 1, budget);
}

/// Integral of the non-negative f over the panels between consecutive `breaks` (ascending), to
/// ~1e-12 relative. The first pass over every panel gives the scale below which a panel does not
/// matter; the breaks must be placed so that this pass cannot step over a narrow feature. At most
/// ~2000 Kronrod steps (30 000 evaluations) are spent.
template <class F>
double integrate(const F& f, const std::vector<double>& breaks) {
    constexpr double                       kRelTol = 1e-12;
    std::vector<std::pair<double, double>> first(breaks.size() - 1);
    double                                 coarse = 0.0;
    for (std::size_t k = 0; k + 1 < breaks.size(); ++k) {
        first[k]  = kronrod15(f, breaks[k], breaks[k + 1]);
        coarse   += first[k].first;
    }
    const double absTol = kRelTol * std::abs(coarse) / static_cast<double>(first.size());
    double       sum    = 0.0;
    int          budget = 2000;
    for (std::size_t k = 0; k + 1 < breaks.size(); ++k) {
        sum +=
            adaptiveKronrod(f, breaks[k], breaks[k + 1], first[k].first, first[k].second, kRelTol, absTol, 40, budget);
    }
    return sum;
}

/// a * b - c * d without cancellation error (Kahan's algorithm with fused multiply-add).
[[nodiscard]] double differenceOfProducts(double a, double b, double c, double d) noexcept {
    const double cd = c * d;
    return std::fma(a, b, -cd) + std::fma(-c, d, cd);
}

/// Principal axes of a 2x2 covariance: variances (major, minor) and the major axis angle.
struct Principal {
    double major{0.0}, minor{0.0}, angle{0.0};
};

[[nodiscard]] Principal principalAxes(double xx, double xy, double yy) noexcept {
    // The major variance is a sum of positive terms. The minor one, as mean - d, would lose
    // (major / minor) rounding units; determinant / major does not.
    const double major = 0.5 * (xx + yy) + std::hypot(0.5 * (xx - yy), xy);
    const double minor = major > 0.0 ? std::max(differenceOfProducts(xx, yy, xy, xy) / major, 0.0) : 0.0;
    return {major, minor, 0.5 * std::atan2(2.0 * xy, xx - yy)};
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
    if (xx < 0.0 || yy < 0.0 || differenceOfProducts(xx, yy, xy, xy) < -1e-12 * scale * scale) {
        throw Error("probabilityInCircle: covariance is not positive semi-definite");
    }
    const Principal p = principalAxes(xx, xy, yy);
    const double    c = std::cos(p.angle), s = std::sin(p.angle);
    const double    u0 = c * mx + s * my;            // mean along the major axis
    const double    v0 = std::abs(-s * mx + c * my); // mean along the minor axis (the sign does not matter)
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
    // h - v0 = (R - v0) - 2 R sin^2(t / 2): formed this way it keeps its relative precision where
    // the chord and the offset nearly cancel, which a narrow minor axis magnifies by 1 / sv.
    const double gap     = radius - v0;
    const auto   overlap = [&](double t) { // h(t) - v0
        const double sh = std::sin(0.5 * t);
        return gap - 2.0 * radius * sh * sh;
    };
    const auto f = [&](double t) {
        const double u = radius * std::sin(t), h = radius * std::cos(t);
        const double z = (u - u0) / su;
        return h * norm_u * std::exp(-0.5 * z * z) * erfInterval((-h - v0) / (kSqrt2 * sv), overlap(t) / (kSqrt2 * sv));
    };

    // Break points. The chord factor depends on t through x = (h - v0) / sv only: 1 for x >> 0, an
    // error-function step around x = 0, and a Gaussian tail exp(-x^2 / 2) for x << 0. When the minor
    // axis is narrow (sv << R) that is a sharp feature in t, which equal panels could step over, so
    // the range is cut where the factor changes character, working down from the longest chord in
    // the range (x = xMax, at the t nearest 0): at x = 8, 0 and -8 around the step, or, when even
    // the longest chord is in the tail (xMax < -8), at the depths where the tail has fallen by
    // exp(-8) and exp(-32). Eight equal panels then resolve the Gaussian along the major axis.
    std::vector<double> breaks{ta, tb};
    const double        tPeak = std::clamp(0.0, ta, tb);
    const double        xMax  = overlap(tPeak) / sv;
    const auto          cutAt = [&](double x) {            // both t where the chord factor's argument is x
        const double q = (gap - x * sv) / radius; // 1 - cos t there
        if (q > 0.0 && q < 1.0) {
            const double t = 2.0 * std::asin(std::sqrt(0.5 * q));
            for (double tc : {t, -t}) {
                if (tc > ta && tc < tb) {
                    breaks.push_back(tc);
                }
            }
        }
    };
    if (xMax < -8.0) {
        cutAt(xMax - 8.0 / -xMax);
        cutAt(xMax - 32.0 / -xMax);
    } else {
        cutAt(8.0);
        cutAt(0.0);
        cutAt(-8.0);
    }
    if (tPeak > ta && tPeak < tb) {
        breaks.push_back(tPeak);
    }
    constexpr int kPanels = 8;
    for (int k = 1; k < kPanels; ++k) {
        breaks.push_back(ta + (tb - ta) * k / kPanels);
    }
    std::sort(breaks.begin(), breaks.end());
    breaks.erase(std::unique(breaks.begin(), breaks.end()), breaks.end());
    return std::clamp(integrate(f, breaks), 0.0, 1.0);
}

CollisionProbability collisionProbability(const CloseApproach& approach,
                                          const Mat3&          covariance1,
                                          const Mat3&          covariance2,
                                          double               hardBodyRadius) {
    if (!(hardBodyRadius > 0.0) || !std::isfinite(hardBodyRadius)) {
        throw Error("collisionProbability: hardBodyRadius must be positive and finite");
    }
    if (!isFinite(approach.first) || !isFinite(approach.second)) {
        throw Error("collisionProbability: non-finite state in the approach");
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
    const double det  = differenceOfProducts(out.xx, out.yy, out.xy, out.xy);
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

UncertainTrajectory::UncertainTrajectory(const Time&        epoch,
                                         const StateVector& state,
                                         const Mat6&        covariance,
                                         TrajectoryModel    model,
                                         const CentralBody& body)
    : trajectory_(epoch, state, model, body) {
    for (const auto& row : covariance) {
        for (double x : row) {
            if (!std::isfinite(x)) {
                throw Error("UncertainTrajectory: covariance is not finite");
            }
        }
    }
    if (model == TrajectoryModel::TwoBody) {
        covarianceAt_ = [epoch, state, covariance, mu = body.mu](const Time& t) {
            return propagateKeplerWithCovariance(state, covariance, t - epoch, mu);
        };
    } else {
        auto p        = std::make_shared<const SecularPropagatorCorrected>(state, body, secularOrder(model));
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
