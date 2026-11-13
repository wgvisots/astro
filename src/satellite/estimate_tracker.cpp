/**********************************************************************
 * estimate_tracker.cpp
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

#include "astro/satellite/estimate_tracker.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>

#include "../core/enum_io.h"
#include "astro/core/angle.h"
#include "astro/core/error.h"
#include "astro/orbit/elements.h"
#include "astro/satellite/orbital_frames.h"

namespace astro {

std::ostream& operator<<(std::ostream& os, ConsistencyTest value) {
    return detail::writeEnum(os, "ConsistencyTest", value);
}

double chiSquareSurvival3(double x) noexcept {
    if (std::isnan(x)) {
        return x;
    }
    if (!(x > 0.0)) {
        return 1.0;
    }
    // Once exp(-x/2) underflows (x > 1490) the tail is below the smallest double. Returning 0 there
    // also keeps a huge or infinite x from giving infinity * 0 = NaN.
    const double e = std::exp(-0.5 * x);
    if (e == 0.0) {
        return 0.0;
    }
    return std::erfc(std::sqrt(0.5 * x)) + std::sqrt(2.0 * x / constants::pi) * e;
}

double chiSquareSurvival6(double x) noexcept {
    if (std::isnan(x)) {
        return x;
    }
    if (!(x > 0.0)) {
        return 1.0;
    }
    const double h = 0.5 * x;
    const double e = std::exp(-h);
    if (e == 0.0) { // as in chiSquareSurvival3: the polynomial overflows from x = 2.7e154
        return 0.0;
    }
    return e * (1.0 + h + 0.5 * h * h);
}

NumericalPropagator makePropagator(const StateEstimate& e, const TrackerOptions& options) {
    const Mat3 pole = poleOfDateMatrix(options.frame, e.epoch);
    return NumericalPropagator(e.epoch,
                               e.state,
                               {std::make_shared<GravityField>(options.J, options.mu, options.re, pole)},
                               options.integrator);
}

namespace {

/// Squared Mahalanobis distance between `next` and the prediction in equinoctial elements, with
/// the predicted covariance mapped through the element Jacobian. Retrograde orbits are first
/// mirrored (180 deg about x) so that the elements stay far from their i = 180 deg singularity.
/// Also fills the semi-major axis residual and sigma.
std::optional<double> elementsDistance(const StateCovariance& pred, StateVector next, double mu, ConsistencyCheck& c) {
    StateVector p = pred.state;
    Mat6        P = pred.covariance;
    if (cross(p.r, p.v).z < 0.0) {
        const Mat3 M{1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0};
        p    = {M * p.r, M * p.v};
        next = {M * next.r, M * next.v};
        P    = rotateCovariance(P, M);
    }
    const Vec3 wn = normalized(cross(next.r, next.v));
    if (wn.z < -0.5) {
        // After the mirroring the prediction has i <= 90 deg; a new state with i > 120 deg is beyond
        // any plausible covariance (and near the element singularity): infinitely far.
        c.semiMajorAxisResidual = 0.0;
        c.semiMajorAxisSigma    = 0.0;
        return std::numeric_limits<double>::infinity();
    }
    // Equinoctial elements exist for ellipses only. A prediction that is not elliptic cannot be tested
    // this way at all: stateToEquinoctial throws for it.
    const EquinoctialElements ep = stateToEquinoctial(p, mu);
    // A new state on an open orbit (v^2 r >= 2 mu) is, like a reversed one, beyond any plausible
    // covariance of an elliptic prediction: infinitely far.
    if (!(norm2(next.v) * norm(next.r) < 2.0 * mu)) {
        c.semiMajorAxisResidual = 0.0;
        c.semiMajorAxisSigma    = 0.0;
        return std::numeric_limits<double>::infinity();
    }
    const EquinoctialElements   en = stateToEquinoctial(next, mu);
    const Mat6                  Pe = congruence(stateToEquinoctialJacobian(p, mu), P);
    const std::array<double, 6> d{en.a - ep.a,
                                  en.h - ep.h,
                                  en.k - ep.k,
                                  en.p - ep.p,
                                  en.q - ep.q,
                                  wrapPi(en.lambda - ep.lambda)};
    c.semiMajorAxisResidual = d[0];
    c.semiMajorAxisSigma    = Pe[0][0] > 0.0 ? std::sqrt(Pe[0][0]) : 0.0;
    return mahalanobis2(d, Pe, 6);
}

ConsistencyCheck compare(const StateCovariance& pred, const StateEstimate& next, double dt, const TrackerOptions& opt) {
    ConsistencyCheck c;
    c.dt                  = dt;
    c.predicted           = pred.state;
    c.predictedCovariance = pred.covariance;

    const Vec3 dr         = next.state.r - pred.state.r;
    const Vec3 dv         = next.state.v - pred.state.v;
    const Mat3 R          = ricMatrix(pred.state);
    c.positionResidualRic = R * dr;
    c.velocityResidualRic = R * dv;
    const Mat6 Pric       = rotateCovariance(pred.covariance, R);
    c.positionSigmaRic    = {std::sqrt(Pric[0][0]), std::sqrt(Pric[1][1]), std::sqrt(Pric[2][2])};
    c.velocitySigmaRic    = {std::sqrt(Pric[3][3]), std::sqrt(Pric[4][4]), std::sqrt(Pric[5][5])};

    const std::array<double, 6> d{dr.x, dr.y, dr.z, dv.x, dv.y, dv.z};
    const auto                  d6 = mahalanobis2(d, pred.covariance, 6);
    const auto                  d3 = mahalanobis2(d, pred.covariance, 3);
    const auto                  de = elementsDistance(pred, next.state, opt.mu, c);
    c.valid                        = d6.has_value() && d3.has_value() && de.has_value();
    if (c.valid) {
        c.mahalanobis2         = *d6;
        c.positionMahalanobis2 = *d3;
        c.elementsMahalanobis2 = *de;
        c.pValue               = chiSquareSurvival6(*d6);
        c.positionPValue       = chiSquareSurvival3(*d3);
        c.elementsPValue       = chiSquareSurvival6(*de);
        const double p         = opt.test == ConsistencyTest::Position         ? c.positionPValue
                               : opt.test == ConsistencyTest::PositionVelocity ? c.pValue
                                                                               : c.elementsPValue;
        c.consistent           = p >= opt.alarmProbability;
    } else {
        c.mahalanobis2 = c.positionMahalanobis2 = c.elementsMahalanobis2 = std::nan("");
        c.pValue = c.positionPValue = c.elementsPValue = std::nan("");
    }
    return c;
}

/// `e` (in frame f) expressed in GCRS.
StateEstimate toGcrs(const StateEstimate& e, Frame f) {
    if (f == Frame::GCRS) {
        return e;
    }
    const Mat3 Mt = gcrsToCelestialMatrix(f, e.epoch).transpose();
    return {
        e.epoch,
        StateVector{Mt * e.state.r, Mt * e.state.v},
        rotateCovariance(e.covariance, Mt)
    };
}

/// A GCRS state and covariance expressed in the frame with GCRS -> frame rotation M.
StateCovariance rotated(const StateCovariance& s, const Mat3& M) {
    return {
        StateVector{M * s.state.r, M * s.state.v},
        rotateCovariance(s.covariance, M)
    };
}

/// GCRS -> frame rotations over [t0, t1] seconds after `epoch`: exact at nodes at most an hour
/// apart, linear in between. The of-date frames turn by under 1e-7 rad per hour, so the blend stays
/// orthonormal to ~1e-14 and accurate to well under a millimetre in LEO.
class FrameSeries {
public:
    FrameSeries(Frame f, const Time& epoch, double t0, double t1)
        : t0_(t0) {
        const double span = std::max(t1 - t0, 1.0);
        const int    n    = static_cast<int>(std::ceil(span / 3600.0));
        h_                = span / n;
        for (int k = 0; k <= n; ++k) {
            nodes_.push_back(gcrsToCelestialMatrix(f, epoch + (t0 + h_ * k)));
        }
    }

    [[nodiscard]] Mat3 at(double t) const {
        const double      u    = (t - t0_) / h_;
        const double      last = static_cast<double>(nodes_.size() - 2);
        const std::size_t k    = static_cast<std::size_t>(std::clamp(std::floor(u), 0.0, last));
        const double      w    = u - static_cast<double>(k);
        return nodes_[k] * (1.0 - w) + nodes_[k + 1] * w;
    }

private:
    double            t0_;
    double            h_{1.0};
    std::vector<Mat3> nodes_;
};

/// Throws astro::Error unless the estimate's state, covariance and epoch are all finite.
void requireFinite(const StateEstimate& e) {
    bool ok = std::isfinite(e.epoch.jd().jd1) && std::isfinite(e.epoch.jd().jd2);
    for (std::size_t i = 0; i < 3; ++i) {
        ok = ok && std::isfinite(e.state.r[i]) && std::isfinite(e.state.v[i]);
    }
    for (const auto& row : e.covariance) {
        for (double x : row) {
            ok = ok && std::isfinite(x);
        }
    }
    if (!ok) {
        throw Error("EstimateTracker: the estimate has a non-finite state, covariance or epoch");
    }
}

void requireCelestial(Frame f) {
    if (isEarthFixed(f)) {
        throw Error("EstimateTracker: TrackerOptions::frame must be a celestial frame (not TIRS, PEF or ITRS)");
    }
}

StateCovariance propagateTo(const NumericalPropagator& prop, const Mat6& P0, double dt, double q) {
    std::vector<StateCovariance> out(1);
    prop.propagateGridWithCovariance(dt, 1.0, P0, out, q); // a single output: integrates exactly to dt
    return out[0];
}

} // namespace

ConsistencyCheck
checkConsistency(const StateEstimate& prior, const StateEstimate& next, const TrackerOptions& options) {
    requireCelestial(options.frame);
    requireFinite(prior);
    requireFinite(next);
    TrackerOptions gcrs          = options;
    gcrs.frame                   = Frame::GCRS;
    const StateEstimate   priorG = toGcrs(prior, options.frame);
    const double          dt     = next.epoch - prior.epoch;
    const StateCovariance pred =
        propagateTo(makePropagator(priorG, gcrs), priorG.covariance, dt, options.accelerationNoisePsd);
    // Compare in the delivery frame at the new epoch.
    const StateCovariance predF =
        options.frame == Frame::GCRS ? pred : rotated(pred, gcrsToCelestialMatrix(options.frame, next.epoch));
    ConsistencyCheck c = compare(predF, next, dt, options);
    c.replaced         = dt >= 0.0;
    return c;
}

EstimateTracker::EstimateTracker(TrackerOptions options)
    : opt_(std::move(options)) {
    requireCelestial(opt_.frame);
}

std::optional<ConsistencyCheck> EstimateTracker::update(const StateEstimate& next) {
    requireFinite(next);
    TrackerOptions gcrs                   = opt_;
    gcrs.frame                            = Frame::GCRS;
    const StateEstimate             nextG = toGcrs(next, opt_.frame);
    std::optional<ConsistencyCheck> check;
    if (current_) {
        const double          dt   = next.epoch - current_->epoch;
        const StateCovariance pred = propagateTo(*prop_, currentGcrs_->covariance, dt, opt_.accelerationNoisePsd);
        // Compare in the delivery frame at the new epoch.
        const StateCovariance predF =
            opt_.frame == Frame::GCRS ? pred : rotated(pred, gcrsToCelestialMatrix(opt_.frame, next.epoch));
        check           = compare(predF, next, dt, opt_);
        check->replaced = dt >= 0.0;
        if (!check->replaced) {
            return check;
        }
    }
    current_     = next;
    currentGcrs_ = nextG;
    prop_        = makePropagator(nextG, gcrs);
    return check;
}

const StateEstimate& EstimateTracker::current() const {
    if (!current_) {
        throw Error("EstimateTracker: no estimate yet");
    }
    return *current_;
}

const NumericalPropagator& EstimateTracker::propagator() const {
    if (!prop_) {
        throw Error("EstimateTracker: no estimate yet");
    }
    return *prop_;
}

StateCovariance EstimateTracker::predict(const Time& t) const {
    const StateCovariance g =
        propagateTo(propagator(), currentGcrs_->covariance, t - current().epoch, opt_.accelerationNoisePsd);
    return opt_.frame == Frame::GCRS ? g : rotated(g, gcrsToCelestialMatrix(opt_.frame, t));
}

void EstimateTracker::predictGrid(double start, double step, std::span<StateCovariance> out) const {
    propagator().propagateGridWithCovariance(start, step, currentGcrs_->covariance, out, opt_.accelerationNoisePsd);
    if (opt_.frame == Frame::GCRS || out.empty()) {
        return;
    }
    const double      end = start + step * static_cast<double>(out.size() - 1);
    const FrameSeries series(opt_.frame, current().epoch, std::min(start, end), std::max(start, end));
    for (std::size_t k = 0; k < out.size(); ++k) {
        out[k] = rotated(out[k], series.at(start + step * static_cast<double>(k)));
    }
}

std::vector<StateCovariance> EstimateTracker::predictGrid(double start, double step, std::size_t count) const {
    std::vector<StateCovariance> out(count);
    predictGrid(start, step, out);
    return out;
}

} // namespace astro
