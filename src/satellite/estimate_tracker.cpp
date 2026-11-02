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

#include <cmath>
#include <memory>

#include "../core/enum_io.h"
#include "astro/core/error.h"
#include "astro/satellite/orbital_frames.h"

namespace astro {

std::ostream& operator<<(std::ostream& os, ConsistencyTest value) {
    return detail::writeEnum(os, "ConsistencyTest", value);
}

double chiSquareSurvival3(double x) noexcept {
    if (!(x > 0.0)) {
        return 1.0;
    }
    return std::erfc(std::sqrt(0.5 * x)) + std::sqrt(2.0 * x / constants::pi) * std::exp(-0.5 * x);
}

double chiSquareSurvival6(double x) noexcept {
    if (!(x > 0.0)) {
        return 1.0;
    }
    const double h = 0.5 * x;
    return std::exp(-h) * (1.0 + h + 0.5 * h * h);
}

NumericalPropagator makePropagator(const StateEstimate& e, const TrackerOptions& options) {
    return NumericalPropagator(e.epoch,
                               e.state,
                               {std::make_shared<GravityField>(options.J, options.mu, options.re)},
                               options.integrator);
}

namespace {

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
    c.valid                        = d6.has_value() && d3.has_value();
    if (c.valid) {
        c.mahalanobis2         = *d6;
        c.positionMahalanobis2 = *d3;
        c.pValue               = chiSquareSurvival6(*d6);
        c.positionPValue       = chiSquareSurvival3(*d3);
        c.consistent = (opt.test == ConsistencyTest::Position ? c.positionPValue : c.pValue) >= opt.alarmProbability;
    } else {
        c.mahalanobis2 = c.positionMahalanobis2 = c.pValue = c.positionPValue = std::nan("");
    }
    return c;
}

StateCovariance propagateTo(const NumericalPropagator& prop, const Mat6& P0, double dt, double q) {
    std::vector<StateCovariance> out(1);
    prop.propagateGridWithCovariance(dt, 1.0, P0, out, q); // a single output: integrates exactly to dt
    return out[0];
}

} // namespace

ConsistencyCheck
checkConsistency(const StateEstimate& prior, const StateEstimate& next, const TrackerOptions& options) {
    const double     dt   = next.epoch - prior.epoch;
    const auto       prop = makePropagator(prior, options);
    ConsistencyCheck c =
        compare(propagateTo(prop, prior.covariance, dt, options.accelerationNoisePsd), next, dt, options);
    c.replaced = dt >= 0.0;
    return c;
}

EstimateTracker::EstimateTracker(TrackerOptions options)
    : opt_(std::move(options)) {}

std::optional<ConsistencyCheck> EstimateTracker::update(const StateEstimate& next) {
    std::optional<ConsistencyCheck> check;
    if (current_) {
        const double dt = next.epoch - current_->epoch;
        check = compare(propagateTo(*prop_, current_->covariance, dt, opt_.accelerationNoisePsd), next, dt, opt_);
        check->replaced = dt >= 0.0;
        if (!check->replaced) {
            return check;
        }
    }
    current_ = next;
    prop_    = makePropagator(next, opt_);
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
    return propagateTo(propagator(), current().covariance, t - current().epoch, opt_.accelerationNoisePsd);
}

void EstimateTracker::predictGrid(double start, double step, std::span<StateCovariance> out) const {
    propagator().propagateGridWithCovariance(start, step, current().covariance, out, opt_.accelerationNoisePsd);
}

std::vector<StateCovariance> EstimateTracker::predictGrid(double start, double step, std::size_t count) const {
    std::vector<StateCovariance> out(count);
    predictGrid(start, step, out);
    return out;
}

} // namespace astro
