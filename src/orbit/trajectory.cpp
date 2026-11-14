/**********************************************************************
 * trajectory.cpp
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

#include "astro/orbit/trajectory.h"

#include <cmath>
#include <memory>
#include <string>
#include <utility>

#include "../core/enum_io.h"
#include "astro/core/error.h"
#include "astro/orbit/secular.h"
#include "astro/orbit/two_body.h"
#include "trajectory_detail.h"

namespace astro {

std::ostream& operator<<(std::ostream& os, TrajectoryModel value) {
    return detail::writeEnum(os, "TrajectoryModel", value);
}

namespace detail {

namespace {

/// "Frame::GCRS", for messages.
std::string frameName(Frame f) {
    return std::string("Frame::").append(enumName(f, "INVALID_VALUE"));
}

} // namespace

void requireFrame(const Trajectory& t, Frame needed, std::string_view function, std::string_view what) {
    if (t.frame() != needed) {
        throw Error(std::string(function) + ": " + std::string(what) + " is in " + frameName(t.frame())
                    + "; it must be in " + frameName(needed));
    }
}

void requireSameFrame(const Trajectory& a,
                      const Trajectory& b,
                      std::string_view  function,
                      std::string_view  whatA,
                      std::string_view  whatB) {
    if (a.frame() != b.frame()) {
        throw Error(std::string(function) + ": " + std::string(whatA) + " is in " + frameName(a.frame()) + " and "
                    + std::string(whatB) + " in " + frameName(b.frame()) + "; both must be in the same frame");
    }
}

void requireInertial(const Trajectory& t, std::string_view function, std::string_view what) {
    if (isEarthFixed(t.frame())) {
        throw Error(std::string(function) + ": " + std::string(what) + " is in " + frameName(t.frame())
                    + ", which rotates with the Earth; it must be in an inertial frame (not Frame::TIRS, "
                      "Frame::PEF or Frame::ITRS)");
    }
}

} // namespace detail

namespace {

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

/// `frame`, after checking that it is not Earth-fixed: the two-body and secular models are
/// inertial. `type` names the class for the message.
Frame inertialFrame(Frame frame, std::string_view type) {
    if (isEarthFixed(frame)) {
        throw Error(std::string(type) + ": the two-body and secular models are inertial, and Frame::"
                    + std::string(detail::enumName(frame, "INVALID_VALUE"))
                    + " rotates with the Earth; give the state in an inertial frame");
    }
    return frame;
}

} // namespace

Trajectory::Trajectory(const Time&        epoch,
                       const StateVector& state,
                       Frame              frame,
                       TrajectoryModel    model,
                       const CentralBody& body)
    : frame_(inertialFrame(frame, "Trajectory"))
    , epoch_(epoch) {
    if (model == TrajectoryModel::TwoBody) {
        fromEpoch_ = twoBody(state, body.mu);
    } else {
        fromEpoch_ = secular(std::make_shared<const SecularPropagatorCorrected>(state, body, secularOrder(model)));
    }
}

Trajectory::Trajectory(const Time&              epoch,
                       const KeplerianElements& elements,
                       Frame                    frame,
                       TrajectoryModel          model,
                       const CentralBody&       body)
    : frame_(inertialFrame(frame, "Trajectory"))
    , epoch_(epoch) {
    if (model == TrajectoryModel::TwoBody) {
        fromEpoch_ = twoBody(elementsToState(elements, body.mu), body.mu);
    } else {
        fromEpoch_ = secular(std::make_shared<const SecularPropagatorCorrected>(elements, body, secularOrder(model)));
    }
}

Trajectory::Trajectory(const Time&                epoch,
                       const EquinoctialElements& elements,
                       Frame                      frame,
                       TrajectoryModel            model,
                       const CentralBody&         body)
    : Trajectory(epoch, equinoctialToKeplerian(elements), frame, model, body) {}

Trajectory::Trajectory(Frame frame, StateFunction stateAt)
    : frame_(frame)
    , stateAt_(std::move(stateAt)) {
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

UncertainTrajectory::UncertainTrajectory(const Time&        epoch,
                                         const StateVector& state,
                                         const Mat6&        covariance,
                                         Frame              frame,
                                         TrajectoryModel    model,
                                         const CentralBody& body)
    : trajectory_(epoch, state, frame, model, body) {
    for (const auto& row : covariance.m) {
        for (double x : row) {
            if (!std::isfinite(x)) {
                throw Error("UncertainTrajectory: covariance is not finite");
            }
        }
    }
    if (model == TrajectoryModel::TwoBody) {
        covarianceAt_ = [epoch, state, covariance, mu = body.mu](const Time& t) {
            return propagateKeplerWithCovariance(state, t - epoch, covariance, mu);
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

Trajectory sgp4Trajectory(const Sgp4& sgp4, const EopProvider& eop, Frame frame) {
    if (frame == Frame::TEME) {
        return {Frame::TEME, [&sgp4](const Time& t) { return sgp4.propagate(t); }};
    }
    // A FrameTransformer per call: the rotations at time t, IAU 2000B nutation for speed.
    return {frame, [&sgp4, &eop, frame](const Time& t) {
                const FrameTransformer ft(t, eop, {.model = PrecessionNutationModel::IAU2006_2000B});
                return ft.transform(sgp4.propagate(t), Frame::TEME, frame);
            }};
}

Trajectory bodyTrajectory(const Ephemeris& eph, Body target, Body center) {
    return {Frame::GCRS, [&eph, target, center](const Time& t) { return eph.state(target, center, t); }};
}

} // namespace astro
