/**********************************************************************
 * trajectory.h
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

#pragma once

#include <concepts>
#include <functional>
#include <iosfwd>
#include <type_traits>
#include <utility>

#include <boost/describe/enum.hpp>

#include "astro/core/mat6.h"
#include "astro/core/vec.h"
#include "astro/ephemeris/ephemeris.h"
#include "astro/frames/eop.h"
#include "astro/frames/transform.h"
#include "astro/orbit/central_body.h"
#include "astro/orbit/elements.h"
#include "astro/orbit/secular.h"
#include "astro/orbit/sgp4.h"
#include "astro/orbit/state_covariance.h"
#include "astro/time/time.h"

/// \file
/// A state as a function of time that carries its reference frame (Trajectory), the same with its
/// covariance (UncertainTrajectory), and helpers that make the usual ones from SGP4 and from an
/// ephemeris. The search functions (findPasses, groundTrack, findEclipses, closestApproach, ...)
/// take a Trajectory and check its frame.
namespace astro {

/// How Trajectory propagates an element set or a state vector.
enum class TrajectoryModel {
    TwoBody,     ///< point mass (propagateKepler): any eccentricity
    SecularJ2,   ///< SecularPropagatorCorrected, first-order J2 (elliptic only; z axis = Earth's pole)
    SecularJ2J4, ///< SecularPropagatorCorrected, J2 + J2^2 + J4 (elliptic only; z axis = Earth's pole)
};

BOOST_DESCRIBE_ENUM(TrajectoryModel, TwoBody, SecularJ2, SecularJ2J4)

/// Writes "TrajectoryModel::<enumerator>", or "TrajectoryModel::INVALID_VALUE" for any other value.
std::ostream& operator<<(std::ostream& os, TrajectoryModel value);

/// The motion of one object: its state at any time, and the frame that state is in.
///
/// Built from a state vector or an element set at an epoch, with a propagation model, or from any
/// function of time (e.g. a NumericalPropagator or an Sgp4; see also sgp4Trajectory and
/// bodyTrajectory). Every constructor takes the frame, and the functions that use a trajectory
/// check it: findPasses and groundTrack need Frame::ITRS, closestApproach needs both objects in
/// the same inertial frame, and a wrong frame throws astro::Error instead of giving a different
/// number. The constructors are not explicit, so `{epoch, elements, Frame::GCRS}` and
/// `{Frame::ITRS, function}` can be passed directly. Instances are immutable and safe to share
/// between threads, provided a wrapped function is.
class Trajectory {
public:
    /// A state (m, m/s) at any time.
    using StateFunction = std::function<StateVector(const Time& t)>;

    /// Propagate `state` (m, m/s, in `frame`) from `epoch` with `model`. body.mu is the
    /// gravitational parameter of every model; the secular models also use body.re, J2 and J4 (J3
    /// and higher do not enter), and take the osculating elements of `state` as mean elements (see
    /// propagateSecular). Their zonal field acts about z, so for them `frame` must have its z axis
    /// along the true pole: Frame::TOD, Frame::TEME or Frame::CIRS. Throws astro::Error for an
    /// Earth-fixed `frame` (see isEarthFixed: the models are inertial), for a secular model in any
    /// other frame (convert with poleOfDateMatrix), and for a secular model and a non-elliptic orbit.
    ///
    /// The secular models use SecularPropagatorCorrected, whose velocity is the time derivative of
    /// the position: the search needs the true rate of change of the distance.
    Trajectory(const Time&        epoch,
               const StateVector& state,
               Frame              frame,
               TrajectoryModel    model = TrajectoryModel::TwoBody,
               const CentralBody& body  = CentralBody::earth());
    /// Propagate Keplerian elements (in `frame`) from `epoch`: osculating for TwoBody, mean for the
    /// secular models. Throws astro::Error for the frames the state constructor refuses, and as
    /// elementsToState or SecularPropagator do.
    Trajectory(const Time&              epoch,
               const KeplerianElements& elements,
               Frame                    frame,
               TrajectoryModel          model = TrajectoryModel::TwoBody,
               const CentralBody&       body  = CentralBody::earth());
    /// Propagate equinoctial elements (in `frame`) from `epoch` (converted with
    /// equinoctialToKeplerian, which throws astro::Error unless they describe an ellipse). Throws
    /// astro::Error for an Earth-fixed `frame`.
    Trajectory(const Time&                epoch,
               const EquinoctialElements& elements,
               Frame                      frame,
               TrajectoryModel            model = TrajectoryModel::TwoBody,
               const CentralBody&         body  = CentralBody::earth());
    /// Any trajectory: `stateAt` returns the state in `frame`, which may be any frame, Earth-fixed
    /// included. Throws astro::Error if `stateAt` is empty.
    Trajectory(Frame frame, StateFunction stateAt);
    /// Any callable taking a Time and returning a StateVector in `frame`, so that a lambda can be
    /// passed wherever a Trajectory is expected:
    /// `{Frame::GCRS, [&prop](const Time& t) { return prop.propagate(t); }}`. Throws astro::Error
    /// if `stateAt` converts to an empty function (a null function pointer, for one).
    template <class F>
        requires(!std::same_as<std::remove_cvref_t<F>, StateFunction>
                 && std::is_invocable_r_v<StateVector, F&, const Time&>)
    Trajectory(Frame frame, F&& stateAt)
        : Trajectory(frame, StateFunction(std::forward<F>(stateAt))) {}

    /// The frame of the states this trajectory returns.
    [[nodiscard]] Frame frame() const noexcept { return frame_; }

    /// Position (m) and velocity (m/s) at time t, in frame(). Throws astro::Error as
    /// propagateKepler does (TwoBody model), and if t cannot be subtracted from the epoch (one of
    /// the two in UT1, or a UTC time before 1960; see Time). An exception thrown by a wrapped
    /// function passes through.
    [[nodiscard]] StateVector state(const Time& t) const;

    /// The trajectory as a function of seconds after `origin`. The search uses this so that an
    /// element set or state costs one propagation per sample, with no Time arithmetic. Throws
    /// astro::Error if `origin` cannot be subtracted from the epoch (see Time). The function
    /// returned throws as state() does.
    [[nodiscard]] std::function<StateVector(double)> secondsAfter(const Time& origin) const;

private:
    Frame                              frame_;
    Time                               epoch_;
    std::function<StateVector(double)> fromEpoch_; // seconds after epoch_ -> state (analytic)
    StateFunction                      stateAt_;   // general trajectories
};

/// The motion and position/velocity uncertainty of one object, for assessConjunction and
/// assessConjunctions.
///
/// Built from a state and its covariance at an epoch, with a propagation model, or from any
/// Trajectory with a function giving the state and covariance at any time (e.g. from a
/// NumericalPropagator or EstimateTracker::predict). The covariance is in the frame of the
/// trajectory. Instances are immutable and safe to share between threads, provided wrapped
/// functions are.
class UncertainTrajectory {
public:
    /// A state (m, m/s) and its covariance (m^2, m^2/s, m^2/s^2) at any time.
    using StateCovarianceFunction = std::function<StateCovariance(const Time& t)>;

    /// Propagate `state` (m, m/s) and its covariance (m^2, m^2/s, m^2/s^2), both in `frame`, from
    /// `epoch` with `model`, as Trajectory does: TwoBody with propagateKeplerWithCovariance
    /// (analytic state transition matrix), the secular models with
    /// SecularPropagatorCorrected::propagateWithCovariance (differenced). Neither adds process noise.
    /// Throws astro::Error for a non-finite covariance, or as Trajectory does (an Earth-fixed
    /// `frame`, or a secular model in a frame whose z is not the true pole).
    UncertainTrajectory(const Time&        epoch,
                        const StateVector& state,
                        const Mat6&        covariance,
                        Frame              frame,
                        TrajectoryModel    model = TrajectoryModel::TwoBody,
                        const CentralBody& body  = CentralBody::earth());
    /// Any trajectory, with a function giving its state and covariance (in the trajectory's frame)
    /// at any time. The search samples `trajectory`; `stateCovarianceAt` is called once per
    /// approach, for its covariance. Throws astro::Error if `stateCovarianceAt` is empty.
    UncertainTrajectory(Trajectory trajectory, StateCovarianceFunction stateCovarianceAt);

    /// The trajectory the search samples.
    [[nodiscard]] const Trajectory& trajectory() const noexcept { return trajectory_; }
    /// The frame of the trajectory and of the covariance.
    [[nodiscard]] Frame frame() const noexcept { return trajectory_.frame(); }
    /// State and covariance at time t. Throws astro::Error as propagateKeplerWithCovariance
    /// (TwoBody model) or SecularPropagator::propagateWithCovariance (secular models) does, and if
    /// t cannot be subtracted from the epoch (see Time). An exception thrown by a wrapped function
    /// passes through.
    [[nodiscard]] StateCovariance stateCovariance(const Time& t) const;

private:
    Trajectory              trajectory_;
    StateCovarianceFunction covarianceAt_;
};

/// The trajectory of an SGP4 satellite in `frame`: Sgp4::propagate (TEME) converted to `frame` with
/// FrameTransformer and the given EOP (the IAU 2000B nutation model is used for speed). In
/// Frame::TEME it is Sgp4::propagate itself, and `eop` is not used.
///
/// The trajectory keeps references to `sgp4` and `eop`: both must outlive it (and its copies).
/// Its state() throws astro::Error for a time at which `eop` has no data or which is before 1960
/// (see FrameTransformer), and as Sgp4::propagate does.
[[nodiscard]] Trajectory sgp4Trajectory(const Sgp4& sgp4, const EopProvider& eop, Frame frame);

/// The geocentric trajectory of `target` from an ephemeris (Ephemeris::state relative to
/// Body::Earth), in Frame::GCRS: geocentric ephemeris axes are the ICRF's, which are the GCRS's.
/// (It is geocentric only because Frame has no barycentric or heliocentric frame to label any
/// other centre with; for those, call Ephemeris::state.)
/// For the Sun in findEclipses: `bodyTrajectory(eph, Body::Sun)`.
///
/// The trajectory keeps a reference to `eph`: it must outlive it (and its copies). Its state()
/// throws astro::Error as Ephemeris::state does.
[[nodiscard]] Trajectory bodyTrajectory(const Ephemeris& eph, Body target);

} // namespace astro
