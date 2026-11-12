/**********************************************************************
 * close_approach.h
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

#include <functional>
#include <iosfwd>
#include <limits>
#include <vector>

#include <boost/describe/enum.hpp>

#include "astro/core/vec.h"
#include "astro/orbit/elements.h"
#include "astro/orbit/secular.h"
#include "astro/time/time.h"

/// \file
/// Closest approach of two objects within a time window, each given by orbital elements, a state
/// vector or any function of time.
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

/// The motion of one object, for closestApproach and closeApproaches.
///
/// Built from a state vector or an element set at an epoch, with a propagation model, or from any
/// function of time (e.g. a NumericalPropagator or an Sgp4). The constructors are not explicit, so
/// `{epoch, elements}` and `{epoch, state}` can be passed directly. Both objects of a search must be
/// in the same inertial frame (for the secular models, one whose z axis is the Earth's pole: TOD,
/// TEME or CIRS). Instances are immutable and safe to share between threads, provided a wrapped
/// function is.
class Trajectory {
public:
    /// A state (m, m/s) at any time.
    using StateFunction = std::function<StateVector(const Time& t)>;

    /// Propagate `state` (m, m/s) from `epoch` with `model`. zp.mu is the gravitational parameter
    /// of every model; the secular models also use zp.re, J2 and J4, and take the osculating
    /// elements of `state` as mean elements (see propagateSecular). Throws astro::Error for a
    /// secular model and a non-elliptic orbit.
    ///
    /// The secular models use SecularPropagatorCorrected, whose velocity is the time derivative of
    /// the position: the search needs the true rate of change of the distance.
    Trajectory(const Time&            epoch,
               const StateVector&     state,
               TrajectoryModel        model = TrajectoryModel::TwoBody,
               const ZonalParameters& zp    = {});
    /// Propagate Keplerian elements from `epoch`: osculating for TwoBody, mean for the secular
    /// models. Throws astro::Error as elementsToState or SecularPropagator do.
    Trajectory(const Time&              epoch,
               const KeplerianElements& elements,
               TrajectoryModel          model = TrajectoryModel::TwoBody,
               const ZonalParameters&   zp    = {});
    /// Propagate equinoctial elements from `epoch` (converted with equinoctialToKeplerian).
    Trajectory(const Time&                epoch,
               const EquinoctialElements& elements,
               TrajectoryModel            model = TrajectoryModel::TwoBody,
               const ZonalParameters&     zp    = {});
    /// Any trajectory, e.g. `[&prop](const Time& t) { return prop.propagate(t); }`. Throws
    /// astro::Error if `stateAt` is empty.
    Trajectory(StateFunction stateAt);

    /// Position (m) and velocity (m/s) at time t.
    [[nodiscard]] StateVector state(const Time& t) const;

    /// The trajectory as a function of seconds after `origin`. The search uses this so that an
    /// element set or state costs one propagation per sample, with no Time arithmetic.
    [[nodiscard]] std::function<StateVector(double)> secondsAfter(const Time& origin) const;

private:
    Time                               epoch_;
    std::function<StateVector(double)> fromEpoch_; // seconds after epoch_ -> state (analytic)
    StateFunction                      stateAt_;   // general trajectories
};

/// Settings for closestApproach and closeApproaches.
struct CloseApproachOptions {
    /// Sampling step, as a fraction of the shorter of the two objects' time scales |r| / |v| at
    /// the current sample (|r| / |v| is a circular orbit's period over 2 pi). The default, 0.05,
    /// samples a circular orbit ~125 times per revolution (45 s in LEO, 11 min in GEO): far finer
    /// than the spacing of the distance's minima and maxima, which for orbits about the same body
    /// are typically a large fraction of a revolution apart. A minimum and maximum closer together
    /// than one step can be missed; 0.1 still found every minimum in a test of 300 random orbit
    /// pairs, 0.2 missed 2 of 1245. Steps are at least 1 ms.
    double stepFraction{0.05};
    /// Time of closest approach tolerance (s). The miss distance is stationary there, so its error
    /// is second order: ~|dv|^2 dt^2 / (2 d), micrometres or less at the default.
    double timeTolerance{1e-6};
};

/// One approach of the two objects.
struct CloseApproach {
    Time        time;               ///< time of the approach (in the time scale of the window start)
    double      distance{0.0};      ///< miss distance |r2 - r1| (m)
    double      relativeSpeed{0.0}; ///< |v2 - v1| (m/s)
    StateVector first;              ///< state of the first object at `time`
    StateVector second;             ///< state of the second object at `time`
    /// The minimum is at the start or the end of the window, where the objects are still closing
    /// (end) or already separating (start), rather than a true closest approach of zero range rate.
    bool atBoundary{false};

    /// Position of the second object relative to the first, r2 - r1 (m).
    [[nodiscard]] Vec3 relativePosition() const noexcept { return second.r - first.r; }
    /// Velocity of the second object relative to the first, v2 - v1 (m/s).
    [[nodiscard]] Vec3 relativeVelocity() const noexcept { return second.v - first.v; }
};

/// The smallest distance between two objects over [start, end].
///
/// The range rate g(t) = (r2 - r1) . (v2 - v1), which is negative while the objects close and
/// positive while they separate, is sampled with the step of `options`. Each change of sign from
/// negative to positive is a local minimum of the distance; it is refined with Brent's method to
/// options.timeTolerance. The result is the smallest of those minima and the distances at the two
/// ends of the window (flagged CloseApproach::atBoundary). A window of zero length gives the
/// distance at `start`.
///
/// Elements and states, or one of each, can be passed directly:
///
///     closestApproach({epoch1, elements1}, {epoch2, state2}, start, end);
///
/// Throws astro::Error if end is before start, for options that are not positive and finite, or
/// if a trajectory returns a non-finite state.
[[nodiscard]] CloseApproach closestApproach(const Trajectory&           first,
                                            const Trajectory&           second,
                                            const Time&                 start,
                                            const Time&                 end,
                                            const CloseApproachOptions& options = {});

/// Every local minimum of the distance strictly inside [start, end] (zero range rate, so never
/// atBoundary) that is no farther than `maxDistance` (m), in time order: a conjunction list. Same
/// search and errors as closestApproach.
[[nodiscard]] std::vector<CloseApproach> closeApproaches(const Trajectory& first,
                                                         const Trajectory& second,
                                                         const Time&       start,
                                                         const Time&       end,
                                                         double maxDistance = std::numeric_limits<double>::infinity(),
                                                         const CloseApproachOptions& options = {});

} // namespace astro
