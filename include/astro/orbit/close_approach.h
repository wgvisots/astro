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

#include "astro/core/mat6.h"
#include "astro/core/vec.h"
#include "astro/orbit/elements.h"
#include "astro/orbit/secular.h"
#include "astro/time/time.h"

/// \file
/// Closest approach of two objects within a time window, each given by orbital elements, a state
/// vector or any function of time, and the probability that they collide.
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

/// Probability that a 2-D Gaussian point with mean (mx, my) and covariance [[xx, xy], [xy, yy]]
/// (m^2) lies within `radius` (m) of the origin.
///
/// The covariance is rotated to its principal axes. Along the minor axis the integral over the
/// circle's chord is exact (error functions, in complementary form for small probabilities); along
/// the major axis it is integrated numerically with u = radius sin t, which removes the square-root
/// end points, by adaptive Gauss-Kronrod (7-15) to ~1e-12 relative. A singular covariance (one or both standard
/// deviations zero) is handled exactly. Probabilities below ~1e-45 come out as 0. Throws
/// astro::Error for a radius that is not positive and finite, or a covariance that is not positive
/// semi-definite.
[[nodiscard]] double probabilityInCircle(double mx, double my, double xx, double xy, double yy, double radius);

/// The collision probability of one approach, and the encounter-plane quantities behind it.
struct CollisionProbability {
    double probability{0.0};  ///< probability that the objects pass within the hard-body radius
    double missDistance{0.0}; ///< miss distance in the encounter plane (m)
    double xx{0.0};           ///< combined position covariance in the encounter plane, xAxis-xAxis (m^2)
    double xy{0.0};           ///< combined position covariance in the encounter plane, xAxis-yAxis (m^2)
    double yy{0.0};           ///< combined position covariance in the encounter plane, yAxis-yAxis (m^2)
    double sigmaMajor{0.0};   ///< 1-sigma semi-major axis of that covariance (m)
    double sigmaMinor{0.0};   ///< 1-sigma semi-minor axis of that covariance (m)
    /// Miss distance in standard deviations, sqrt(m^T C^-1 m); infinite for a nonzero miss with a
    /// singular covariance.
    double mahalanobis{0.0};
    Vec3   xAxis; ///< encounter frame x: along the miss vector (any direction perpendicular to zAxis if 0)
    Vec3   yAxis; ///< encounter frame y: zAxis x xAxis
    Vec3   zAxis; ///< encounter frame z: along the relative velocity v2 - v1
};

/// Probability that two objects collide at an approach, by the short-term encounter model (Foster
/// & Estes 1992; Akella & Alfriend 2000; Alfano 2005).
///
/// The objects are taken as spheres whose radii add up to `hardBodyRadius` (m). `covariance1` and
/// `covariance2` are their position covariances (m^2) at approach.time, in the frame of the
/// approach's states, and their errors are independent, so the relative position has covariance
/// C = C1 + C2. During the encounter the relative motion is taken as a straight line at constant
/// velocity, with C constant: the probability is then the integral of the relative position's
/// Gaussian, projected onto the plane perpendicular to the relative velocity, over the
/// hard-body circle (probabilityInCircle).
///
/// That holds when the encounter is short: the objects cross the region of their uncertainty in a
/// small fraction of an orbit, as for relative speeds of km/s. For slow, nearly co-orbital
/// encounters (m/s) the relative motion curves during the encounter and the result is not
/// reliable. Use the approach's own states (a true closest approach, not atBoundary) for the
/// encounter plane. Throws astro::Error for a radius that is not positive and finite, a zero
/// relative velocity, or a covariance that is not finite and positive semi-definite in the plane.
[[nodiscard]] CollisionProbability collisionProbability(const CloseApproach& approach,
                                                        const Mat3&          covariance1,
                                                        const Mat3&          covariance2,
                                                        double               hardBodyRadius);

/// As above, from 6x6 position/velocity covariances (only the position blocks are used), e.g. from
/// NumericalPropagator::propagateGridWithCovariance at the time of the approach.
[[nodiscard]] CollisionProbability collisionProbability(const CloseApproach& approach,
                                                        const Mat6&          covariance1,
                                                        const Mat6&          covariance2,
                                                        double               hardBodyRadius);

} // namespace astro
