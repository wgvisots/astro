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

#include <concepts>
#include <functional>
#include <iosfwd>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

#include <boost/describe/enum.hpp>

#include "astro/core/mat6.h"
#include "astro/core/vec.h"
#include "astro/orbit/elements.h"
#include "astro/orbit/secular.h"
#include "astro/orbit/state_covariance.h"
#include "astro/time/time.h"

/// \file
/// Closest approach of two objects within a time window, each given by orbital elements, a state
/// vector or any function of time, and the probability that they collide; assessConjunction does
/// both from states and covariances.
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
/// `{epoch, elements}` and `{epoch, state}` can be passed directly. Both objects of a search must
/// be in the same inertial frame (for the secular models, one whose z axis is the Earth's pole:
/// TOD, TEME or CIRS). Instances are immutable and safe to share between threads, provided a
/// wrapped function is.
class Trajectory {
public:
    /// A state (m, m/s) at any time.
    using StateFunction = std::function<StateVector(const Time& t)>;

    /// Propagate `state` (m, m/s) from `epoch` with `model`. body.mu is the gravitational parameter
    /// of every model; the secular models also use body.re, J2 and J4 (J3 and higher do not enter),
    /// and take the osculating elements of `state` as mean elements (see propagateSecular). Throws
    /// astro::Error for a secular model and a non-elliptic orbit.
    ///
    /// The secular models use SecularPropagatorCorrected, whose velocity is the time derivative of
    /// the position: the search needs the true rate of change of the distance.
    Trajectory(const Time&        epoch,
               const StateVector& state,
               TrajectoryModel    model = TrajectoryModel::TwoBody,
               const CentralBody& body  = CentralBody::earth());
    /// Propagate Keplerian elements from `epoch`: osculating for TwoBody, mean for the secular
    /// models. Throws astro::Error as elementsToState or SecularPropagator do.
    Trajectory(const Time&              epoch,
               const KeplerianElements& elements,
               TrajectoryModel          model = TrajectoryModel::TwoBody,
               const CentralBody&       body  = CentralBody::earth());
    /// Propagate equinoctial elements from `epoch` (converted with equinoctialToKeplerian, which
    /// throws astro::Error unless they describe an ellipse).
    Trajectory(const Time&                epoch,
               const EquinoctialElements& elements,
               TrajectoryModel            model = TrajectoryModel::TwoBody,
               const CentralBody&         body  = CentralBody::earth());
    /// Any trajectory. Throws astro::Error if `stateAt` is empty.
    Trajectory(StateFunction stateAt);
    /// Any callable taking a Time and returning a StateVector, so that a lambda can be passed
    /// wherever a Trajectory is expected: `[&prop](const Time& t) { return prop.propagate(t); }`.
    /// Throws astro::Error if `stateAt` converts to an empty function (a null function pointer,
    /// for one).
    template <class F>
        requires(!std::same_as<std::remove_cvref_t<F>, Trajectory>
                 && !std::same_as<std::remove_cvref_t<F>, StateFunction>
                 && std::is_invocable_r_v<StateVector, F&, const Time&>)
    Trajectory(F&& stateAt)
        : Trajectory(StateFunction(std::forward<F>(stateAt))) {}

    /// Position (m) and velocity (m/s) at time t. Throws astro::Error as propagateKepler does
    /// (TwoBody model), and if t cannot be subtracted from the epoch (one of the two in UT1, or a
    /// UTC time before 1960; see Time). An exception thrown by a wrapped function passes through.
    [[nodiscard]] StateVector state(const Time& t) const;

    /// The trajectory as a function of seconds after `origin`. The search uses this so that an
    /// element set or state costs one propagation per sample, with no Time arithmetic. Throws
    /// astro::Error if `origin` cannot be subtracted from the epoch (see Time). The function
    /// returned throws as state() does.
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
    /// pairs, 0.2 missed 2 of 1245 (astro_study_conjunctions). Steps are at least 1 ms.
    double stepFraction{0.05};
    /// Time of closest approach tolerance (s). The miss distance d is stationary there, so its
    /// error is second order, ~|dv|^2 dt^2 / (2 d): below a micrometre at the default for misses
    /// over ~100 m at 15 km/s. It can never exceed |dv| dt (15 mm at 15 km/s), which is what a
    /// near-zero miss gets.
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
/// negative to positive is a local minimum of the distance (as in Alfano & Negron 1993); it is
/// refined with Brent's method (Brent 1973) to options.timeTolerance. The result is the smallest of
/// those minima and the distances at the two ends of the window (flagged
/// CloseApproach::atBoundary). A window of zero length gives the distance at `start`.
///
/// Elements and states, or one of each, can be passed directly:
///
///     closestApproach({epoch1, elements1}, {epoch2, state2}, start, end);
///
/// Throws astro::Error if end is before start, for options that are not positive and finite, or
/// if a trajectory returns a non-finite state. Also throws as Trajectory::state does, for either
/// trajectory, and if `end - start` does (see Time).
[[nodiscard]] CloseApproach closestApproach(const Trajectory&           first,
                                            const Trajectory&           second,
                                            const Time&                 start,
                                            const Time&                 end,
                                            const CloseApproachOptions& options = {});

/// Every local minimum of the distance in (start, end] that is no farther than `maxDistance` (m),
/// in time order: a conjunction list. Each is a point of zero range rate reached while closing, so
/// none is flagged atBoundary; one that falls exactly on `end` is included, one exactly on `start`
/// is not. Objects that keep a constant distance (the same orbit, a fixed phase apart) have no
/// minimum: a change of sign of the range rate within its rounding error is ignored. Same search
/// as closestApproach. Throws astro::Error as closestApproach does, and if `maxDistance` is NaN.
[[nodiscard]] std::vector<CloseApproach> closeApproaches(const Trajectory& first,
                                                         const Trajectory& second,
                                                         const Time&       start,
                                                         const Time&       end,
                                                         double maxDistance = std::numeric_limits<double>::infinity(),
                                                         const CloseApproachOptions& options = {});

/// Probability that a 2-D Gaussian point with mean (mx, my) and covariance [[xx, xy], [xy, yy]]
/// (m^2) lies within `radius` (m) of the origin.
///
/// The covariance is rotated to its principal axes (the minor variance as determinant / major
/// variance, which stays accurate however elongated the ellipse). Along the minor axis the integral
/// over the circle's chord is exact (error functions, in complementary form for small
/// probabilities); along the major axis it is integrated numerically with u = radius sin t, which
/// removes the square-root end points, by adaptive Gauss-Kronrod (7-15; Piessens et al. 1983,
/// QUADPACK): the 1-D form of Alfano (2005a). The range is cut where the chord meets the minor
/// axis's offset, so that a covariance far narrower than the circle is not stepped over.
///
/// Relative accuracy ~1e-12, and ~1e-11 for a near-singular covariance whose axis grazes the
/// circle; 7-20 us, at most ~1 ms. A singular covariance (one or both standard deviations zero) is
/// handled exactly. A mean more than 15 standard deviations outside the circle along the major axis
/// gives 0 (the probability is below ~1e-49). Throws astro::Error for a radius that is not positive
/// and finite, a non-finite argument, or a covariance that is not positive semi-definite.
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
/// & Estes 1992; Akella & Alfriend 2000; Alfano 2005a; links in docs/references.md).
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
/// encounters (m/s) the relative motion curves during the encounter and the result is not reliable
/// (Patera 2003 treats curved relative motion). Use the approach's own states (a true closest
/// approach, not atBoundary) for the encounter plane. Throws astro::Error for a radius that is not
/// positive and finite, a non-finite state, a zero relative velocity, or a covariance that is not
/// finite and positive semi-definite in the plane.
[[nodiscard]] CollisionProbability collisionProbability(const CloseApproach& approach,
                                                        const Mat3&          covariance1,
                                                        const Mat3&          covariance2,
                                                        double               hardBodyRadius);

/// As above, from 6x6 position/velocity covariances (only the position blocks are used), e.g. from
/// NumericalPropagator::propagateGridWithCovariance at the time of the approach. Throws
/// astro::Error as the overload above does.
[[nodiscard]] CollisionProbability collisionProbability(const CloseApproach& approach,
                                                        const Mat6&          covariance1,
                                                        const Mat6&          covariance2,
                                                        double               hardBodyRadius);

/// The motion and position/velocity uncertainty of one object, for assessConjunction and
/// assessConjunctions.
///
/// Built from a state and its covariance at an epoch, with a propagation model, or from any
/// Trajectory with a function giving the state and covariance at any time (e.g. from a
/// NumericalPropagator or EstimateTracker::predict). Instances are immutable and safe to share
/// between threads, provided wrapped functions are.
class UncertainTrajectory {
public:
    /// A state (m, m/s) and its covariance (m^2, m^2/s, m^2/s^2) at any time.
    using StateCovarianceFunction = std::function<StateCovariance(const Time& t)>;

    /// Propagate `state` (m, m/s) and its covariance (m^2, m^2/s, m^2/s^2, in the frame of the
    /// state) from `epoch` with `model`, as Trajectory does: TwoBody with propagateKeplerWithCovariance
    /// (analytic state transition matrix), the secular models with
    /// SecularPropagatorCorrected::propagateWithCovariance (differenced). Neither adds process noise.
    /// Throws astro::Error for a non-finite covariance, or as Trajectory does.
    UncertainTrajectory(const Time&        epoch,
                        const StateVector& state,
                        const Mat6&        covariance,
                        TrajectoryModel    model = TrajectoryModel::TwoBody,
                        const CentralBody& body  = CentralBody::earth());
    /// Any trajectory, with a function giving its state and covariance at any time. The search
    /// samples `trajectory`; `stateCovarianceAt` is called once per approach, for its covariance.
    /// Throws astro::Error if `stateCovarianceAt` is empty.
    UncertainTrajectory(Trajectory trajectory, StateCovarianceFunction stateCovarianceAt);

    /// The trajectory the search samples.
    [[nodiscard]] const Trajectory& trajectory() const noexcept { return trajectory_; }
    /// State and covariance at time t. Throws astro::Error as propagateKeplerWithCovariance
    /// (TwoBody model) or SecularPropagator::propagateWithCovariance (secular models) does, and if
    /// t cannot be subtracted from the epoch (see Time). An exception thrown by a wrapped function
    /// passes through.
    [[nodiscard]] StateCovariance stateCovariance(const Time& t) const;

private:
    Trajectory              trajectory_;
    StateCovarianceFunction covarianceAt_;
};

/// One approach of two objects with its collision probability.
struct ConjunctionAssessment {
    CloseApproach        approach;         ///< the approach (states at approach.time)
    Mat6                 firstCovariance;  ///< covariance of the first object at approach.time
    Mat6                 secondCovariance; ///< covariance of the second object at approach.time
    CollisionProbability collision;        ///< collisionProbability of the approach
};

/// The closest approach of two objects over [start, end] (closestApproach), with both covariances
/// propagated to it and its collision probability (collisionProbability) for spheres whose radii
/// add up to `hardBodyRadius` (m):
///
///     const ConjunctionAssessment c = assessConjunction({epoch1, state1, cov1}, {epoch2, state2, cov2},
///                                                       start, end, 10.0);
///
/// The closest approach is not always the most probable collision; assessConjunctions gives every
/// approach. If it is at an end of the window (approach.atBoundary) it is not a true closest
/// approach, and the probability is that of the straight-line encounter through it. Throws
/// astro::Error as closestApproach, UncertainTrajectory::stateCovariance and collisionProbability
/// do.
[[nodiscard]] ConjunctionAssessment assessConjunction(const UncertainTrajectory&  first,
                                                      const UncertainTrajectory&  second,
                                                      const Time&                 start,
                                                      const Time&                 end,
                                                      double                      hardBodyRadius,
                                                      const CloseApproachOptions& options = {});

/// Every local minimum of the distance strictly inside [start, end] no farther than `maxDistance`
/// (m), as closeApproaches, each with its covariances and collision probability, in time order.
/// Throws astro::Error as closeApproaches, UncertainTrajectory::stateCovariance and
/// collisionProbability do.
[[nodiscard]] std::vector<ConjunctionAssessment>
assessConjunctions(const UncertainTrajectory&  first,
                   const UncertainTrajectory&  second,
                   const Time&                 start,
                   const Time&                 end,
                   double                      hardBodyRadius,
                   double                      maxDistance = std::numeric_limits<double>::infinity(),
                   const CloseApproachOptions& options     = {});

} // namespace astro
