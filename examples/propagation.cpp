/**********************************************************************
 * propagation.cpp
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

// Propagation: two-body and Lambert's problem, analytic J2/J4 (secular), numerical J2/J4 on an output grid, a
// 12-component (position, velocity, acceleration, jerk) state with an unmodeled force, and the
// closest approach of two objects with their collision probability, from states and covariances.

#include <iomanip>
#include <iostream>

#include "astro/core/angle.h"
#include "astro/orbit/close_approach.h"
#include "astro/orbit/elements.h"
#include "astro/orbit/numerical.h"
#include "astro/orbit/secular.h"
#include "astro/orbit/two_body.h"
#include "astro/satellite/orbital_frames.h"

using namespace astro;

int main() {
    std::cout << std::fixed << std::setprecision(1);
    const Time        epoch = Time::fromIso8601("2026-03-01T00:00:00Z");
    const StateVector s0    = elementsToState({6778e3, 0.0005, deg2rad(51.6), 0.3, 0.5, 0.7}); // GCRS, m, m/s

    //! [numerical]
    // Numerical J2 + J4, 4 h at 10 s. Say which frame the state is in: the zonal field is applied
    // about the true pole (in GCRS ~0.15 deg off z). relTol 1e-9 keeps LEO within ~0.2 m over 4 h.
    AdaptiveOptions opt;
    opt.relTol                          = 1e-9;
    opt.absTol                          = 1e-6;
    const NumericalPropagator      prop = NumericalPropagator::j2j4(epoch, s0, Frame::GCRS, opt);
    const std::vector<StateVector> grid = prop.propagateGrid(0.0, 10.0, 1441); // t = 0, 10, ..., 14400 s
    std::cout << "numerical, 4 h: " << grid.back().r << " m\n";
    //! [numerical]

    //! [two_body]
    const StateVector kepler = propagateKepler(s0, 14400.0); // universal variables, point mass only
    std::cout << "two-body differs from J2/J4 by " << norm(kepler.r - grid.back().r) / 1e3 << " km\n";
    //! [two_body]

    //! [lambert]
    // Lambert's problem: the orbit that joins two positions in a given time. These two are 40 minutes
    // apart on the two-body orbit above, so the solution gives back its velocity.
    const StateVector     later    = propagateKepler(s0, 2400.0);
    const LambertSolution transfer = lambert(s0.r, later.r, 2400.0);
    std::cout << "Lambert recovers the velocity to " << norm(transfer.v1 - s0.v) * 1e6 << " um/s\n";
    // With 5 hours to spare, other orbits make complete revolutions on the way: two for each number.
    const std::vector<LambertSolution> all = lambertIzzo(s0.r, later.r, 18000.0, constants::gmEarth, true, 3);
    std::cout << all.size() << " orbits join the same points in 5 h with up to 3 revolutions\n";
    //! [lambert]

    //! [secular]
    // Analytic J2/J4 mean-element drift. Fast (~40 ns per output) but the state's osculating
    // elements are taken as mean elements, so expect kilometres of error: see the propagation guide.
    const SecularPropagator        secular(s0); // z must be the pole: fine for this illustration only
    const std::vector<StateVector> sec = secular.propagateGrid(0.0, 10.0, 1441);
    std::cout << "secular differs from numerical by " << norm(sec.back().r - grid.back().r) / 1e3 << " km\n";
    //! [secular]

    //! [state12]
    // A 12-component state: the input acceleration and jerk differ from gravity by a small thrust,
    // which UnmodeledForce keeps acting (da0 + dj0 t) during the propagation.
    StateVector12 s12     = prop.withDerivatives(0.0, s0); // r, v and the model's a and j
    s12.a                += normalized(s0.v) * 1e-3;       // plus 1 mm/s^2 along-track
    const auto thrusting  = NumericalPropagator::j2j4(epoch, s12, Frame::GCRS, HigherOrderMode::UnmodeledForce, opt);
    std::cout << std::setprecision(6) << norm(thrusting.accelerationResidual()) << " m/s^2 unmodeled; after 1 h "
              << std::setprecision(1) << norm(thrusting.propagate(3600.0).r - prop.propagate(3600.0).r) / 1e3
              << " km from the unthrusted orbit\n";
    std::cout << "mode: " << HigherOrderMode::UnmodeledForce << '\n'; // enums print as Enum::Value
    //! [state12]

    //! [close_approach]
    // Closest approach within a day of a state vector and an element set (either may be either, or
    // any function of time, e.g. [&prop](const Time& t) { return prop.propagate(t); }). Both must be
    // in the same inertial frame; the default model is two-body.
    const KeplerianElements debris{6780e3, 0.001, deg2rad(97.6), 0.25, 0.0, 1.22};
    const Time              end = epoch + 86400.0;
    const CloseApproach     ca  = closestApproach({epoch, s0}, {epoch, debris}, epoch, end);
    std::cout << "closest approach " << ca.distance / 1e3 << " km at " << ca.time.iso8601(0) << ", "
              << ca.relativeSpeed / 1e3 << " km/s\n";
    // Every approach within 10 km: a conjunction list.
    for (const CloseApproach& c : closeApproaches({epoch, s0}, {epoch, debris}, epoch, end, 10e3)) {
        std::cout << "  " << c.time.iso8601(0) << "  " << c.distance / 1e3 << " km\n";
    }
    //! [close_approach]

    //! [collision_probability]
    // Probability of collision at that approach, for a 5 m and a 1 m object. Position covariances at
    // the time of the approach: 400 m radial, 4 km in-track, 400 m cross-track (1 sigma) in each
    // object's RIC frame, rotated into the inertial frame. (From an epoch covariance,
    // NumericalPropagator::propagateWithCovariance gives them at any time.)
    // clang-format off
    const Mat3 ric{400.0 * 400.0,             0.0,           0.0,
                             0.0, 4000.0 * 4000.0,           0.0,
                             0.0,             0.0, 400.0 * 400.0};
    // clang-format on
    const Mat3                 m1 = ricMatrix(ca.first), m2 = ricMatrix(ca.second);
    const CollisionProbability pc = collisionProbability(ca, m1.transpose() * ric * m1, m2.transpose() * ric * m2, 6.0);
    std::cout << "collision probability " << std::scientific << std::setprecision(2) << pc.probability << std::fixed
              << std::setprecision(1) << ", miss " << pc.mahalanobis << " sigma\n";
    //! [collision_probability]

    //! [conjunction_assessment]
    // The same from each object's state and covariance at the epoch, in one call: the covariances
    // are propagated to every approach within 10 km with the trajectory's model (two-body here;
    // TrajectoryModel::SecularJ2J4 needs a frame with z along the pole). 1 sigma per axis: 100 m and
    // 10 cm/s for the first object, 1 km and 1 m/s for the debris.
    Mat6 P1{}, P2{};
    for (std::size_t i = 0; i < 3; ++i) {
        P1[i][i]         = 100.0 * 100.0;
        P1[i + 3][i + 3] = 0.1 * 0.1;
        P2[i][i]         = 1000.0 * 1000.0;
        P2[i + 3][i + 3] = 1.0 * 1.0;
    }
    const UncertainTrajectory first{epoch, s0, P1}, second{epoch, elementsToState(debris), P2};
    for (const ConjunctionAssessment& c : assessConjunctions(first, second, epoch, end, 6.0, 10e3)) {
        std::cout << "  " << c.approach.time.iso8601(0) << "  " << c.approach.distance / 1e3 << " km, Pc "
                  << std::scientific << std::setprecision(2) << c.collision.probability << std::fixed
                  << std::setprecision(1) << " (" << c.collision.mahalanobis << " sigma)\n";
    }
    //! [conjunction_assessment]
    return 0;
}
