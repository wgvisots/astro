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

// Propagation: two-body, analytic J2/J4 (secular), numerical J2/J4 on an output grid, a
// 12-component (position, velocity, acceleration, jerk) state with an unmodeled force, and the
// closest approach of two objects.

#include <iomanip>
#include <iostream>

#include "astro/core/angle.h"
#include "astro/orbit/close_approach.h"
#include "astro/orbit/elements.h"
#include "astro/orbit/numerical.h"
#include "astro/orbit/secular.h"
#include "astro/orbit/two_body.h"

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
    return 0;
}
