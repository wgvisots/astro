/**********************************************************************
 * two_body.h
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

#include <vector>

#include "astro/core/constants.h"
#include "astro/core/vec.h"
#include "astro/orbit/state12.h"

/// Two-body (Keplerian) propagation and Lambert's problem.
namespace astro {

/// Stumpff functions c2(psi), c3(psi).
void stumpff(double psi, double& c2, double& c3) noexcept;

/// Propagate a two-body state by dt seconds using universal variables
/// (Vallado's KEPLER algorithm). Valid for elliptic, parabolic and hyperbolic orbits.
[[nodiscard]] StateVector propagateKepler(const StateVector& s0, double dt, double mu = constants::gmEarth);

/// Two-body propagation of a 12-component state (position, velocity, acceleration, jerk) by dt
/// seconds. The output acceleration and jerk are those of the point mass (plus the unmodeled
/// force, if any) at the new time.
///  - ModelOnly: universal variables from position and velocity.
///  - UnmodeledForce: the residuals against the point mass are applied as a force linear in time,
///    integrated numerically with RKF7(8). If both residuals are exactly zero it uses universal
///    variables.
[[nodiscard]] Propagation12 propagateKepler(const StateVector12& s0,
                                            double               dt,
                                            double               mu   = constants::gmEarth,
                                            HigherOrderMode      mode = HigherOrderMode::ModelOnly);

struct LambertSolution {
    Vec3 v1;                 ///< departure velocity
    Vec3 v2;                 ///< arrival velocity
    int  revolutions{0};     ///< complete revolutions
    bool rightBranch{false}; ///< for multi-revolution solutions: right (long-period) branch
};

/// Lambert's problem by Izzo's (2015) method: all solutions for 0..maxRevolutions
/// complete revolutions transferring from r1 to r2 in time `tof`.
/// The transfer direction is prograde (counter-clockwise about +z) unless
/// `prograde` is false. Multi-revolution solutions come in pairs.
[[nodiscard]] std::vector<LambertSolution> lambertIzzo(const Vec3& r1,
                                                       const Vec3& r2,
                                                       double      tof,
                                                       double      mu             = constants::gmEarth,
                                                       bool        prograde       = true,
                                                       int         maxRevolutions = 0);

/// The zero-revolution Lambert solution.
[[nodiscard]] LambertSolution
lambert(const Vec3& r1, const Vec3& r2, double tof, double mu = constants::gmEarth, bool prograde = true);

} // namespace astro
