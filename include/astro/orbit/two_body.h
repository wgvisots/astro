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
#include "astro/orbit/state_covariance.h"

/// \file
/// Two-body (Keplerian) propagation and Lambert's problem.
namespace astro {

/// The two Stumpff functions of one argument, as stumpff() returns them.
struct Stumpff {
    double c2{0.0}; ///< c2(psi), dimensionless
    double c3{0.0}; ///< c3(psi), dimensionless
};

/// Stumpff functions c2(psi) = (1 - cos sqrt(psi)) / psi and c3(psi) = (sqrt(psi) - sin sqrt(psi)) / psi^(3/2),
/// continued to psi <= 0 (hyperbolic functions). `psi` is dimensionless (chi^2 alpha in the universal
/// variables). Accurate to a few rounding errors for every psi: by their power series for |psi| < 1, where the
/// closed forms cancel.
[[nodiscard]] Stumpff stumpff(double psi) noexcept;

/// Propagate a two-body state by dt seconds using universal variables
/// (Vallado's KEPLER algorithm). Valid for elliptic, parabolic and hyperbolic orbits, any dt of
/// either sign, and any step size down to rounding.
///
/// The universal Kepler equation is solved by the Laguerre-Conway iteration inside a bracket, so it
/// always ends; the result is within a few rounding errors of the equation's own conditioning
/// (1e-13 of the state for steps up to a few revolutions; see docs/accuracy.md).
///
/// Throws astro::Error if the state or dt is not finite, the position is zero, or mu is not positive.
/// With dt = 0 the state is returned as given, unchecked.
[[nodiscard]] StateVector propagateKepler(const StateVector& s0, double dt, double mu = constants::gmEarth);

/// propagateKepler with the state transition matrix Phi(dt) = d(state at dt) / d(s0), in closed
/// form (the quantity of Shepperd 1985 and Battin 1999): the universal-variable Lagrange
/// coefficients differentiated with respect to r0, v0 and, implicitly through Kepler's equation,
/// the universal anomaly over the full dt (whole revolutions included, which carry the
/// along-track drift). Exact to rounding and valid for any conic: equal to the integrated
/// variational equations to ~1e-11, and symplectic to 6e-9 after a day in LEO and 1e-7 after a week
/// (astro_study_propagation). The state is
/// propagateKepler's, bit for bit. Phi(0) is the identity. Throws astro::Error as propagateKepler does.
[[nodiscard]] StateTransition propagateKeplerStm(const StateVector& s0, double dt, double mu = constants::gmEarth);

/// Two-body state and covariance dt seconds after s0: P = Phi P0 Phi^T, Phi from propagateKeplerStm.
/// P0 is the covariance of s0 (m^2, m^2/s, m^2/s^2) in the frame of s0, and mu the gravitational
/// parameter of the central body (m^3/s^2). No process noise; for that,
/// or for any force beyond the point mass, use NumericalPropagator::propagateWithCovariance.
/// Throws astro::Error as propagateKepler does.
[[nodiscard]] StateCovariance
propagateKeplerWithCovariance(const StateVector& s0, double dt, const Mat6& P0, double mu = constants::gmEarth);

/// Two-body propagation of a 12-component state (position, velocity, acceleration, jerk) by dt
/// seconds. The output acceleration and jerk are those of the point mass (plus the unmodeled
/// force, if any) at the new time.
///  - ModelOnly: universal variables from position and velocity.
///  - UnmodeledForce: the residuals against the point mass are applied as a force linear in time,
///    integrated numerically with RKF7(8). If both residuals are exactly zero it uses universal
///    variables.
///
/// Throws astro::Error as the 6-component propagateKepler does or, where the unmodeled force is
/// integrated, as NumericalPropagator::propagate does.
[[nodiscard]] Propagation12 propagateKepler(const StateVector12& s0,
                                            double               dt,
                                            double               mu   = constants::gmEarth,
                                            HigherOrderMode      mode = HigherOrderMode::ModelOnly);

/// One solution of Lambert's problem; velocities (m/s) in the frame of r1 and r2.
struct LambertSolution {
    Vec3 v1;                 ///< departure velocity
    Vec3 v2;                 ///< arrival velocity
    int  revolutions{0};     ///< complete revolutions
    bool rightBranch{false}; ///< for multi-revolution solutions: the one with the larger semi-major axis
};

/// Lambert's problem: the two-body orbits from r1 to r2 (m) in the time `tof` (s), by the method of
/// Izzo (2015). Returns the solution without a complete revolution first, then, for each number of
/// revolutions M = 1 .. maxRevolutions that the time allows, the two solutions of M revolutions
/// (rightBranch false, then true).
///
/// The method iterates on the Lancaster-Blanchard variable x with Householder's third-order scheme
/// from Izzo's starting values; the non-dimensional time of flight comes from his eq. 18, and from
/// Battin's hypergeometric series near the parabola. Each root is checked against the time of
/// flight (to 1e-10) and a bracketed search stands behind the iteration. Typical accuracy: the
/// arrival position is reproduced to ~1e-12 relative, less where the problem itself is sensitive
/// (see below and docs/accuracy.md).
///
/// `prograde` selects the direction of motion: counter-clockwise seen from +z (the orbit's angular
/// momentum has a positive z component) or, if false, clockwise. If the plane of r1 and r2 contains
/// the z axis (the z component of r1 x r2 is zero to within its rounding error) the two directions
/// cannot be told apart that way, and the shorter way round is taken.
///
/// Sensitive cases, by the nature of the problem: r1 and r2 within ~1e-6 rad of opposite (the orbit
/// plane is fixed only to rounding / sin(angle)); flights so long that the orbit's period is not
/// resolved by a double-precision velocity (non-dimensional times above ~1e5); and the long way
/// round in a very short time, which passes the centre of attraction closely.
///
/// Throws astro::Error if tof or mu is not positive and finite, if r1 or r2 is zero or not finite,
/// or if r1 and r2 are collinear (sine of the angle between them below 1e-12: the plane is
/// undefined). Also throws if no root is found for a number of revolutions that should have one: a
/// guard against a failed iteration.
[[nodiscard]] std::vector<LambertSolution> lambertIzzo(const Vec3& r1,
                                                       const Vec3& r2,
                                                       double      tof,
                                                       double      mu             = constants::gmEarth,
                                                       bool        prograde       = true,
                                                       int         maxRevolutions = 0);

/// The Lambert solution without a complete revolution (the first of lambertIzzo's). Throws
/// astro::Error as lambertIzzo does.
[[nodiscard]] LambertSolution
lambert(const Vec3& r1, const Vec3& r2, double tof, double mu = constants::gmEarth, bool prograde = true);

} // namespace astro
