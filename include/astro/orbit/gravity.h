/**********************************************************************
 * gravity.h
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

#include "astro/core/constants.h"
#include "astro/core/vec.h"
#include "astro/orbit/central_body.h"

/// \file
/// Gravitational accelerations. Positions are in a frame whose z axis is the
/// body's rotation axis (for the Earth: CIRS, TOD, TEME or ITRS; rotate GCRS or J2000 vectors with
/// poleOfDateMatrix first: GCRS z is ~0.15 deg from the true pole).
namespace astro {

/// Point-mass acceleration -mu r / |r|^3.
[[nodiscard]] Vec3 pointMassAcceleration(const Vec3& r, double mu = constants::gmEarth) noexcept;

/// Point-mass jerk (time derivative of pointMassAcceleration along the motion):
/// -mu / |r|^3 (v - 3 (r.v) r / |r|^2).
[[nodiscard]] Vec3 pointMassJerk(const Vec3& r, const Vec3& v, double mu = constants::gmEarth) noexcept;

/// Acceleration from the zonal harmonics J2..Jn of `body` (point mass excluded): every entry of
/// body.J from degree 2 up, with body.mu and body.re. Does not allocate.
[[nodiscard]] Vec3 zonalAcceleration(const Vec3& r, const CentralBody& body = {}) noexcept;

/// Gradient d(pointMassAcceleration)/dr = mu / |r|^3 (3 r_hat r_hat^T - I).
[[nodiscard]] Mat3 pointMassGradient(const Vec3& r, double mu = constants::gmEarth) noexcept;

/// Gradient d(zonalAcceleration)/dr (symmetric, exact).
[[nodiscard]] Mat3 zonalGradient(const Vec3& r, const CentralBody& body = {}) noexcept;

/// Time derivative of zonalAcceleration along the motion (r moving with velocity v). Exact.
[[nodiscard]] Vec3 zonalJerk(const Vec3& r, const Vec3& v, const CentralBody& body = {}) noexcept;

/// zonalAcceleration and zonalJerk together, sharing the Legendre recursion.
void zonalAccelerationAndJerk(const Vec3&        r,
                              const Vec3&        v,
                              Vec3&              acceleration,
                              Vec3&              jerk,
                              const CentralBody& body = {}) noexcept;

/// Closed-form J2 acceleration (for reference/testing): body.zonal(2) with body.mu and body.re;
/// the other harmonics of `body` are ignored.
[[nodiscard]] Vec3 j2Acceleration(const Vec3& r, const CentralBody& body = {}) noexcept;

/// Point-mass third-body perturbation on a satellite at r (geocentric) from a body at s (geocentric).
[[nodiscard]] Vec3 thirdBodyAcceleration(const Vec3& r, const Vec3& s, double muBody) noexcept;

} // namespace astro
