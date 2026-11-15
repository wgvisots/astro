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
#include "astro/export.h"
#include "astro/orbit/central_body.h"

/// \file
/// Gravitational accelerations. Positions are in a frame whose z axis is the
/// body's rotation axis (for the Earth: CIRS, TOD, TEME or ITRS; rotate GCRS or J2000 vectors with
/// poleOfDateMatrix first: GCRS z is ~0.15 deg from the true pole).
namespace astro {

/// Point-mass acceleration -mu r / |r|^3.
[[nodiscard]] ASTRO_EXPORT Vec3 pointMassAcceleration(const Vec3& r, double mu = constants::gmEarth) noexcept;

/// Acceleration from the zonal harmonics J2..Jn of `body` (point mass excluded): every entry of
/// body.J from degree 2 up, with body.mu and body.re. Does not allocate.
[[nodiscard]] ASTRO_EXPORT Vec3 zonalAcceleration(const Vec3&        r,
                                                  const CentralBody& body = CentralBody::earth()) noexcept;

/// Point-mass third-body perturbation on a satellite at r (geocentric) from a body at s (geocentric).
[[nodiscard]] ASTRO_EXPORT Vec3 thirdBodyAcceleration(const Vec3& r, const Vec3& s, double muBody) noexcept;

} // namespace astro
