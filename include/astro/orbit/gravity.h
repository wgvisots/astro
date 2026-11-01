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

#include <vector>

#include "astro/core/constants.h"
#include "astro/core/vec.h"

/// Gravitational accelerations. Positions are in a frame whose z axis is the
/// body's rotation axis (for the Earth, e.g. CIRS/TOD/ITRS, or GCRS as an approximation).
namespace astro {

/// Point-mass acceleration -mu r / |r|^3.
[[nodiscard]] Vec3 pointMassAcceleration(const Vec3& r, double mu = constants::gmEarth) noexcept;

/// Acceleration from the zonal harmonics J2..Jn (point mass excluded).
/// `J` holds unnormalized zonal coefficients indexed by degree: J[2] = J2, J[3] = J3, ...
/// (entries 0 and 1 are ignored).
[[nodiscard]] Vec3 zonalAcceleration(const Vec3& r, const std::vector<double>& J, double mu = constants::gmEarth,
                                     double re = constants::earthRadiusEgm2008) noexcept;

/// Closed-form J2 acceleration (for reference/testing).
[[nodiscard]] Vec3 j2Acceleration(const Vec3& r, double J2 = constants::J2, double mu = constants::gmEarth,
                                  double re = constants::earthRadiusEgm2008) noexcept;

/// Point-mass third-body perturbation on a satellite at r (geocentric) from a body at s (geocentric).
[[nodiscard]] Vec3 thirdBodyAcceleration(const Vec3& r, const Vec3& s, double muBody) noexcept;

} // namespace astro
