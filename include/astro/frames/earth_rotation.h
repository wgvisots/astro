/**********************************************************************
 * earth_rotation.h
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

#include "astro/core/vec.h"
#include "astro/frames/precession_nutation.h"
#include "astro/time/julian_date.h"

/// Earth rotation angle, sidereal time and polar motion. Angles in radians [0, 2pi).
namespace astro {

/// Nominal Earth rotation rate implied by the ERA definition, rad/s (SI, UT1 seconds).
inline constexpr double kEarthRotationRateEra = 7.292115146706979e-5;

/// Earth Rotation Angle (IAU 2000), from UT1 as a two-part Julian Date.
[[nodiscard]] double earthRotationAngle(const JulianDate& ut1) noexcept;

/// Greenwich Mean Sidereal Time, IAU 2006 (consistent with IAU 2006 precession).
[[nodiscard]] double gmst06(const JulianDate& ut1, const JulianDate& tt) noexcept;

/// Equation of the origins, from a bias-precession-nutation matrix and the CIO locator s.
[[nodiscard]] double equationOfOrigins(const Mat3& npb, double s) noexcept;

/// Greenwich Apparent Sidereal Time, IAU 2006/2000A (or 2000B), CIO based: GAST = ERA - EO.
[[nodiscard]] double gast06(const JulianDate&       ut1,
                            const JulianDate&       tt,
                            PrecessionNutationModel model = PrecessionNutationModel::IAU2006_2000A) noexcept;

/// Greenwich Mean Sidereal Time, IAU 1982 model (used to define TEME for SGP4).
[[nodiscard]] double gmst82(const JulianDate& ut1) noexcept;

/// Polar motion matrix W^T: TIRS -> ITRS, R1(-yp) R2(-xp) R3(s').
[[nodiscard]] Mat3 polarMotionMatrix(double xp, double yp, double sp) noexcept;

} // namespace astro
