/**********************************************************************
 * local.h
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
#include "astro/geodesy/geodetic.h"

/// \file
/// Local topocentric frames: East-North-Up (ENU), North-East-Down (NED), and
/// azimuth / elevation / range (AER). Azimuth is measured clockwise from north.
namespace astro {

/// Rotation matrix taking ECEF vector components to ENU at (lat, lon).
[[nodiscard]] Mat3 ecefToEnuMatrix(double lat, double lon) noexcept;
/// Rotation matrix taking ECEF vector components to NED at (lat, lon).
[[nodiscard]] Mat3 ecefToNedMatrix(double lat, double lon) noexcept;

/// ECEF position to ENU relative to a geodetic origin.
[[nodiscard]] Vec3
ecefToEnu(const Vec3& ecef, const Geodetic& origin, const Ellipsoid& ell = Ellipsoid::wgs84()) noexcept;
[[nodiscard]] Vec3
enuToEcef(const Vec3& enu, const Geodetic& origin, const Ellipsoid& ell = Ellipsoid::wgs84()) noexcept;
[[nodiscard]] Vec3
ecefToNed(const Vec3& ecef, const Geodetic& origin, const Ellipsoid& ell = Ellipsoid::wgs84()) noexcept;
[[nodiscard]] Vec3
nedToEcef(const Vec3& ned, const Geodetic& origin, const Ellipsoid& ell = Ellipsoid::wgs84()) noexcept;

/// Azimuth (rad, [0, 2pi), clockwise from north), elevation (rad), range (m).
struct Aer {
    double az{0.0};
    double el{0.0};
    double range{0.0};
};

[[nodiscard]] Aer  enuToAer(const Vec3& enu) noexcept;
[[nodiscard]] Vec3 aerToEnu(const Aer& aer) noexcept;
/// Look angles from an observer to an ECEF target.
[[nodiscard]] Aer
ecefToAer(const Vec3& target, const Geodetic& observer, const Ellipsoid& ell = Ellipsoid::wgs84()) noexcept;
[[nodiscard]] Vec3
aerToEcef(const Aer& aer, const Geodetic& observer, const Ellipsoid& ell = Ellipsoid::wgs84()) noexcept;

} // namespace astro
