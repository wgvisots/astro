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

#include <limits>

#include "astro/core/vec.h"
#include "astro/geodesy/geodetic.h"

/// \file
/// Local topocentric frames: East-North-Up (ENU), North-East-Down (NED), and look angles
/// (azimuth, elevation, range and range rate). Azimuth is measured clockwise from north.
namespace astro {

/// Rotation matrix taking ECEF vector components to ENU at (lat, lon).
[[nodiscard]] Mat3 ecefToEnuMatrix(double lat, double lon) noexcept;
/// Rotation matrix taking ECEF vector components to NED at (lat, lon).
[[nodiscard]] Mat3 ecefToNedMatrix(double lat, double lon) noexcept;

/// ECEF position to ENU relative to a geodetic origin.
[[nodiscard]] Vec3
ecefToEnu(const Vec3& ecef, const Geodetic& origin, const Ellipsoid& ell = Ellipsoid::wgs84()) noexcept;
/// ENU position relative to a geodetic origin to ECEF (m).
[[nodiscard]] Vec3
enuToEcef(const Vec3& enu, const Geodetic& origin, const Ellipsoid& ell = Ellipsoid::wgs84()) noexcept;
/// ECEF position to NED relative to a geodetic origin (m).
[[nodiscard]] Vec3
ecefToNed(const Vec3& ecef, const Geodetic& origin, const Ellipsoid& ell = Ellipsoid::wgs84()) noexcept;
/// NED position relative to a geodetic origin to ECEF (m).
[[nodiscard]] Vec3
nedToEcef(const Vec3& ned, const Geodetic& origin, const Ellipsoid& ell = Ellipsoid::wgs84()) noexcept;

/// Topocentric look angles from an observer to a target, relative to the observer's geodetic
/// horizon (geometric: no refraction or light time, unless the function that returns them says so).
/// The range rate is NaN unless the function was given a velocity: a silent 0 would read as "not
/// moving".
struct LookAngles {
    double az{0.0};    ///< azimuth, rad, [0, 2pi), clockwise from north
    double el{0.0};    ///< elevation above the local horizontal, rad, [-pi/2, pi/2]
    double range{0.0}; ///< m (0 when unknown, as for a star)
    /// range rate, m/s, positive receding; NaN if not computed (no velocity given)
    double rangeRate{std::numeric_limits<double>::quiet_NaN()};
};

/// Azimuth, elevation and range of an ENU vector `enu` (m); the range rate is NaN.
[[nodiscard]] LookAngles enuToLookAngles(const Vec3& enu) noexcept;
/// ENU vector (m) from the azimuth, elevation and range of `look`; its range rate is not used.
[[nodiscard]] Vec3 lookAnglesToEnu(const LookAngles& look) noexcept;
/// Look angles from an observer to an ECEF target position `target` (m); the range rate is NaN.
[[nodiscard]] LookAngles
ecefToLookAngles(const Vec3& target, const Geodetic& observer, const Ellipsoid& ell = Ellipsoid::wgs84()) noexcept;
/// Look angles, with the range rate, from an observer to a target whose state `target` is given in
/// the ITRS (ECEF; m, m/s). For a target exactly at the observer (zero range) the direction and the
/// range rate do not exist: the range rate is NaN.
[[nodiscard]] LookAngles ecefToLookAngles(const StateVector& target,
                                          const Geodetic&    observer,
                                          const Ellipsoid&   ell = Ellipsoid::wgs84()) noexcept;
/// ECEF position (m) of a target seen at the azimuth, elevation and range of `look` from an
/// observer; its range rate is not used.
[[nodiscard]] Vec3
lookAnglesToEcef(const LookAngles& look, const Geodetic& observer, const Ellipsoid& ell = Ellipsoid::wgs84()) noexcept;

} // namespace astro
