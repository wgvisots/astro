/**********************************************************************
 * spherical.h
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

/// \file
/// Celestial spherical coordinates and the classical coordinate systems.
namespace astro {

/// Right ascension [0, 2pi), declination [-pi/2, pi/2] (rad), and distance (any unit).
struct RaDec {
    double ra{0.0};       ///< right ascension, rad, [0, 2pi)
    double dec{0.0};      ///< declination, rad, [-pi/2, pi/2]
    double distance{1.0}; ///< distance (any unit; 1 for a direction)
};

/// RA/Dec and length of a vector; RA is 0 on the z axis.
[[nodiscard]] RaDec vectorToRaDec(const Vec3& v) noexcept;
/// Vector of length `rd.distance` in the direction (ra, dec).
[[nodiscard]] Vec3 raDecToVector(const RaDec& rd) noexcept;
/// Unit vector in the direction (ra, dec), rad.
[[nodiscard]] inline Vec3 raDecToUnitVector(double ra, double dec) noexcept {
    return raDecToVector({ra, dec, 1.0});
}

/// Angular separation between two directions (rad), robust for all angles.
[[nodiscard]] double angularSeparation(double ra1, double dec1, double ra2, double dec2) noexcept;
/// Position angle of point 2 from point 1, east of north (rad, [0, 2pi)).
[[nodiscard]] double positionAngle(double ra1, double dec1, double ra2, double dec2) noexcept;

/// Equatorial -> ecliptic rotation for obliquity eps: v_ecl = R1(eps) v_eq.
[[nodiscard]] Mat3 equatorialToEclipticMatrix(double eps) noexcept;
/// Ecliptic coordinates of an equatorial vector (lon, lat as RaDec fields ra/dec).
[[nodiscard]] RaDec equatorialToEcliptic(const Vec3& eq, double eps) noexcept;
/// Equatorial vector from ecliptic coordinates (lon, lat, distance as RaDec
/// fields ra/dec/distance): v_eq = R1(-eps) v_ecl.
[[nodiscard]] Vec3 eclipticToEquatorial(const RaDec& ecl, double eps) noexcept;

/// ICRS -> Galactic rotation matrix (Hipparcos definition, ESA 1997).
[[nodiscard]] const Mat3& icrsToGalacticMatrix() noexcept;
/// Galactic longitude/latitude (as RaDec ra/dec fields) from ICRS RA/Dec.
[[nodiscard]] RaDec icrsToGalactic(double ra, double dec) noexcept;
/// ICRS RA/Dec (as RaDec ra/dec fields, unit distance) from Galactic longitude
/// `l` and latitude `b` (rad).
[[nodiscard]] RaDec galacticToIcrs(double l, double b) noexcept;

} // namespace astro
