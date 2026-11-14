/**********************************************************************
 * constants.h
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

#include <numbers>

/// Physical, astronomical and geodetic constants. SI units throughout.
/// \file
/// Physical, geodetic and astronomical constants (WGS84, GRS80, IERS 2010, EGM2008, DE440). SI units.
namespace astro::constants {

// ------------------------------------------------------------ math
inline constexpr double pi     = std::numbers::pi;       ///< pi
inline constexpr double twoPi  = 2.0 * std::numbers::pi; ///< 2 pi (one turn, rad)
inline constexpr double halfPi = 0.5 * std::numbers::pi; ///< pi / 2 (90 degrees, rad)
/// @name Angle factors
/// The `...ToRad` and `radTo...` constants are factors to multiply an angle by. The functions of the
/// same names in namespace astro (astro::degToRad(), astro::radToDeg(), ...) apply them to an angle.
/// @{
inline constexpr double degToRad    = pi / 180.0;           ///< Degrees to radians (a factor, not a function).
inline constexpr double radToDeg    = 180.0 / pi;           ///< Radians to degrees (a factor, not a function).
inline constexpr double arcsecToRad = pi / 648000.0;        ///< Arcseconds to radians (a factor, not a function).
inline constexpr double masToRad    = arcsecToRad / 1000.0; ///< Milliarcseconds to radians (a factor, not a function).
inline constexpr double uasToRad    = arcsecToRad / 1.0e6;  ///< Microarcseconds to radians (a factor, not a function).
/// @}
/// Turns (revolutions) to radians.
inline constexpr double turnToRad = twoPi;

// ------------------------------------------------------------ time
inline constexpr double secondsPerDay           = 86400.0;  ///< SI seconds in a day (86400 s).
inline constexpr double daysPerJulianYear       = 365.25;   ///< Days in a Julian year.
inline constexpr double daysPerJulianCentury    = 36525.0;  ///< Days in a Julian century.
inline constexpr double daysPerJulianMillennium = 365250.0; ///< Days in a Julian millennium.
/// Julian Date of the J2000.0 epoch (2000-01-01 12:00:00 TT).
inline constexpr double jdJ2000 = 2451545.0;
/// MJD = JD - mjdOffset.
inline constexpr double mjdOffset = 2400000.5;
/// MJD of the J2000.0 epoch.
inline constexpr double mjdJ2000 = 51544.5;
/// Julian Date of the Unix epoch, 1970-01-01T00:00:00 UTC.
inline constexpr double jdUnixEpoch = 2440587.5;
/// Julian Date of the GPS epoch, 1980-01-06T00:00:00 UTC.
inline constexpr double jdGpsEpoch = 2444244.5;
/// TT - TAI, seconds (exact).
inline constexpr double ttMinusTai = 32.184;
/// TAI - GPS, seconds (exact). Also applies to Galileo System Time.
inline constexpr double taiMinusGps = 19.0;
/// TAI - BDT (BeiDou Time), seconds (exact).
inline constexpr double taiMinusBdt = 33.0;
/// IAU 2000 Resolution B1.9: 1 - d(TT)/d(TCG).
inline constexpr double LG = 6.969290134e-10;
/// IAU 2006 Resolution B3: 1 - d(TDB)/d(TCB).
inline constexpr double LB = 1.550519768e-8;
/// IAU 2006 Resolution B3: TDB - TCB at T0, seconds.
inline constexpr double TDB0 = -6.55e-5;
/// Julian date in TT of the instant 1977-01-01T00:00:00 TAI = 1977-01-01T00:00:32.184 TT, at which
/// TT, TCG and TCB read the same. As one double it is good to 14 microseconds only.
inline constexpr double jdT0 = 2443144.5003725;

// ------------------------------------------------------------ physics
/// Speed of light in vacuum, m/s (exact).
inline constexpr double c = 299792458.0;
/// Astronomical unit, m (IAU 2012, exact).
inline constexpr double au = 149597870700.0;
/// Heliocentric gravitational constant, m^3/s^2 (TDB-compatible, DE440).
inline constexpr double gmSun = 1.32712440041279419e20;
/// Selenocentric gravitational constant, m^3/s^2 (DE440).
inline constexpr double gmMoon = 4.902800118457549e12;
/// Earth/Moon mass ratio (DE440).
inline constexpr double earthMoonMassRatio = 81.3005682214972154;

// ------------------------------------------------------------ Earth
/// Geocentric gravitational constant, m^3/s^2: the value of IERS Conventions (2010) Table 1.1,
/// which is TCG-compatible, and of WGS 84. The TT-compatible value, which is also the GM of
/// EGM2008, is 3.986004415e14: 7.5e-10 smaller, about 2 cm along the track per revolution in low
/// Earth orbit.
inline constexpr double gmEarth = 3.986004418e14;
/// Nominal mean angular velocity of the Earth, rad/s (IERS 2010 / WGS84).
inline constexpr double earthRotationRate = 7.292115e-5;
/// Earth equatorial radius used by EGM2008, m.
inline constexpr double earthRadiusEgm2008 = 6378136.3;
/// Earth equatorial radius, IERS 2010 conventional, m.
inline constexpr double earthRadiusIers2010 = 6378136.6;
/// Unnormalized zonal harmonics from EGM2008 (tide-free), Jn = -sqrt(2n+1) * Cbar_n0.
inline constexpr double J2 = 1.0826261738522227e-3;
inline constexpr double J3 = -2.5324105185677225e-6; ///< Unnormalized zonal harmonic J3 (EGM2008, tide-free).
inline constexpr double J4 = -1.6198975999169731e-6; ///< Unnormalized zonal harmonic J4 (EGM2008, tide-free).
inline constexpr double J5 = -2.2775359073083618e-7; ///< Unnormalized zonal harmonic J5 (EGM2008, tide-free).
inline constexpr double J6 = 5.406665762838132e-7;   ///< Unnormalized zonal harmonic J6 (EGM2008, tide-free).

// ------------------------------------------------------------ ellipsoids
inline constexpr double wgs84A       = 6378137.0;          ///< WGS84 semi-major axis, m.
inline constexpr double wgs84InvF    = 298.257223563;      ///< WGS84 inverse flattening 1/f.
inline constexpr double wgs84F       = 1.0 / wgs84InvF;    ///< WGS84 flattening f.
inline constexpr double grs80A       = 6378137.0;          ///< GRS80 semi-major axis, m.
inline constexpr double grs80InvF    = 298.257222101;      ///< GRS80 inverse flattening 1/f.
inline constexpr double grs80F       = 1.0 / grs80InvF;    ///< GRS80 flattening f.
inline constexpr double iers2010A    = 6378136.6;          ///< IERS 2010 conventional ellipsoid semi-major axis, m.
inline constexpr double iers2010InvF = 298.25642;          ///< IERS 2010 conventional ellipsoid inverse flattening 1/f.
inline constexpr double iers2010F    = 1.0 / iers2010InvF; ///< IERS 2010 conventional ellipsoid flattening f.

// ------------------------------------------------------------ bodies
/// IAU 2015 nominal solar radius, m.
inline constexpr double sunRadius = 6.957e8;
/// IAU mean lunar radius, m.
inline constexpr double moonRadius = 1.7374e6;

} // namespace astro::constants
