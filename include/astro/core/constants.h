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
namespace astro::constants {

// ------------------------------------------------------------ math
inline constexpr double pi = std::numbers::pi;
inline constexpr double twoPi = 2.0 * std::numbers::pi;
inline constexpr double halfPi = 0.5 * std::numbers::pi;
inline constexpr double degToRad = pi / 180.0;
inline constexpr double radToDeg = 180.0 / pi;
inline constexpr double arcsecToRad = pi / 648000.0;
inline constexpr double masToRad = arcsecToRad / 1000.0;
inline constexpr double uasToRad = arcsecToRad / 1.0e6;
/// Turns (revolutions) to radians.
inline constexpr double turnToRad = twoPi;

// ------------------------------------------------------------ time
inline constexpr double secondsPerDay = 86400.0;
inline constexpr double daysPerJulianYear = 365.25;
inline constexpr double daysPerJulianCentury = 36525.0;
inline constexpr double daysPerJulianMillennium = 365250.0;
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
/// JD (TAI) of 1977-01-01T00:00:00 TAI, i.e. 1977-01-01T00:00:32.184 TT.
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
/// Geocentric gravitational constant, m^3/s^2 (IERS 2010, TT-compatible; same as WGS84 EGM).
inline constexpr double gmEarth = 3.986004418e14;
/// Nominal mean angular velocity of the Earth, rad/s (IERS 2010 / WGS84).
inline constexpr double earthRotationRate = 7.292115e-5;
/// Earth equatorial radius used by EGM2008, m.
inline constexpr double earthRadiusEgm2008 = 6378136.3;
/// Earth equatorial radius, IERS 2010 conventional, m.
inline constexpr double earthRadiusIers2010 = 6378136.6;
/// Unnormalized zonal harmonics from EGM2008 (tide-free), Jn = -sqrt(2n+1) * Cbar_n0.
inline constexpr double J2 = 1.0826261738522227e-3;
inline constexpr double J3 = -2.5324105185677225e-6;
inline constexpr double J4 = -1.6198975999169731e-6;
inline constexpr double J5 = -2.2775359073083618e-7;
inline constexpr double J6 = 5.406665762838132e-7;

// ------------------------------------------------------------ ellipsoids
inline constexpr double wgs84A = 6378137.0;
inline constexpr double wgs84InvF = 298.257223563;
inline constexpr double wgs84F = 1.0 / wgs84InvF;
inline constexpr double grs80A = 6378137.0;
inline constexpr double grs80InvF = 298.257222101;
inline constexpr double grs80F = 1.0 / grs80InvF;
inline constexpr double iers2010A = 6378136.6;
inline constexpr double iers2010InvF = 298.25642;
inline constexpr double iers2010F = 1.0 / iers2010InvF;

// ------------------------------------------------------------ bodies
/// IAU 2015 nominal solar radius, m.
inline constexpr double sunRadius = 6.957e8;
/// IAU mean lunar radius, m.
inline constexpr double moonRadius = 1.7374e6;

} // namespace astro::constants
