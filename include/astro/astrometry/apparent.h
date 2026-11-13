/**********************************************************************
 * apparent.h
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

#include "astro/astrometry/refraction.h"
#include "astro/astrometry/spherical.h"
#include "astro/ephemeris/ephemeris.h"
#include "astro/frames/eop.h"
#include "astro/frames/transform.h"
#include "astro/geodesy/local.h"

/// \file
/// Apparent and observed places of stars and solar-system bodies.
///
/// Pipeline (as in ERFA, the open-source derivative of SOFA): catalog ICRS place -> space motion -> parallax
/// -> gravitational light deflection by the Sun -> aberration (annual and
/// diurnal, from the observer's barycentric velocity) -> GCRS apparent
/// direction -> ITRS -> local horizon -> optional refraction.
namespace astro {

/// A catalog star (ICRS).
struct CatalogStar {
    double ra{0.0};                        ///< right ascension at the catalog epoch, rad
    double dec{0.0};                       ///< declination at the catalog epoch, rad
    double pmRaCosDec{0.0};                ///< proper motion in RA times cos(dec), rad / Julian year
    double pmDec{0.0};                     ///< rad / Julian year
    double parallax{0.0};                  ///< rad (0 = unknown / very distant)
    double radialVelocity{0.0};            ///< m/s, positive receding
    double epochJdTdb{constants::jdJ2000}; ///< catalog epoch (e.g. 2457389.0 for Gaia DR3's J2016.0)
};

/// Barycentric state and heliocentric position of an observer (m, m/s, ICRF axes).
struct ObserverEphemeris {
    StateVector barycentric;  ///< position and velocity relative to the solar-system barycenter
    Vec3        heliocentric; ///< position relative to the Sun
};

/// Ephemeris of the geocenter, plus an optional geocentric observer offset
/// (GCRS position/velocity). If the ephemeris cannot supply the solar-system
/// barycenter (e.g. LowPrecisionEphemeris), heliocentric is used in its place
/// (~13 m/s error in velocity, ~0.01 arcsec in aberration). Throws
/// astro::Error if `eph` cannot give the Earth relative to the Sun at `t`.
[[nodiscard]] ObserverEphemeris
observerEphemeris(const Time& t, const Ephemeris& eph, const StateVector& geocentricOffset = {});

/// Gravitational deflection of direction p by a body of `massSolar` solar masses.
/// q: unit vector from the body to the source; e: unit vector from the body to the
/// observer; emAu: body-observer distance (au). (ERFA eraLd.)
[[nodiscard]] Vec3
lightDeflection(const Vec3& p, const Vec3& q, const Vec3& e, double emAu, double massSolar = 1.0) noexcept;

/// Stellar aberration: natural direction p -> proper direction, given observer
/// barycentric velocity in units of c and Sun-observer distance (au). (ERFA eraAb.)
[[nodiscard]] Vec3 aberration(const Vec3& p, const Vec3& vOverC, double sunDistanceAu) noexcept;

/// Star position (ICRS, at time t) after space motion: barycentric unit vector or, if
/// the parallax is known, barycentric position in meters. Throws astro::Error if `t`
/// cannot be converted to TDB: a UTC time before 1960, or a UT1 time.
[[nodiscard]] Vec3 starBarycentricPosition(const CatalogStar& star, const Time& t);

/// Apparent (proper) direction of a star for an observer, GCRS axes (unit vector).
/// Throws astro::Error as starBarycentricPosition does.
[[nodiscard]] Vec3 starApparentDirection(const CatalogStar& star, const Time& t, const ObserverEphemeris& obs);
/// Geocentric apparent RA/Dec in the GCRS. Throws astro::Error as
/// starBarycentricPosition and observerEphemeris do.
[[nodiscard]] RaDec starApparentGcrs(const CatalogStar& star, const Time& t, const Ephemeris& eph);

/// Observed azimuth / elevation of a star from a site. With `atm` null, the
/// elevation is unrefracted; range is set to 0 (unknown) for stars. Throws
/// astro::Error if `eop` has no data at `t` or `t` is before 1960 (see
/// FrameTransformer), and as starApparentGcrs does.
[[nodiscard]] Aer starObserved(const CatalogStar& star,
                               const Time&        t,
                               const Geodetic&    site,
                               const EopProvider& eop,
                               const Ephemeris&   eph,
                               const Atmosphere*  atm = nullptr);

/// Observed azimuth / elevation / range of a solar-system body from a site,
/// including light time and aberration (light deflection is neglected, < 2 arcsec).
/// Throws astro::Error if `eop` has no data at `t` or `t` is before 1960 (see
/// FrameTransformer), as observerEphemeris does, and if `eph` has no data for
/// `body` at `t` less the light time.
[[nodiscard]] Aer bodyObserved(Body               body,
                               const Time&        t,
                               const Geodetic&    site,
                               const EopProvider& eop,
                               const Ephemeris&   eph,
                               const Atmosphere*  atm = nullptr);

/// Geometric topocentric az/el/range of a GCRS position (no light time or aberration).
[[nodiscard]] Aer topocentricFromGcrs(const Vec3&             rGcrs,
                                      const FrameTransformer& ft,
                                      const Geodetic&         site,
                                      const Ellipsoid&        ell = Ellipsoid::wgs84());

/// Local mean / apparent sidereal time (rad, [0, 2pi)) at east longitude `lon`
/// (rad): GMST (IAU 2006) / GAST (IAU 2006/2000A) plus `lon`. `t` may be in any
/// time scale; throws astro::Error if `eop` has no data at `t`, or if `t` is
/// before 1960 (see FrameTransformer).
/// @{
[[nodiscard]] double localMeanSiderealTime(const Time& t, const EopProvider& eop, double lon);
[[nodiscard]] double localApparentSiderealTime(const Time& t, const EopProvider& eop, double lon);
/// @}
/// Hour angle (rad, [-pi, pi)) from local sidereal time and right ascension.
[[nodiscard]] double hourAngle(double lst, double ra) noexcept;

} // namespace astro
