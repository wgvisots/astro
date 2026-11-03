/**********************************************************************
 * geodetic.h
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
#include "astro/geodesy/ellipsoid.h"

/// \file
/// Geodetic latitude/longitude/height <-> Earth-centred Earth-fixed coordinates, and auxiliary latitudes.
namespace astro {

/// Geodetic coordinates: latitude and longitude (rad), height above the ellipsoid (m).
struct Geodetic {
    double lat{0.0};
    double lon{0.0};
    double h{0.0};
};

/// Geodetic (lat, lon, h) to Earth-centered Earth-fixed Cartesian (m).
[[nodiscard]] Vec3 geodeticToEcef(const Geodetic& g, const Ellipsoid& ell = Ellipsoid::wgs84()) noexcept;

/// ECEF (m) to geodetic. Closed-form (Vermeille 2004) for all points outside the
/// ellipsoid's evolute (i.e. everywhere farther than ~43 km from the Earth's
/// center for WGS84), with a Newton fallback inside it. Longitude in (-pi, pi];
/// at the poles longitude is 0.
[[nodiscard]] Geodetic ecefToGeodetic(const Vec3& r, const Ellipsoid& ell = Ellipsoid::wgs84()) noexcept;

/// Geocentric latitude of a point on the ellipsoid surface with the given geodetic latitude.
[[nodiscard]] double geodeticToGeocentricLatitude(double           geodeticLat,
                                                  const Ellipsoid& ell = Ellipsoid::wgs84()) noexcept;
/// Geodetic latitude of a surface point with the given geocentric latitude.
[[nodiscard]] double geocentricToGeodeticLatitude(double           geocentricLat,
                                                  const Ellipsoid& ell = Ellipsoid::wgs84()) noexcept;
/// Reduced (parametric) latitude beta: tan(beta) = (1 - f) tan(phi).
[[nodiscard]] double geodeticToReducedLatitude(double geodeticLat, const Ellipsoid& ell = Ellipsoid::wgs84()) noexcept;
[[nodiscard]] double reducedToGeodeticLatitude(double reducedLat, const Ellipsoid& ell = Ellipsoid::wgs84()) noexcept;

/// Spherical (geocentric) coordinates of a vector: latitude, longitude, radius.
struct Spherical {
    double lat{0.0};
    double lon{0.0};
    double r{0.0};
};
[[nodiscard]] Spherical cartesianToSpherical(const Vec3& v) noexcept;
[[nodiscard]] Vec3      sphericalToCartesian(const Spherical& s) noexcept;

} // namespace astro
