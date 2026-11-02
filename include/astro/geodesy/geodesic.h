/**********************************************************************
 * geodesic.h
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

#include "astro/geodesy/ellipsoid.h"

/// Geodesics on the ellipsoid (Vincenty 1975). Sub-millimeter accurate; the
/// inverse problem may fail to converge for nearly antipodal points, in which
/// case astro::Error is thrown.
namespace astro {

struct GeodesicInverseResult {
    double distance{0.0}; ///< m
    double azimuth1{0.0}; ///< forward azimuth at point 1, rad [0, 2pi)
    double azimuth2{0.0}; ///< forward azimuth at point 2 (direction of travel), rad [0, 2pi)
};

struct GeodesicDirectResult {
    double lat{0.0};      ///< rad
    double lon{0.0};      ///< rad, (-pi, pi]
    double azimuth2{0.0}; ///< forward azimuth at the end point, rad [0, 2pi)
};

/// Distance and azimuths between two points (lat/lon in radians).
[[nodiscard]] GeodesicInverseResult
geodesicInverse(double lat1, double lon1, double lat2, double lon2, const Ellipsoid& ell = Ellipsoid::wgs84());

/// Destination after travelling `distance` meters from (lat1, lon1) along initial azimuth `azimuth1`.
[[nodiscard]] GeodesicDirectResult
geodesicDirect(double lat1, double lon1, double azimuth1, double distance, const Ellipsoid& ell = Ellipsoid::wgs84());

/// Great-circle (haversine) distance on a sphere of radius r.
[[nodiscard]] double haversineDistance(double lat1, double lon1, double lat2, double lon2, double r) noexcept;

} // namespace astro
