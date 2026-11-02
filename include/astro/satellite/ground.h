/**********************************************************************
 * ground.h
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

#include <functional>
#include <vector>

#include "astro/core/vec.h"
#include "astro/geodesy/geodetic.h"
#include "astro/time/time.h"

/// Ground track, sub-satellite point, coverage footprint and line of sight.
namespace astro {

/// Geodetic sub-satellite point (the point on the ellipsoid along the local
/// normal through the satellite) of an ITRS position; h is the altitude.
[[nodiscard]] Geodetic subSatellitePoint(const Vec3& satItrs, const Ellipsoid& ell = Ellipsoid::wgs84()) noexcept;

/// Sub-satellite points at times start, start+step, ..., end.
[[nodiscard]] std::vector<Geodetic> groundTrack(const std::function<Vec3(const Time&)>& satItrs,
                                                const Time&                             start,
                                                const Time&                             end,
                                                double                                  step,
                                                const Ellipsoid&                        ell = Ellipsoid::wgs84());

/// Earth central angle (rad) from the sub-satellite point to the edge of
/// coverage for a satellite at `altitude` and a minimum elevation, spherical
/// Earth of radius `radius`.
[[nodiscard]] double
coverageHalfAngle(double altitude, double minElevation = 0.0, double radius = constants::wgs84A) noexcept;

/// Slant range (m) to the edge of coverage.
[[nodiscard]] double
maxSlantRange(double altitude, double minElevation = 0.0, double radius = constants::wgs84A) noexcept;

/// Coverage footprint boundary: `n` points at central angle `halfAngle` around
/// a sub-satellite point (spherical geometry; lat/lon in rad, h = 0).
[[nodiscard]] std::vector<Geodetic> footprint(double lat, double lon, double halfAngle, int n = 72);

/// True if the straight segment between two positions does not intersect the
/// ellipsoid (optionally inflated by `margin` meters, e.g. for the atmosphere).
[[nodiscard]] bool
lineOfSight(const Vec3& a, const Vec3& b, const Ellipsoid& ell = Ellipsoid::wgs84(), double margin = 0.0) noexcept;

} // namespace astro
