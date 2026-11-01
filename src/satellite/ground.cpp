/**********************************************************************
 * ground.cpp
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

#include "astro/satellite/ground.h"

#include <algorithm>
#include <cmath>

#include "astro/core/angle.h"

namespace astro {

Geodetic subSatellitePoint(const Vec3& satItrs, const Ellipsoid& ell) noexcept {
    return ecefToGeodetic(satItrs, ell);
}

std::vector<Geodetic> groundTrack(const std::function<Vec3(const Time&)>& satItrs, const Time& start, const Time& end,
                                  double step, const Ellipsoid& ell) {
    std::vector<Geodetic> out;
    const double span = end - start;
    if (step <= 0.0 || span < 0.0) return out;
    const auto n = static_cast<int>(std::floor(span / step + 1e-9));
    for (int i = 0; i <= n; ++i) out.push_back(subSatellitePoint(satItrs(start + i * step), ell));
    if (span - n * step > 1e-6) out.push_back(subSatellitePoint(satItrs(end), ell));
    return out;
}

double coverageHalfAngle(double altitude, double minElevation, double radius) noexcept {
    // Law of sines in the Earth-center / site / satellite triangle.
    const double eta = std::asin(std::clamp(radius * std::cos(minElevation) / (radius + altitude), -1.0, 1.0));
    return constants::halfPi - minElevation - eta;
}

double maxSlantRange(double altitude, double minElevation, double radius) noexcept {
    const double lambda = coverageHalfAngle(altitude, minElevation, radius);
    const double rs = radius + altitude;
    return std::sqrt(radius * radius + rs * rs - 2.0 * radius * rs * std::cos(lambda));
}

std::vector<Geodetic> footprint(double lat, double lon, double halfAngle, int n) {
    std::vector<Geodetic> out;
    const double sl = std::sin(lat), cl = std::cos(lat), sd = std::sin(halfAngle), cd = std::cos(halfAngle);
    for (int i = 0; i < n; ++i) {
        const double az = constants::twoPi * i / n;
        const double lat2 = std::asin(std::clamp(sl * cd + cl * sd * std::cos(az), -1.0, 1.0));
        const double lon2 = lon + std::atan2(std::sin(az) * sd * cl, cd - sl * std::sin(lat2));
        out.push_back({lat2, wrapPi(lon2), 0.0});
    }
    return out;
}

bool lineOfSight(const Vec3& a, const Vec3& b, const Ellipsoid& ell, double margin) noexcept {
    // Scale to a unit sphere and test the minimum distance of the segment from the origin.
    const double ra = ell.a + margin, rb = ell.b() + margin;
    const Vec3 p{a.x / ra, a.y / ra, a.z / rb};
    const Vec3 q{b.x / ra, b.y / ra, b.z / rb};
    const Vec3 d = q - p;
    const double dd = norm2(d);
    const double t = dd > 0.0 ? std::clamp(-dot(p, d) / dd, 0.0, 1.0) : 0.0;
    return norm(p + t * d) >= 1.0;
}

} // namespace astro
