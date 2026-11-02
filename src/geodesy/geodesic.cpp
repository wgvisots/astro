/**********************************************************************
 * geodesic.cpp
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

#include "astro/geodesy/geodesic.h"

#include <algorithm>
#include <cmath>

#include "astro/core/angle.h"
#include "astro/core/error.h"

namespace astro {

namespace {
struct AB {
    double A, B;
};
AB vincentyAB(double cos2Alpha, const Ellipsoid& ell) {
    const double a = ell.a, b = ell.b();
    const double u2 = cos2Alpha * (a * a - b * b) / (b * b);
    const double A  = 1.0 + u2 / 16384.0 * (4096.0 + u2 * (-768.0 + u2 * (320.0 - 175.0 * u2)));
    const double B  = u2 / 1024.0 * (256.0 + u2 * (-128.0 + u2 * (74.0 - 47.0 * u2)));
    return {A, B};
}
double deltaSigma(double B, double sinS, double cosS, double cos2Sm) {
    const double c2 = cos2Sm * cos2Sm;
    return B * sinS
         * (cos2Sm
            + B / 4.0 * (cosS * (-1.0 + 2.0 * c2) - B / 6.0 * cos2Sm * (-3.0 + 4.0 * sinS * sinS) * (-3.0 + 4.0 * c2)));
}
} // namespace

GeodesicInverseResult geodesicInverse(double lat1, double lon1, double lat2, double lon2, const Ellipsoid& ell) {
    const double f   = ell.f;
    const double L   = wrapPi(lon2 - lon1);
    const double U1  = std::atan((1.0 - f) * std::tan(lat1));
    const double U2  = std::atan((1.0 - f) * std::tan(lat2));
    const double sU1 = std::sin(U1), cU1 = std::cos(U1), sU2 = std::sin(U2), cU2 = std::cos(U2);

    double lambda = L, sinS = 0.0, cosS = 0.0, sigma = 0.0, cos2A = 0.0, cos2Sm = 0.0, sinA = 0.0;
    int    it = 0;
    for (; it < 200; ++it) {
        const double sl = std::sin(lambda), cl = std::cos(lambda);
        sinS = std::hypot(cU2 * sl, cU1 * sU2 - sU1 * cU2 * cl);
        if (sinS == 0.0) {
            return {0.0, 0.0, 0.0}; // coincident points
        }
        cosS              = sU1 * sU2 + cU1 * cU2 * cl;
        sigma             = std::atan2(sinS, cosS);
        sinA              = cU1 * cU2 * sl / sinS;
        cos2A             = 1.0 - sinA * sinA;
        cos2Sm            = (cos2A != 0.0) ? cosS - 2.0 * sU1 * sU2 / cos2A : 0.0;
        const double C    = f / 16.0 * cos2A * (4.0 + f * (4.0 - 3.0 * cos2A));
        const double prev = lambda;
        lambda = L + (1.0 - C) * f * sinA * (sigma + C * sinS * (cos2Sm + C * cosS * (-1.0 + 2.0 * cos2Sm * cos2Sm)));
        if (std::abs(lambda - prev) < 1e-13) {
            break;
        }
    }
    if (it >= 200) {
        throw Error("geodesicInverse: Vincenty iteration did not converge (nearly antipodal points)");
    }
    const auto [A, B] = vincentyAB(cos2A, ell);
    const double s    = ell.b() * A * (sigma - deltaSigma(B, sinS, cosS, cos2Sm));
    const double sl = std::sin(lambda), cl = std::cos(lambda);
    const double a1 = std::atan2(cU2 * sl, cU1 * sU2 - sU1 * cU2 * cl);
    const double a2 = std::atan2(cU1 * sl, -sU1 * cU2 + cU1 * sU2 * cl);
    return {s, wrapTwoPi(a1), wrapTwoPi(a2)};
}

GeodesicDirectResult geodesicDirect(double lat1, double lon1, double azimuth1, double distance, const Ellipsoid& ell) {
    const double f = ell.f, b = ell.b();
    const double sA1 = std::sin(azimuth1), cA1 = std::cos(azimuth1);
    const double U1  = std::atan((1.0 - f) * std::tan(lat1));
    const double sU1 = std::sin(U1), cU1 = std::cos(U1);
    const double sigma1 = std::atan2(std::tan(U1), cA1);
    const double sinA   = cU1 * sA1;
    const double cos2A  = 1.0 - sinA * sinA;
    const auto [A, B]   = vincentyAB(cos2A, ell);

    double sigma = distance / (b * A), sinS = 0.0, cosS = 0.0, cos2Sm = 0.0;
    for (int it = 0; it < 200; ++it) {
        cos2Sm            = std::cos(2.0 * sigma1 + sigma);
        sinS              = std::sin(sigma);
        cosS              = std::cos(sigma);
        const double prev = sigma;
        sigma             = distance / (b * A) + deltaSigma(B, sinS, cosS, cos2Sm);
        if (std::abs(sigma - prev) < 1e-14) {
            break;
        }
    }
    cos2Sm              = std::cos(2.0 * sigma1 + sigma);
    sinS                = std::sin(sigma);
    cosS                = std::cos(sigma);
    const double tmp    = sU1 * sinS - cU1 * cosS * cA1;
    const double lat2   = std::atan2(sU1 * cosS + cU1 * sinS * cA1, (1.0 - f) * std::hypot(sinA, tmp));
    const double lambda = std::atan2(sinS * sA1, cU1 * cosS - sU1 * sinS * cA1);
    const double C      = f / 16.0 * cos2A * (4.0 + f * (4.0 - 3.0 * cos2A));
    const double L =
        lambda - (1.0 - C) * f * sinA * (sigma + C * sinS * (cos2Sm + C * cosS * (-1.0 + 2.0 * cos2Sm * cos2Sm)));
    const double a2 = std::atan2(sinA, -tmp);
    return {lat2, wrapPi(lon1 + L), wrapTwoPi(a2)};
}

double haversineDistance(double lat1, double lon1, double lat2, double lon2, double r) noexcept {
    const double sdl = std::sin(0.5 * (lat2 - lat1));
    const double sdn = std::sin(0.5 * (lon2 - lon1));
    const double h   = sdl * sdl + std::cos(lat1) * std::cos(lat2) * sdn * sdn;
    return 2.0 * r * std::asin(std::min(1.0, std::sqrt(h)));
}

} // namespace astro
