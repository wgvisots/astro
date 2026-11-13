/**********************************************************************
 * geodetic.cpp
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

#include "astro/geodesy/geodetic.h"

#include <cmath>
#include <limits>

#include "astro/core/rootfind.h"

namespace astro {

Vec3 geodeticToEcef(const Geodetic& g, const Ellipsoid& ell) noexcept {
    const double sl = std::sin(g.lat), cl = std::cos(g.lat);
    const double N   = ell.a / std::sqrt(1.0 - ell.e2() * sl * sl);
    const double rho = (N + g.h) * cl;
    return {rho * std::cos(g.lon), rho * std::sin(g.lon), (N * (1.0 - ell.e2()) + g.h) * sl};
}

namespace {

// Height above the ellipsoid for a given geodetic latitude (robust at all latitudes).
double heightFor(double rho, double z, double lat, const Ellipsoid& ell) {
    const double s = std::sin(lat), c = std::cos(lat);
    return rho * c + z * s - ell.a * std::sqrt(1.0 - ell.e2() * s * s);
}

// Inside the evolute up to four normals to the ellipse pass through the point.
// Find every root of  f(phi) = rho sin(phi) - z cos(phi) - N(phi) e^2 sin(phi) cos(phi)
// in [0, pi/2] (first quadrant, rho > 0, z >= 0) and keep the nearest foot point (largest h).
double latitudeNearCenter(double rho, double z, const Ellipsoid& ell) {
    const double e2 = ell.e2();
    auto         f  = [&](double phi) {
        const double s = std::sin(phi), c = std::cos(phi);
        return rho * s - z * c - ell.a / std::sqrt(1.0 - e2 * s * s) * e2 * s * c;
    };
    constexpr int n       = 512;
    double        bestPhi = 0.0, bestH = -std::numeric_limits<double>::infinity();
    auto          consider = [&](double phi) {
        const double h = heightFor(rho, z, phi, ell);
        if (h > bestH) {
            bestH   = h;
            bestPhi = phi;
        }
    };
    double prevPhi = 0.0, prevF = f(0.0);
    if (prevF == 0.0) {
        consider(0.0);
    }
    for (int i = 1; i <= n; ++i) {
        const double phi = constants::halfPi * i / n;
        const double fi  = f(phi);
        if (fi == 0.0) {
            consider(phi);
        } else if (prevF != 0.0 && (fi > 0.0) != (prevF > 0.0)) {
            consider(brentRoot(f, prevPhi, phi, 1e-16));
        }
        prevPhi = phi;
        prevF   = fi;
    }
    // The scan ends at the double nearest pi/2, which is 6e-17 rad short of the pole. At the pole
    // itself f is exactly rho, and rho > 0 here. So if f is still negative at the last sample, a foot
    // point lies in that gap (the point is within about 5e-12 m of the polar axis): it is nearer the
    // pole than a double can express, and the last sample is the latitude.
    if (prevF < 0.0) {
        consider(constants::halfPi);
    }
    return bestPhi;
}

} // namespace

Geodetic ecefToGeodetic(const Vec3& r, const Ellipsoid& ell) noexcept {
    const double x = r.x, y = r.y, z = r.z;
    const double rho = std::hypot(x, y);
    Geodetic     g;
    // Longitude: half-angle form is accurate everywhere except the -x axis.
    g.lon = (rho == 0.0) ? 0.0 : (x >= 0.0 ? 2.0 * std::atan2(y, rho + x) : std::atan2(y, x));
    if (ell.f == 0.0) {
        g.lat = std::atan2(z, rho);
        g.h   = std::hypot(rho, z) - ell.a;
        return g;
    }
    const double a = ell.a, e2 = ell.e2(), e4 = e2 * e2;
    if (rho == 0.0) { // on the polar axis
        g.lat = z >= 0.0 ? constants::halfPi : -constants::halfPi;
        g.h   = std::abs(z) - ell.b();
        return g;
    }
    // Vermeille (2002), "Direct transformation from geocentric coordinates to geodetic coordinates",
    // J. Geodesy 76, 451-454. (The longitude above is the half-angle form of Vermeille 2004,
    // "Computing geodetic coordinates from geocentric coordinates", J. Geodesy 78.)
    const double p  = (rho * rho) / (a * a);
    const double q  = (1.0 - e2) * z * z / (a * a);
    const double rr = (p + q - e4) / 6.0;
    if (rr > 0.0) {
        const double s      = e4 * p * q / (4.0 * rr * rr * rr);
        const double t      = std::cbrt(1.0 + s + std::sqrt(s * (2.0 + s)));
        const double u      = rr * (1.0 + t + 1.0 / t);
        const double v      = std::sqrt(u * u + e4 * q);
        const double w      = e2 * (u + v - q) / (2.0 * v);
        const double k      = std::sqrt(u + v + w * w) - w;
        const double D      = k * rho / (k + e2);
        const double sqrtDZ = std::hypot(D, z);
        g.lat               = 2.0 * std::atan2(z, D + sqrtDZ);
        g.h                 = (k + e2 - 1.0) / k * sqrtDZ;
        return g;
    }
    // Inside the evolute (within ~43 km of the center for WGS84).
    const double phi = latitudeNearCenter(rho, std::abs(z), ell);
    g.lat            = z < 0.0 ? -phi : phi;
    g.h              = heightFor(rho, z, g.lat, ell);
    return g;
}

double geodeticToGeocentricLatitude(double lat, const Ellipsoid& ell) noexcept {
    return std::atan2((1.0 - ell.e2()) * std::sin(lat), std::cos(lat));
}

double geocentricToGeodeticLatitude(double lat, const Ellipsoid& ell) noexcept {
    return std::atan2(std::sin(lat), (1.0 - ell.e2()) * std::cos(lat));
}

double geodeticToReducedLatitude(double lat, const Ellipsoid& ell) noexcept {
    return std::atan2((1.0 - ell.f) * std::sin(lat), std::cos(lat));
}

double reducedToGeodeticLatitude(double beta, const Ellipsoid& ell) noexcept {
    return std::atan2(std::sin(beta), (1.0 - ell.f) * std::cos(beta));
}

Spherical cartesianToSpherical(const Vec3& v) noexcept {
    const double rho = std::hypot(v.x, v.y);
    return {std::atan2(v.z, rho), (rho == 0.0) ? 0.0 : std::atan2(v.y, v.x), std::hypot(rho, v.z)};
}

Vec3 sphericalToCartesian(const Spherical& s) noexcept {
    const double cl = std::cos(s.lat);
    return {s.r * cl * std::cos(s.lon), s.r * cl * std::sin(s.lon), s.r * std::sin(s.lat)};
}

} // namespace astro
