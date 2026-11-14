/**********************************************************************
 * geodesy_test.cpp
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

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <random>

#include "astro/core/angle.h"
#include "astro/core/error.h"
#include "astro/geodesy/geodesic.h"
#include "astro/geodesy/geodetic.h"
#include "astro/geodesy/helmert.h"
#include "astro/geodesy/local.h"
#include "test_util.h"

using namespace astro;
using astro::test::dataPath;
using astro::test::expectMatNear;
using astro::test::expectOrthonormal;
using astro::test::expectVecNear;
using astro::test::readCsv;

TEST(Ellipsoid, DerivedQuantities) {
    const Ellipsoid w = Ellipsoid::wgs84();
    EXPECT_NEAR(w.b(), 6356752.314245, 1e-6);
    EXPECT_NEAR(w.e2(), 6.69437999014e-3, 1e-14);
    EXPECT_NEAR(w.primeVerticalRadius(0.0), w.a, 1e-9);
    EXPECT_NEAR(w.meridionalRadius(constants::halfPi), w.a * w.a / w.b(), 1e-6);
    EXPECT_NEAR(w.geocentricRadius(0.0), w.a, 1e-9);
    EXPECT_NEAR(w.geocentricRadius(constants::halfPi), w.b(), 1e-6);
}

TEST(Ellipsoid, PublishedAndClosedFormValues) {
    // WGS84 derived constants published in NIMA TR8350.2 (3rd edition, 2000), table 3.3, to their
    // printed precision: second eccentricity squared 6.73949674228e-3, mean radius R1 6371008.7714 m.
    const Ellipsoid w = Ellipsoid::wgs84();
    EXPECT_NEAR(w.ep2(), 6.73949674228e-3, 0.5e-14);
    EXPECT_NEAR(w.meanRadius(), 6371008.7714, 0.5e-4);
    // Third flattening n = (a - b) / (a + b), from the axes rather than from f. a - b = 21 km loses the
    // rounding of b (1e-9 m), 4e-14 of it, so 1e-16 on n (measured 1.6e-17).
    EXPECT_NEAR(w.n(), (w.a - w.b()) / (w.a + w.b()), 1e-16);
    // The Gaussian radius sqrt(M N) is a sqrt(1 - e^2) = b on the equator (M = a (1 - e^2), N = a) and
    // a^2 / b at the poles (M = N); between, it lies between the two radii it is the mean of.
    EXPECT_NEAR(w.gaussianRadius(0.0), w.b(), 1e-8); // a few units in the last place of 6.4e6 m
    EXPECT_NEAR(w.gaussianRadius(constants::halfPi), w.a * w.a / w.b(), 1e-8);
    const double lat = 0.8;
    EXPECT_GT(w.gaussianRadius(lat), w.meridionalRadius(lat));
    EXPECT_LT(w.gaussianRadius(lat), w.primeVerticalRadius(lat));
    // IERS Conventions (2010), table 1.1: a = 6378136.6 m, 1/f = 298.25642.
    const Ellipsoid iers = Ellipsoid::iers2010();
    EXPECT_EQ(iers.a, 6378136.6);
    EXPECT_NEAR(1.0 / iers.f, 298.25642, 1e-9);
    // GRS80 (Moritz 2000): 1/f = 298.257222101.
    EXPECT_NEAR(1.0 / Ellipsoid::grs80().f, 298.257222101, 1e-9);
}

TEST(Ellipsoid, SphereHasEveryRadiusEqualToItsOwn) {
    const double    r = 6371000.0;
    const Ellipsoid s = Ellipsoid::sphere(r);
    EXPECT_EQ(s.a, r);
    EXPECT_EQ(s.f, 0.0);
    EXPECT_EQ(s.b(), r);
    EXPECT_EQ(s.e2(), 0.0);
    EXPECT_EQ(s.ep2(), 0.0);
    EXPECT_EQ(s.n(), 0.0);
    EXPECT_EQ(s.meanRadius(), r);
    for (double lat : {0.0, 0.7, constants::halfPi}) {
        EXPECT_EQ(s.primeVerticalRadius(lat), r);
        EXPECT_EQ(s.meridionalRadius(lat), r);
        EXPECT_EQ(s.gaussianRadius(lat), r);
        EXPECT_NEAR(s.geocentricRadius(lat), r, 1e-9);
    }
}

TEST(Geodetic, OnASphereLatitudeIsGeocentricAndHeightIsRadialDistance) {
    // With no flattening the normal through a point passes through the centre: the geodetic latitude is
    // atan2(z, rho) and the height |r| - R, for points outside, on and inside the sphere. This is the
    // f = 0 branch of ecefToGeodetic().
    const double    r   = 6371000.0;
    const Ellipsoid sph = Ellipsoid::sphere(r);
    for (const Vec3& p : {
             Vec3{7.0e6,  -2.0e6, 3.0e6 },
             Vec3{-1.0e6, 1.0e5,  -6.0e6},
             Vec3{3.0e5,  4.0e5,  1.0e5 }
    }) {
        const Geodetic g = ecefToGeodetic(p, sph);
        EXPECT_NEAR(g.lat, std::atan2(p.z, std::hypot(p.x, p.y)), 1e-15);
        EXPECT_NEAR(g.lon, std::atan2(p.y, p.x), 1e-15);
        EXPECT_NEAR(g.h, norm(p) - r, 1e-8); // rounding of terms of size r
        expectVecNear(geodeticToEcef(g, sph), p, 1e-8, "sphere round trip");
    }
}

TEST(Geodetic, KnownPoint) {
    // Equator / prime meridian / poles.
    expectVecNear(geodeticToEcef({0, 0, 0}), {constants::wgs84A, 0, 0}, 1e-9);
    expectVecNear(geodeticToEcef({constants::halfPi, 0, 100}), {0, 0, Ellipsoid::wgs84().b() + 100}, 1e-6);
    const Geodetic g = ecefToGeodetic({0, 0, -Ellipsoid::wgs84().b() - 5.0});
    EXPECT_NEAR(g.lat, -constants::halfPi, 1e-15);
    EXPECT_NEAR(g.h, 5.0, 1e-8);
    const Geodetic c = ecefToGeodetic({0, 0, 0});
    EXPECT_NEAR(c.h, -Ellipsoid::wgs84().b(), 1e-6);
}

TEST(Geodetic, MatchesErfaAndRoundTrips) {
    const auto csv = readCsv(dataPath("geodetic.csv"));
    for (const auto& r : csv.rows) {
        const double lat = std::stod(r[0]), lon = std::stod(r[1]), h = std::stod(r[2]);
        const Vec3   xyz{std::stod(r[3]), std::stod(r[4]), std::stod(r[5])};
        const double tol = 1e-8 + 1e-15 * norm(xyz); // a few ulp
        expectVecNear(geodeticToEcef({lat, lon, h}), xyz, tol, "gd->ecef");
        const Geodetic g = ecefToGeodetic(xyz);
        // The generating (lat, h) is the unique answer outside the evolute; ERFA (Fukushima 2006) itself
        // degrades to ~1e-7 rad deep inside the Earth and ~1e-11 rad at GEO, so compare to it only near the
        // surface.
        EXPECT_NEAR(g.lat, lat, 2e-15) << "lat, h=" << h;
        EXPECT_NEAR(wrapPi(g.lon - lon), 0.0, 1e-15) << "lon";
        EXPECT_NEAR(g.h, h, 1e-8 * std::max(1.0, std::abs(h) / 1e6)) << "h, h=" << h;
        if (std::abs(h) < 1e5) {
            // Measured: 3.3e-16 rad and 2.5e-9 m from ERFA's gc2gd. Tolerances: three units in the last
            // place of the latitude, and the library's own 1e-8 m bound for the height.
            EXPECT_NEAR(g.lat, std::stod(r[6]), 1e-15) << "lat vs ERFA, h=" << h;
            EXPECT_NEAR(g.h, std::stod(r[8]), 1e-8) << "h vs ERFA, h=" << h;
        }
        expectVecNear(geodeticToEcef(g), xyz, tol, "round trip");
    }
}

namespace {

/// Geodetic colatitude (pi/2 - latitude) and height of a point of a meridian plane.
struct MeridianReference {
    long double colatitude; ///< rad
    long double height;     ///< m
};

/// Reference for ecefToGeodetic() at the point (rho, z), rho >= 0 and z >= 0, of a meridian plane:
/// extended precision, and a different route from the library's. With g the parametric (reduced)
/// colatitude, the distance from (rho, z) to the point (a sin g, b cos g) of the ellipse is stationary
/// where
///     a rho cos g - b z sin g - (a^2 - b^2) sin g cos g = 0.
/// Measuring g from the pole makes the pole exactly g = 0, so a root 1e-18 rad from it is found as
/// surely as any other. Every root in the quadrant is bracketed on a fine scan and refined by
/// bisection; the pole is one more when rho = 0, and the equator when z = 0. The nearest of these
/// points is the answer, and tan(colatitude) = (b / a) tan g.
MeridianReference meridianReference(long double rho, long double z, const Ellipsoid& ell) {
    const long double a      = ell.a;
    const long double b      = a * (1.0L - static_cast<long double>(ell.f));
    const long double halfPi = std::numbers::pi_v<long double> / 2.0L;
    const auto        slope  = [&](long double g) {
        return a * rho * std::cos(g) - b * z * std::sin(g) - (a * a - b * b) * std::sin(g) * std::cos(g);
    };
    const auto distance2 = [&](long double g) {
        const long double dRho = rho - a * std::sin(g), dZ = z - b * std::cos(g);
        return dRho * dRho + dZ * dZ;
    };
    long double best = 0.0L, bestDistance2 = std::numeric_limits<long double>::infinity();
    const auto  consider = [&](long double g) {
        if (distance2(g) < bestDistance2) {
            bestDistance2 = distance2(g);
            best          = g;
        }
    };
    if (rho == 0.0L) {
        consider(0.0L);
    }
    if (z == 0.0L) {
        consider(halfPi);
    }
    constexpr int n = 4096;
    for (int i = 0; i < n; ++i) {
        long double lo = halfPi * i / n, hi = halfPi * (i + 1) / n;
        if ((slope(lo) > 0.0L) == (slope(hi) > 0.0L)) {
            continue;
        }
        for (int it = 0; it < 100; ++it) {
            const long double mid                                 = 0.5L * (lo + hi);
            ((slope(mid) > 0.0L) == (slope(lo) > 0.0L) ? lo : hi) = mid;
        }
        consider(0.5L * (lo + hi));
    }
    const bool inside = (rho / a) * (rho / a) + (z / b) * (z / b) < 1.0L;
    return {std::atan2(b * std::sin(best), a * std::cos(best)), (inside ? -1.0L : 1.0L) * std::sqrt(bestDistance2)};
}

} // namespace

TEST(Geodetic, OnAndNextToThePolarAxisMatchesExtendedPrecision) {
    // Distances from the polar axis from zero up to a metre, at heights inside the evolute (it reaches
    // 42.84 km along the axis for WGS84), just outside it, and at the surface. Next to the axis the
    // nearest point of the ellipsoid is next to the pole: the latitude must approach +-90 deg, and the
    // height |z| - b, continuously as the distance goes to zero.
    const Ellipsoid   ell    = Ellipsoid::wgs84();
    const long double halfPi = std::numbers::pi_v<long double> / 2.0L;
    const double      lon    = 2.0;
    for (double z : {0.0, 1.0, 1000.0, 20000.0, 42000.0, 42800.0, 42900.0, 1.0e5, 6.4e6}) {
        for (double rho : {0.0, 1e-300, 1e-13, 1e-12, 3e-12, 1e-11, 1e-9, 1e-3, 1.0}) {
            const MeridianReference ref = meridianReference(rho, z, ell);
            for (double sign : {1.0, -1.0}) {
                const Geodetic g = ecefToGeodetic({rho * std::cos(lon), rho * std::sin(lon), sign * z}, ell);
                SCOPED_TRACE(testing::Message() << "rho = " << rho << " m, z = " << sign * z << " m");
                // 1e-15 rad: a few units in the last place of pi/2 (2.2e-16).
                const long double latError = std::abs(static_cast<long double>(g.lat)) - (halfPi - ref.colatitude);
                EXPECT_NEAR(static_cast<double>(latError), 0.0, 1e-15);
                if (z > 0.0) {
                    EXPECT_EQ(std::signbit(g.lat), sign < 0.0);
                }
                // 1e-8 m: the library's own bound for the height (rounding of terms of size a).
                EXPECT_NEAR(static_cast<double>(static_cast<long double>(g.h) - ref.height), 0.0, 1e-8);
                // The longitude is that of the point, however close to the axis; on the axis it is 0.
                EXPECT_NEAR(g.lon, rho > 0.0 ? lon : 0.0, 1e-15);
            }
        }
    }
    // The reviewer's point: 1e-13 m from the axis, 1 km north of the center.
    const Geodetic g = ecefToGeodetic({1e-13, 0.0, 1000.0}, ell);
    EXPECT_NEAR(g.lat, constants::halfPi, 1e-15);
    EXPECT_NEAR(g.h, 1000.0 - ell.b(), 1e-8);
}

TEST(Geodetic, InsideTheEvoluteMatchesExtendedPrecision) {
    // Within about 43 km of the centre (inside the evolute of the meridian ellipse) a point has several
    // normals to the ellipsoid and the library falls back to a bracketed root search; the answer is the
    // nearest foot of a normal, found by meridianReference() in long double.
    const Ellipsoid                        ell = Ellipsoid::wgs84();
    std::mt19937                           rng(3);
    std::uniform_real_distribution<double> u(-40000.0, 40000.0);
    for (int i = 0; i < 40; ++i) {
        const Vec3              p{u(rng), u(rng), u(rng)};
        const MeridianReference ref = meridianReference(std::hypot(p.x, p.y), std::abs(p.z), ell);
        const Geodetic          g   = ecefToGeodetic(p, ell);
        SCOPED_TRACE(testing::Message() << "point " << p.x << ", " << p.y << ", " << p.z);
        const double lat = static_cast<double>(std::numbers::pi_v<long double> / 2.0L - ref.colatitude);
        EXPECT_NEAR(g.lat, p.z < 0.0 ? -lat : lat, 1e-15); // a few units in the last place of pi/2
        EXPECT_NEAR(g.h, static_cast<double>(ref.height), 1e-8);
        EXPECT_NEAR(g.lon, std::atan2(p.y, p.x), 1e-15);
        expectVecNear(geodeticToEcef(g, ell), p, 1e-8, "round trip");
    }
}

TEST(Geodetic, Latitudes) {
    const double phi = deg2rad(45.0);
    const double psi = geodeticToGeocentricLatitude(phi);
    EXPECT_NEAR(rad2deg(phi - psi) * 3600.0, 692.7, 0.1); // ~11.5 arcmin at 45 deg
    EXPECT_NEAR(geocentricToGeodeticLatitude(psi), phi, 1e-15);
    EXPECT_NEAR(reducedToGeodeticLatitude(geodeticToReducedLatitude(phi)), phi, 1e-15);
    // The definitions, on WGS84 and on a strongly flattened ellipsoid: tan(geocentric) = (1 - f)^2
    // tan(geodetic) and tan(reduced) = (1 - f) tan(geodetic). The point (a cos(beta), b sin(beta)) of
    // the meridian ellipse has geocentric latitude atan2(b sin(beta), a cos(beta)).
    for (const Ellipsoid& ell : {
             Ellipsoid::wgs84(),
             Ellipsoid{1.0, 0.3}
    }) {
        const double k = 1.0 - ell.f;
        for (double lat : {-1.2, -0.3, 0.4, 1.0, 1.5}) {
            SCOPED_TRACE(testing::Message() << "f = " << ell.f << ", lat = " << lat);
            const double beta = std::atan(k * std::tan(lat));
            EXPECT_NEAR(geodeticToGeocentricLatitude(lat, ell), std::atan(k * k * std::tan(lat)), 1e-15);
            EXPECT_NEAR(geodeticToReducedLatitude(lat, ell), beta, 1e-15);
            EXPECT_NEAR(reducedToGeodeticLatitude(beta, ell), lat, 1e-15);
            EXPECT_NEAR(geocentricToGeodeticLatitude(std::atan2(ell.b() * std::sin(beta), ell.a * std::cos(beta)), ell),
                        lat,
                        1e-15);
        }
    }
    const Spherical s = cartesianToSpherical({1, 1, std::sqrt(2.0)});
    EXPECT_NEAR(s.lat, deg2rad(45.0), 1e-15);
    EXPECT_NEAR(s.lon, deg2rad(45.0), 1e-15);
    EXPECT_NEAR(s.r, 2.0, 1e-15);
    expectVecNear(sphericalToCartesian(s), {1, 1, std::sqrt(2.0)}, 1e-15);
}

TEST(LocalFrames, EnuNedAer) {
    // "Up" and "north" from geometry written here, not from the matrices under test: up is the outward
    // normal of the ellipsoid (x^2 + y^2) / a^2 + z^2 / b^2 = 1 at the origin's foot, the gradient
    // (x / a^2, y / a^2, z / b^2); north is the direction of the polar axis within the horizontal plane;
    // east = north x up. Both WGS84 and a more flattened ellipsoid (1/f = 150).
    for (const Ellipsoid& ell : {
             Ellipsoid::wgs84(),
             Ellipsoid{6.4e6, 1.0 / 150.0}
    }) {
        SCOPED_TRACE(testing::Message() << "f = " << ell.f);
        const Geodetic o{deg2rad(40.0), deg2rad(-105.0), 1600.0};
        const Vec3     foot  = geodeticToEcef({o.lat, o.lon, 0.0}, ell);
        const Vec3     base  = geodeticToEcef(o, ell);
        const double   b     = ell.b();
        const Vec3     up    = normalized(Vec3{foot.x / (ell.a * ell.a), foot.y / (ell.a * ell.a), foot.z / (b * b)});
        const Vec3     north = normalized(Vec3{0, 0, 1} - up.z * up);
        const Vec3     east  = cross(north, up);
        // 1e-15: rounding in unit vectors (the normal is computed from coordinates of 6e6 m).
        expectMatNear(ecefToEnuMatrix(o.lat, o.lon), Mat3::fromRows(east, north, up), 1e-15, "ENU rows");
        expectMatNear(ecefToNedMatrix(o.lat, o.lon), Mat3::fromRows(north, east, -up), 1e-15, "NED rows");
        expectVecNear(base, foot + o.h * up, 1e-8, "height along the normal");

        // An arbitrary offset, its components along the three axes, and its look angles.
        const Vec3 d{-3456.0, 2345.0, 1234.0};
        const Vec3 enu{dot(d, east), dot(d, north), dot(d, up)};
        // 1e-8 m: rounding of coordinates of 6e6 m (1e-9 m) subtracted from one another.
        expectVecNear(ecefToEnu(base + d, o, ell), enu, 1e-8, "ENU");
        expectVecNear(ecefToNed(base + d, o, ell), {enu.y, enu.x, -enu.z}, 1e-8, "NED");
        expectVecNear(enuToEcef(enu, o, ell), base + d, 1e-8, "ENU to ECEF");
        expectVecNear(nedToEcef({enu.y, enu.x, -enu.z}, o, ell), base + d, 1e-8, "NED to ECEF");
        const Aer a = ecefToAer(base + d, o, ell);
        EXPECT_NEAR(a.az, wrapTwoPi(std::atan2(enu.x, enu.y)), 1e-12); // 1e-8 m seen from 4 km
        EXPECT_NEAR(a.el, std::asin(enu.z / norm(d)), 1e-12);
        EXPECT_NEAR(a.range, norm(d), 1e-8);
        expectVecNear(aerToEcef(a, o, ell), base + d, 1e-8, "AER to ECEF");

        // Due north, east, south and west on the horizon, and straight up.
        const double r = 5000.0;
        EXPECT_NEAR(wrapPi(ecefToAer(base + r * north, o, ell).az), 0.0, 1e-12);
        EXPECT_NEAR(ecefToAer(base + r * east, o, ell).az, constants::halfPi, 1e-12);
        EXPECT_NEAR(ecefToAer(base - r * north, o, ell).az, constants::pi, 1e-12);
        EXPECT_NEAR(ecefToAer(base - r * east, o, ell).az, 1.5 * constants::pi, 1e-12);
        EXPECT_NEAR(ecefToAer(base + r * north, o, ell).el, 0.0, 1e-12);
        EXPECT_NEAR(ecefToAer(base + r * up, o, ell).el, constants::halfPi, 1e-12);
    }
    // enuToAer and aerToEnu with the angles written out.
    EXPECT_NEAR(enuToAer({10, 0, 0}).az, constants::halfPi, 1e-15);
    EXPECT_NEAR(enuToAer({-10, 0, 0}).az, 1.5 * constants::pi, 1e-15);
    const Aer x{deg2rad(123.0), deg2rad(17.0), 1.2e6};
    expectVecNear(aerToEnu(x),
                  {x.range * std::cos(x.el) * std::sin(x.az),
                   x.range * std::cos(x.el) * std::cos(x.az),
                   x.range * std::sin(x.el)},
                  1e-9,
                  "AER to ENU");
}

TEST(LocalFrames, AtAPoleNorthPointsAlongTheOriginsMeridian) {
    // At the north pole every direction is south; the ENU axes there are the limit of those of a point
    // just short of the pole on the meridian `lon`: north is -(cos lon, sin lon, 0), east is
    // (-sin lon, cos lon, 0), up is z. At the south pole north is +(cos lon, sin lon, 0) and up is -z.
    const double lon = 0.7;
    const Vec3   east{-std::sin(lon), std::cos(lon), 0.0};
    const Vec3   alongMeridian{std::cos(lon), std::sin(lon), 0.0};
    // 1e-15: cos(pi/2) is 6e-17, not 0.
    expectMatNear(ecefToEnuMatrix(constants::halfPi, lon), Mat3::fromRows(east, -alongMeridian, {0, 0, 1}), 1e-15);
    expectMatNear(ecefToEnuMatrix(-constants::halfPi, lon), Mat3::fromRows(east, alongMeridian, {0, 0, -1}), 1e-15);
}

TEST(Geodesic, VincentyFlindersPeakToBuninyong) {
    // The Flinders Peak to Buninyong sample line of the Geocentric Datum of Australia technical
    // manual (Intergovernmental Committee on Surveying and Mapping), on GRS80: 54 972.271 m,
    // forward azimuth 306 deg 52' 05.37". Vincenty's formulae are the method under test.
    const Ellipsoid grs80 = Ellipsoid::grs80();
    const double    lat1 = -dmsToRad(37, 57, 3.72030), lon1 = dmsToRad(144, 25, 29.52440);
    const double    lat2 = -dmsToRad(37, 39, 10.15610), lon2 = dmsToRad(143, 55, 35.38390);
    const auto      inv = geodesicInverse(lat1, lon1, lat2, lon2, grs80);
    EXPECT_NEAR(inv.distance, 54972.271, 1e-3);
    EXPECT_NEAR(rad2arcsec(inv.azimuth1 - dmsToRad(306, 52, 5.37)), 0.0, 0.01);
    EXPECT_NEAR(rad2arcsec(inv.azimuth2 - dmsToRad(307, 10, 25.07)), 0.0, 0.01);
    const auto dir = geodesicDirect(lat1, lon1, inv.azimuth1, inv.distance, grs80);
    EXPECT_NEAR(dir.lat, lat2, 1e-11);
    EXPECT_NEAR(dir.lon, lon2, 1e-11);
    EXPECT_NEAR(dir.azimuth2, inv.azimuth2, 1e-10);
}

namespace {

/// End of a geodesic found by integrating its differential equation (geodesicByIntegration()).
struct IntegratedGeodesic {
    long double lat;        ///< geodetic latitude, rad
    long double lon;        ///< longitude, rad
    long double azimuth;    ///< forward azimuth, rad, clockwise from north
    long double clairaut;   ///< largest change along the line of Clairaut's constant (rho sin(azimuth)), m
    long double offSurface; ///< largest |(x^2 + y^2) / a^2 + z^2 / b^2 - 1| along the line
};

/// Reference for geodesicDirect() and geodesicInverse() that shares nothing with Vincenty's series: the
/// geodesic as a curve r(s) of the surface F(r) = (x^2 + y^2) / a^2 + z^2 / b^2 - 1 = 0 whose
/// acceleration is along the surface normal (it does not turn within the surface),
///     r'' = -lambda g,   g = (x / a^2, y / a^2, z / b^2),
///     lambda = ((v_x^2 + v_y^2) / a^2 + v_z^2 / b^2) / |g|^2,
/// where lambda keeps r on the surface (the second derivative of F along the curve is zero). Integrated
/// in Cartesian coordinates, which have no singularity at the poles, with the classical fourth-order
/// Runge-Kutta method in long double (64-bit significand) and `steps` equal steps. Clairaut's
/// constant (the distance from the axis times the sine of the azimuth) and the distance from the surface
/// are tracked as checks of the integration itself.
IntegratedGeodesic
geodesicByIntegration(double lat1, double lon1, double azimuth1, double distance, const Ellipsoid& ell, int steps) {
    using LD    = long double;
    const LD a  = ell.a;
    const LD b  = a * (1.0L - static_cast<LD>(ell.f));
    const LD a2 = a * a, b2 = b * b;
    struct Point {
        LD x, y, z, vx, vy, vz;
    };
    const auto derivative = [&](const Point& p) {
        const LD gx = p.x / a2, gy = p.y / a2, gz = p.z / b2;
        const LD q      = (p.vx * p.vx + p.vy * p.vy) / a2 + p.vz * p.vz / b2;
        const LD lambda = q / (gx * gx + gy * gy + gz * gz);
        return Point{p.vx, p.vy, p.vz, -lambda * gx, -lambda * gy, -lambda * gz};
    };
    const auto axpy = [](const Point& p, LD h, const Point& d) {
        return Point{p.x + h * d.x, p.y + h * d.y, p.z + h * d.z, p.vx + h * d.vx, p.vy + h * d.vy, p.vz + h * d.vz};
    };
    // Unit vectors east and north of the surface normal at geodetic latitude `lat`, longitude `lon`.
    const auto horizon = [](LD lat, LD lon, LD(&east)[3], LD(&north)[3]) {
        east[0]  = -std::sin(lon);
        east[1]  = std::cos(lon);
        east[2]  = 0.0L;
        north[0] = -std::sin(lat) * std::cos(lon);
        north[1] = -std::sin(lat) * std::sin(lon);
        north[2] = std::cos(lat);
    };
    // Geodetic latitude, longitude and azimuth of a point of the surface: the latitude is that of the
    // normal g.
    const auto describe = [&](const Point& p, LD& lat, LD& lon, LD& az) {
        lat = std::atan2(p.z / b2, std::hypot(p.x, p.y) / a2);
        lon = std::atan2(p.y, p.x);
        LD east[3], north[3];
        horizon(lat, lon, east, north);
        az = std::atan2(p.vx * east[0] + p.vy * east[1] + p.vz * east[2],
                        p.vx * north[0] + p.vy * north[1] + p.vz * north[2]);
    };
    const LD phi = lat1, lam = lon1, alpha = azimuth1;
    const LD e2 = 1.0L - b2 / a2;
    const LD N  = a / std::sqrt(1.0L - e2 * std::sin(phi) * std::sin(phi));
    LD       east[3], north[3];
    horizon(phi, lam, east, north);
    Point p{N * std::cos(phi) * std::cos(lam),
            N * std::cos(phi) * std::sin(lam),
            N * (1.0L - e2) * std::sin(phi),
            std::cos(alpha) * north[0] + std::sin(alpha) * east[0],
            std::cos(alpha) * north[1] + std::sin(alpha) * east[1],
            std::cos(alpha) * north[2] + std::sin(alpha) * east[2]};
    // Clairaut's constant as the z component of r x v (rho sin(azimuth), conserved by the axial symmetry).
    const auto clairaut   = [](const Point& q) { return q.x * q.vy - q.y * q.vx; };
    const LD   clairaut0  = clairaut(p);
    LD         clairautCh = 0.0L, offSurface = 0.0L;
    const LD   h = static_cast<LD>(distance) / steps;
    for (int i = 0; i < steps; ++i) {
        const Point k1  = derivative(p);
        const Point k2  = derivative(axpy(p, h / 2, k1));
        const Point k3  = derivative(axpy(p, h / 2, k2));
        const Point k4  = derivative(axpy(p, h, k3));
        p.x            += h / 6 * (k1.x + 2 * k2.x + 2 * k3.x + k4.x);
        p.y            += h / 6 * (k1.y + 2 * k2.y + 2 * k3.y + k4.y);
        p.z            += h / 6 * (k1.z + 2 * k2.z + 2 * k3.z + k4.z);
        p.vx           += h / 6 * (k1.vx + 2 * k2.vx + 2 * k3.vx + k4.vx);
        p.vy           += h / 6 * (k1.vy + 2 * k2.vy + 2 * k3.vy + k4.vy);
        p.vz           += h / 6 * (k1.vz + 2 * k2.vz + 2 * k3.vz + k4.vz);
        clairautCh      = std::max(clairautCh, std::abs(clairaut(p) - clairaut0));
        offSurface      = std::max(offSurface, std::abs((p.x * p.x + p.y * p.y) / a2 + p.z * p.z / b2 - 1.0L));
    }
    IntegratedGeodesic out{};
    describe(p, out.lat, out.lon, out.azimuth);
    out.clairaut   = clairautCh;
    out.offSurface = offSurface;
    return out;
}

/// A test line: start, azimuth and length.
struct GeodesicLine {
    const char* what;
    double      lat1, lon1, azimuth1, distance; // rad, rad, rad, m
};

// clang-format off
const GeodesicLine kGeodesicLines[] = {
    {"1 km, mid-latitude",                 0.7,   0.3,   0.5,                          1.0e3 },
    {"100 km, mid-latitude",               0.7,   0.3,   0.5,                          1.0e5 },
    {"1000 km, southern, westward",       -0.6,   2.9,   4.0,                          1.0e6 },
    {"5000 km, across the antimeridian",   0.2,   3.0,   1.2,                          5.0e6 },
    {"10000 km, high latitude",            1.2,  -1.0,   2.5,                          1.0e7 },
    {"15000 km",                          -0.3,   0.1,   0.9,                          1.5e7 },
    {"19000 km",                           0.4,  -2.0,   2.0,                          1.9e7 },
    {"near-meridional, over the pole",     0.5,   1.0,   1e-3,                         1.0e7 },
    {"near-meridional, southward",        -0.1,  -0.5,   constants::pi + 1e-4,         1.2e7 },
    {"near-equatorial, eastward",          1e-4,  0.0,   constants::halfPi - 1e-4,     1.0e7 },
    {"near-equatorial, westward",         -1e-4,  1.0,   1.5 * constants::pi + 1e-3,   1.7e7 },
};
// clang-format on

} // namespace

namespace {

/// Distance (m) on the surface between two nearby points given by geodetic latitude and longitude:
/// the meridional radius times the latitude difference and the radius of the parallel times the
/// longitude difference.
double surfaceOffset(double lat, double dLat, double dLon, const Ellipsoid& ell) {
    return std::hypot(ell.meridionalRadius(lat) * dLat, ell.primeVerticalRadius(lat) * std::cos(lat) * dLon);
}

int stepsFor(double distance) {
    // 2 km steps: halving them changes the end point by less than 3e-16 rad (measured), so the reference
    // is a thousand times more accurate than the tolerances below.
    return std::max(200, static_cast<int>(distance / 2000.0));
}

// Vincenty's series are truncated in the flattening; the header promises sub-millimetre results for
// the Earth. Measured against the integration: end points within 0.085 mm (the 19000 km line),
// distances within 0.071 mm, azimuths within 3.3e-11 rad. Tolerances: 0.1 mm and 5e-11 rad.
constexpr double kVincentyLength  = 1e-4;
constexpr double kVincentyAzimuth = 5e-11;

} // namespace

TEST(Geodesic, DirectMatchesTheIntegratedGeodesicEquation) {
    // geodesicDirect() against the curve found by integrating the geodesic equation (see
    // geodesicByIntegration()): short and long lines, near-meridional and near-equatorial ones.
    const Ellipsoid ell = Ellipsoid::wgs84();
    for (const GeodesicLine& line : kGeodesicLines) {
        SCOPED_TRACE(line.what);
        const IntegratedGeodesic ref =
            geodesicByIntegration(line.lat1, line.lon1, line.azimuth1, line.distance, ell, stepsFor(line.distance));
        // Checks of the integration: the true curve keeps Clairaut's constant (measured change 1.7e-11 m
        // of rounding, on 6e6 m) and stays on the surface (measured 1.7e-17).
        EXPECT_LT(static_cast<double>(ref.clairaut), 1e-10);
        EXPECT_LT(static_cast<double>(ref.offSurface), 1e-16);
        const GeodesicDirectResult dir  = geodesicDirect(line.lat1, line.lon1, line.azimuth1, line.distance, ell);
        const double               dLat = static_cast<double>(dir.lat - ref.lat);
        const double               dLon = wrapPi(static_cast<double>(dir.lon - ref.lon));
        EXPECT_LT(surfaceOffset(dir.lat, dLat, dLon, ell), kVincentyLength);
        EXPECT_NEAR(wrapPi(static_cast<double>(dir.azimuth2 - ref.azimuth)), 0.0, kVincentyAzimuth);
    }
}

TEST(Geodesic, InverseMatchesTheIntegratedGeodesicEquation) {
    // geodesicInverse() between the ends of the integrated curves must give back their length and the
    // azimuths at both ends. Every line here is the shortest between its ends (all are shorter than half
    // the meridian's or the equator's length) and none is nearly antipodal, so the iteration converges.
    const Ellipsoid ell = Ellipsoid::wgs84();
    for (const GeodesicLine& line : kGeodesicLines) {
        SCOPED_TRACE(line.what);
        const IntegratedGeodesic ref =
            geodesicByIntegration(line.lat1, line.lon1, line.azimuth1, line.distance, ell, stepsFor(line.distance));
        const GeodesicInverseResult inv =
            geodesicInverse(line.lat1, line.lon1, static_cast<double>(ref.lat), static_cast<double>(ref.lon), ell);
        EXPECT_NEAR(inv.distance, line.distance, kVincentyLength);
        EXPECT_NEAR(wrapPi(inv.azimuth1 - line.azimuth1), 0.0, kVincentyAzimuth);
        EXPECT_NEAR(wrapPi(static_cast<double>(inv.azimuth2 - ref.azimuth)), 0.0, kVincentyAzimuth);
        EXPECT_GE(inv.azimuth1, 0.0); // [0, 2 pi), as documented
        EXPECT_LT(inv.azimuth1, constants::twoPi);
    }
}

TEST(Geodesic, OnASphereTheGeodesicIsTheGreatCircle) {
    // With no flattening the geodesic is a great circle, whose spherical trigonometry is written out
    // here: distance r * sigma with cos(sigma) from the law of cosines, azimuth from the four-part
    // formula, and the destination from the sine and cosine rules.
    const double    r    = 6371000.0;
    const Ellipsoid sph  = Ellipsoid::sphere(r);
    const double    lat1 = 0.4, lon1 = -1.1, lat2 = -0.7, lon2 = 1.9;
    const double    dLon = lon2 - lon1;
    const double sigma = std::acos(std::sin(lat1) * std::sin(lat2) + std::cos(lat1) * std::cos(lat2) * std::cos(dLon));
    const double az1   = std::atan2(std::sin(dLon) * std::cos(lat2),
                                  std::cos(lat1) * std::sin(lat2) - std::sin(lat1) * std::cos(lat2) * std::cos(dLon));
    const GeodesicInverseResult inv = geodesicInverse(lat1, lon1, lat2, lon2, sph);
    // The arc is 2.8 rad, near half a circle, where acos magnifies the rounding of its argument by
    // 1 / sin(sigma) = 3 and the haversine's asin by 1 / cos(sigma / 2) = 6: 1.3e-15 rad, times r, is
    // 8e-9 m. Tolerance 5e-8 m (measured 1.1e-8 m for the haversine).
    EXPECT_NEAR(inv.distance, r * sigma, 5e-8);
    EXPECT_NEAR(wrapPi(inv.azimuth1 - az1), 0.0, 1e-14);
    EXPECT_NEAR(haversineDistance(lat1, lon1, lat2, lon2, r), r * sigma, 5e-8);

    const double s = 9.0e6, az = 2.2, c = s / r;
    const double latEnd = std::asin(std::sin(lat1) * std::cos(c) + std::cos(lat1) * std::sin(c) * std::cos(az));
    const double lonEnd =
        lon1 + std::atan2(std::sin(az) * std::sin(c) * std::cos(lat1), std::cos(c) - std::sin(lat1) * std::sin(latEnd));
    const GeodesicDirectResult dir = geodesicDirect(lat1, lon1, az, s, sph);
    EXPECT_NEAR(dir.lat, latEnd, 1e-14);
    EXPECT_NEAR(wrapPi(dir.lon - lonEnd), 0.0, 1e-14);
}

TEST(Geodesic, CoincidentPointsAndHaversine) {
    EXPECT_EQ(geodesicInverse(0.3, 0.4, 0.3, 0.4).distance, 0.0);
    // Half a great circle, and a quarter of one along the equator.
    EXPECT_NEAR(haversineDistance(0, 0, 0, constants::pi, 1.0), constants::pi, 1e-15);
    EXPECT_NEAR(haversineDistance(0, 0, 0, constants::halfPi, 2.0), constants::pi, 1e-15);
}

TEST(Geodesic, NearlyAntipodalOrNotFiniteInverseThrows) {
    // On the equator Vincenty's iteration does not converge within about 0.6 deg of the antipode
    // (from 179.40 deg of longitude); the header promises astro::Error, and the same for input that
    // is not finite.
    EXPECT_THROW((void)geodesicInverse(0.0, 0.0, 0.0, deg2rad(179.7)), Error);
    EXPECT_THROW((void)geodesicInverse(0.0, 0.0, 0.0, std::numeric_limits<double>::quiet_NaN()), Error);
    EXPECT_THROW((void)geodesicInverse(std::numeric_limits<double>::infinity(), 0.0, 0.1, 0.2), Error);
    EXPECT_NO_THROW((void)geodesicInverse(0.0, 0.0, 0.0, deg2rad(179.0)));
}

TEST(Helmert, TransformAndInverse) {
    // Pure translation + scale + small rotation about z.
    HelmertParams p{
        {1.0, -2.0, 3.0              },
        {0.0, 0.0,  arcsec2rad(0.001)},
        1e-9
    };
    const Vec3   x{4.0e6, 3.0e6, 3.5e6};
    const Vec3   y  = helmertTransform(x, p);
    const double rz = arcsec2rad(0.001);
    expectVecNear(y,
                  {x.x + 1.0 + 1e-9 * x.x - rz * x.y, x.y - 2.0 + 1e-9 * x.y + rz * x.x, x.z + 3.0 + 1e-9 * x.z},
                  1e-9);
    // inverse() negates the parameters, which undoes the transformation to first order; what is left
    // is second order (s t, r t, s^2 x, r^2 x, r s x), about 2e-8 m here (measured 1.1e-8 m).
    expectVecNear(helmertTransform(y, inverse(p)), x, 5e-8);
    // Coordinate-frame convention flips rotation sign only.
    HelmertParams cf = p;
    cf.convention    = HelmertConvention::CoordinateFrame;
    expectVecNear(helmertTransform(x, cf),
                  {x.x + 1.0 + 1e-9 * x.x + rz * x.y, x.y - 2.0 + 1e-9 * x.y - rz * x.x, x.z + 3.0 + 1e-9 * x.z},
                  1e-9);
    // Time-dependent parameters.
    HelmertTimeDependent h{
        p,
        {0.001, 0, 0},
        {},
        0.0,
        2010.0
    };
    EXPECT_NEAR(h.at(2020.0).t.x, 1.01, 1e-12);
    const StateVector sv = helmertTransform(
        StateVector{
            x,
            {0.01, 0, 0}
    },
        h,
        2020.0);
    // A velocity of 0.01 m/s gains the translation rate, 0.001 m per Julian year, and the scale term.
    EXPECT_NEAR(sv.v.x, 0.01 + 0.001 / (365.25 * 86400.0) + 1e-9 * 0.01, 1e-16);
}

TEST(Helmert, VelocityIsTheTimeDerivativeOfThePositionInMetresPerSecond) {
    // A station moving at centimetres per year, and a 14-parameter transformation with every parameter
    // and every rate non-zero. Positions are in m, velocities in m/s, rates per Julian year.
    constexpr double           secondsPerYear = 365.25 * 86400.0;
    const HelmertTimeDependent h{
        {{0.1, -0.2, 0.3}, {mas2rad(1.0), mas2rad(-2.0), mas2rad(3.0)}, 5e-9        },
        {0.01,             -0.02,                                       0.03        },
        {mas2rad(0.1),     mas2rad(-0.2),                               mas2rad(0.3)},
        1e-10,
        2010.0
    };
    const StateVector x{
        {4.0e6, 3.0e6, 3.5e6},
        Vec3{0.02,  -0.01, 0.015}
        / secondsPerYear
    };
    const double year = 2020.0;
    // Truth, written out from the definition X2 = X1 + T + s X1 + r x X1 (position-vector convention):
    // the transformed position `t` seconds after `year`, with the station moving and the parameters
    // changing at their rates.
    const auto position = [&](double t) {
        const double years = (year - h.epoch) + t / secondsPerYear;
        const Vec3   x1    = x.r + x.v * t;
        const Vec3   T     = h.p.t + h.tRate * years;
        const Vec3   r     = h.p.r + h.rRate * years;
        const double s     = h.p.s + h.sRate * years;
        return x1 + T + s * x1 + cross(r, x1);
    };
    const StateVector y = helmertTransform(x, h, year);
    expectVecNear(y.r, position(0.0), 5e-9, "position"); // a few units in the last place of 4e6 m
    // The position is a quadratic in t, so the central difference over +-1 year is its derivative
    // exactly, but for rounding: positions of 5e6 m carry 1e-9 m, which over 6.3e7 s is 2e-17 m/s.
    // The velocity itself is 7e-10 m/s, and the smallest rate term (the scale rate) 1e-11 m/s.
    const Vec3 velocity = (position(secondsPerYear) - position(-secondsPerYear)) / (2.0 * secondsPerYear);
    expectVecNear(y.v, velocity, 1e-16, "velocity");

    // The coordinate-frame convention with the opposite rotation signs is the same transformation.
    HelmertTimeDependent cf = h;
    cf.p.convention         = HelmertConvention::CoordinateFrame;
    cf.p.r                  = -h.p.r;
    cf.rRate                = -h.rRate;
    EXPECT_EQ(helmertTransform(x, cf, year).r, y.r);
    EXPECT_EQ(helmertTransform(x, cf, year).v, y.v);

    // Without rates the position is the 7-parameter result. EPSG Guidance Note 7-2, position vector
    // transformation, WGS 72 to WGS 84: tZ = 4.5 m, rZ = 0.554 arcsec, dS = 0.219 ppm. The result is
    // published to 1 cm.
    const HelmertTimeDependent wgs72{
        {{0.0, 0.0, 4.5}, {0.0, 0.0, arcsec2rad(0.554)}, 0.219e-6},
        {},
        {},
        0.0,
        1984.0
    };
    const StateVector w = helmertTransform(
        StateVector{
            {3657660.66, 255768.55, 5201382.11},
            {}
    },
        wgs72,
        2000.0);
    expectVecNear(w.r, {3657660.78, 255778.43, 5201387.75}, 0.01, "EPSG");
    EXPECT_EQ(w.v, Vec3{});
}
