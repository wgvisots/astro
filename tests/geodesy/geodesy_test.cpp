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
            EXPECT_NEAR(g.lat, std::stod(r[6]), 2e-12) << "lat vs ERFA, h=" << h;
            EXPECT_NEAR(g.h, std::stod(r[8]), 1e-5) << "h vs ERFA, h=" << h;
        }
        expectVecNear(geodeticToEcef(g), xyz, tol, "round trip");
    }
}

TEST(Geodetic, InsideEvoluteFallback) {
    // Points within ~43 km of the center exercise the Newton branch.
    std::mt19937                           rng(3);
    std::uniform_real_distribution<double> u(-40000.0, 40000.0);
    for (int i = 0; i < 200; ++i) {
        const Vec3     p{u(rng), u(rng), u(rng)};
        const Geodetic g = ecefToGeodetic(p);
        expectVecNear(geodeticToEcef(g), p, 1e-6, "evolute");
    }
}

TEST(Geodetic, Latitudes) {
    const double phi = deg2rad(45.0);
    const double psi = geodeticToGeocentricLatitude(phi);
    EXPECT_NEAR(rad2deg(phi - psi) * 3600.0, 692.7, 0.1); // ~11.5 arcmin at 45 deg
    EXPECT_NEAR(geocentricToGeodeticLatitude(psi), phi, 1e-15);
    EXPECT_NEAR(reducedToGeodeticLatitude(geodeticToReducedLatitude(phi)), phi, 1e-15);
    const Spherical s = cartesianToSpherical({1, 1, std::sqrt(2.0)});
    EXPECT_NEAR(s.lat, deg2rad(45.0), 1e-15);
    EXPECT_NEAR(s.lon, deg2rad(45.0), 1e-15);
    EXPECT_NEAR(s.r, 2.0, 1e-15);
    expectVecNear(sphericalToCartesian(s), {1, 1, std::sqrt(2.0)}, 1e-15);
}

TEST(LocalFrames, EnuNedAer) {
    const Geodetic o{deg2rad(40.0), deg2rad(-105.0), 1600.0};
    expectOrthonormal(ecefToEnuMatrix(o.lat, o.lon));
    expectOrthonormal(ecefToNedMatrix(o.lat, o.lon));
    const Vec3 base = geodeticToEcef(o);
    // Up is along the ellipsoid normal.
    const Vec3 up = ecefToEnuMatrix(o.lat, o.lon).row(2);
    expectVecNear(ecefToEnu(base + 1000.0 * up, o), {0, 0, 1000}, 1e-6);
    expectVecNear(ecefToNed(base + 1000.0 * up, o), {0, 0, -1000}, 1e-6);
    // A point due north at the horizon has azimuth 0.
    const Vec3 north = ecefToEnuMatrix(o.lat, o.lon).row(1);
    const Aer  a     = ecefToAer(base + 5000.0 * north, o);
    EXPECT_NEAR(wrapPi(a.az), 0.0, 1e-12);
    EXPECT_NEAR(a.el, 0.0, 1e-12);
    EXPECT_NEAR(a.range, 5000.0, 1e-6);
    // East is 90 deg; round trips.
    EXPECT_NEAR(enuToAer({10, 0, 0}).az, constants::halfPi, 1e-15);
    EXPECT_NEAR(enuToAer({-10, 0, 0}).az, 1.5 * constants::pi, 1e-15);
    const Aer x{deg2rad(123.0), deg2rad(17.0), 1.2e6};
    const Aer y = ecefToAer(aerToEcef(x, o), o);
    EXPECT_NEAR(y.az, x.az, 1e-12);
    EXPECT_NEAR(y.el, x.el, 1e-12);
    EXPECT_NEAR(y.range, x.range, 1e-6);
    const Vec3 p{1234.0, -567.0, 89.0};
    expectVecNear(nedToEcef(ecefToNed(enuToEcef(p, o), o), o), enuToEcef(p, o), 1e-6);
}

TEST(Geodesic, VincentyFlindersPeakToBuninyong) {
    // Vincenty (1975) test case on the GRS80 / ANS-like ellipsoid; published values for GRS80.
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

TEST(Geodesic, RandomRoundTrips) {
    std::mt19937                           rng(11);
    std::uniform_real_distribution<double> la(-1.5, 1.5), lo(-3.14, 3.14), az(0.0, 6.28), d(1.0, 1.5e7);
    for (int i = 0; i < 300; ++i) {
        const double          lat1 = la(rng), lon1 = lo(rng), a1 = az(rng), s = d(rng);
        const auto            dir = geodesicDirect(lat1, lon1, a1, s);
        GeodesicInverseResult inv;
        try {
            inv = geodesicInverse(lat1, lon1, dir.lat, dir.lon);
        } catch (const Error&) {
            continue; // nearly antipodal
        }
        EXPECT_NEAR(inv.distance, s, 1e-4 + 1e-12 * s);
        EXPECT_NEAR(wrapPi(inv.azimuth1 - a1), 0.0, 1e-9);
    }
    EXPECT_NEAR(haversineDistance(0, 0, 0, constants::pi, 1.0), constants::pi, 1e-15);
    EXPECT_DOUBLE_EQ(geodesicInverse(0.3, 0.4, 0.3, 0.4).distance, 0.0);
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
    expectVecNear(helmertTransform(y, inverse(p)), x, 1e-4);
    // Coordinate-frame convention flips rotation sign only.
    HelmertParams cf = p;
    cf.convention    = HelmertConvention::CoordinateFrame;
    const Vec3 yc    = helmertTransform(x, cf);
    EXPECT_NEAR(yc.x, x.x + 1.0 + 1e-9 * x.x + rz * x.y, 1e-9);
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
    EXPECT_NEAR(sv.v.x, 0.011 + 1e-9 * 0.01, 1e-12);
}
