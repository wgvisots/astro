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
    // Points within ~43 km of the center exercise the bracketed root search inside the evolute.
    std::mt19937                           rng(3);
    std::uniform_real_distribution<double> u(-40000.0, 40000.0);
    for (int i = 0; i < 200; ++i) {
        const Vec3     p{u(rng), u(rng), u(rng)};
        const Geodetic g = ecefToGeodetic(p);
        expectVecNear(geodeticToEcef(g), p, 1e-6, "evolute");
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
