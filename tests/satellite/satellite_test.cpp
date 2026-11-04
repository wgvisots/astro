/**********************************************************************
 * satellite_test.cpp
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

#include "astro/core/angle.h"
#include "astro/ephemeris/low_precision.h"
#include "astro/frames/transform.h"
#include "astro/geodesy/local.h"
#include "astro/orbit/elements.h"
#include "astro/orbit/two_body.h"
#include "astro/satellite/eclipse.h"
#include "astro/satellite/ground.h"
#include "astro/satellite/look.h"
#include "astro/satellite/orbital_frames.h"
#include "test_util.h"

using namespace astro;
using astro::test::expectMatNear;
using astro::test::expectOrthonormal;
using astro::test::expectVecNear;

namespace {
constexpr double km   = 1000.0;
const Tle        kIss = Tle::parse("1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927",
                            "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537");
} // namespace

TEST(Look, OverheadAndRangeRate) {
    const Geodetic site{deg2rad(30.0), deg2rad(-90.0), 0.0};
    const Vec3     up   = ecefToEnuMatrix(site.lat, site.lon).row(2);
    const Vec3     east = ecefToEnuMatrix(site.lat, site.lon).row(0);
    // Directly overhead, moving east: elevation 90, range rate 0.
    const StateVector sat{geodeticToEcef(site) + 500 * km * up, 7 * km * east};
    const LookAngles  la = lookAngles(sat, site);
    EXPECT_NEAR(la.el, constants::halfPi, 1e-9);
    EXPECT_NEAR(la.range, 500 * km, 1e-6);
    EXPECT_NEAR(la.rangeRate, 0.0, 1e-9);
    // Receding straight up at 1 km/s.
    EXPECT_NEAR(lookAngles({sat.r, 1 * km * up}, site).rangeRate, 1 * km, 1e-9);
    // Doppler: approaching raises the received frequency.
    EXPECT_NEAR(dopplerShift(437e6, -7000.0), 437e6 * 7000.0 / constants::c, 1e-9);
    // They differ by the second-order term f b^2 / 2 (~0.12 Hz here).
    const double b = 7000.0 / constants::c;
    EXPECT_NEAR(dopplerShiftRelativistic(437e6, -7000.0) - dopplerShift(437e6, -7000.0), 437e6 * b * b / 2.0, 1e-3);
}

TEST(Passes, IssOverSiteMatchesBruteForce) {
    const Sgp4        sgp4(kIss);
    const ConstantEop eop;
    const auto        sat = sgp4ItrsFunction(sgp4, eop);
    const Geodetic    site{deg2rad(40.0), deg2rad(-105.0), 1600.0};
    const Time        start  = kIss.epoch;
    const Time        end    = start + 86400.0;
    const double      mask   = deg2rad(10.0);
    const auto        passes = findPasses(sat, site, start, end, mask, 30.0);
    ASSERT_GE(passes.size(), 2u);
    for (const Pass& p : passes) {
        if (p.startsBeforeWindow || p.endsAfterWindow) {
            continue;
        }
        EXPECT_LT(p.aos, p.tca);
        EXPECT_LT(p.tca, p.los);
        EXPECT_LT(p.los - p.aos, 15 * 60.0); // ISS passes last < 15 minutes
        EXPECT_NEAR(lookAngles(sat(p.aos), site).el, mask, 1e-5);
        EXPECT_NEAR(lookAngles(sat(p.los), site).el, mask, 1e-5);
        // Brute-force maximum elevation: 1 s scan, then 10 ms around the best sample (high passes peak sharply).
        double best = -1.0, bestS = 0.0;
        for (double s = 0.0; s <= p.los - p.aos; s += 1.0) {
            const double el = lookAngles(sat(p.aos + s), site).el;
            if (el > best) {
                best  = el;
                bestS = s;
            }
        }
        for (double s = bestS - 1.0; s <= bestS + 1.0; s += 0.01) {
            best = std::max(best, lookAngles(sat(p.aos + s), site).el);
        }
        EXPECT_NEAR(p.maxElevation, best, deg2rad(0.001));
        EXPECT_GE(p.maxElevation, best - 1e-9);
    }
}

TEST(Eclipse, ShadowGeometry) {
    const Vec3   sun{constants::au, 0, 0};
    const double r = 7000 * km;
    // Sunward: lit. Directly behind: umbra.
    EXPECT_EQ(shadow({r, 0, 0}, sun).state, ShadowState::Sunlit);
    EXPECT_EQ(shadow({-r, 0, 0}, sun).state, ShadowState::Umbra);
    EXPECT_EQ(shadow({-r, 0, 0}, sun, ShadowModel::Cylindrical).state, ShadowState::Umbra);
    // Illumination is continuous and monotone across the penumbra.
    double prev = 0.0;
    for (double y = 6000 * km; y < 6800 * km; y += 2 * km) {
        const ShadowResult s = shadow({-r, y, 0}, sun);
        EXPECT_GE(s.illumination, prev - 1e-12);
        prev = s.illumination;
    }
    EXPECT_NEAR(prev, 1.0, 1e-12);
    // Penumbra width near LEO is ~2 * r_behind * (solar angular radius) ~ 60 km.
    int penumbra = 0;
    for (double y = 6000 * km; y < 6800 * km; y += 1 * km) {
        penumbra += shadow({-r, y, 0}, sun).state == ShadowState::Penumbra;
    }
    EXPECT_GT(penumbra, 40);
    EXPECT_LT(penumbra, 90);
    // Far beyond the umbra cone (~1.4e6 km), the Earth only partially covers the Sun: annular.
    EXPECT_EQ(shadow({-3.0e9, 0, 0}, sun).state, ShadowState::Annular);
    // Beta angle.
    const StateVector polar{
        {r, 0, 0       },
        {0, 0, 7.5 * km}
    }; // orbit normal along -y
    EXPECT_NEAR(betaAngle(polar, {0, -constants::au, 0}), constants::halfPi, 1e-12);
    EXPECT_NEAR(betaAngle(polar, sun), 0.0, 1e-12);
}

TEST(Eclipse, CylindricalDurationMatchesAnalytic) {
    const double      r      = constants::wgs84A + 500 * km;
    const StateVector s0     = elementsToState({r, 0.0, 0.0, 0.0, 0.0, 0.0});
    const double      period = orbitalPeriod(r);
    const Time        t0;
    auto              sat = [&](const Time& t) { return propagateKepler(s0, t - t0).r; };
    auto              sun = [](const Time&) { return Vec3{constants::au, 0, 0}; };
    const auto        ecl = findEclipses(sat, sun, t0, t0 + 2.0 * period, 20.0, ShadowModel::Cylindrical);
    ASSERT_EQ(ecl.size(), 2u);
    const double expected = 2.0 * std::asin(constants::wgs84A / r) / constants::twoPi * period;
    EXPECT_NEAR(ecl[0].penumbraExit - ecl[0].penumbraEntry, expected, 0.01);
    // The conical model adds short penumbra phases around a slightly shorter umbra.
    const auto con = findEclipses(sat, sun, t0, t0 + 2.0 * period, 20.0, ShadowModel::Conical);
    ASSERT_EQ(con.size(), 2u);
    EXPECT_TRUE(con[0].umbral);
    EXPECT_LT(con[0].umbraExit - con[0].umbraEntry, expected);
    EXPECT_GT(con[0].penumbraExit - con[0].penumbraEntry, expected);
    EXPECT_NEAR(con[0].umbraEntry - con[0].penumbraEntry, 8.0, 3.0); // ~8 s penumbra at 500 km
}

TEST(Eclipse, IssWithRealSun) {
    const Sgp4                  sgp4(kIss);
    const ConstantEop           eop;
    const LowPrecisionEphemeris eph;
    auto                        sat = [&](const Time& t) {
        return FrameTransformer(t, eop, {.model = PrecessionNutationModel::IAU2006_2000B})
            .transform(sgp4.propagate(t), Frame::TEME, Frame::GCRS)
            .r;
    };
    auto       sun   = [&](const Time& t) { return sunPositionLowPrecision(t); };
    const Time start = kIss.epoch;
    const auto ecl   = findEclipses(sat, sun, start, start + 86400.0);
    // ~15.7 orbits/day; one eclipse per orbit unless beta is large.
    EXPECT_GE(ecl.size(), 14u);
    EXPECT_LE(ecl.size(), 17u);
    for (const auto& e : ecl) {
        if (e.penumbraEntry == start || e.penumbraExit == start + 86400.0) {
            continue;
        }
        const double dur = e.penumbraExit - e.penumbraEntry;
        EXPECT_GT(dur, 20 * 60.0);
        EXPECT_LT(dur, 40 * 60.0);
        EXPECT_EQ(shadow(sat(e.umbraEntry + 0.5 * (e.umbraExit - e.umbraEntry)),
                         sun(e.umbraEntry + 0.5 * (e.umbraExit - e.umbraEntry)))
                      .state,
                  ShadowState::Umbra);
    }
}

TEST(Ground, SubpointFootprintAndLineOfSight) {
    // GEO over lon 0.
    const Geodetic g = subSatellitePoint({42164.17 * km, 0, 0});
    EXPECT_NEAR(g.lat, 0.0, 1e-12);
    EXPECT_NEAR(g.lon, 0.0, 1e-12);
    EXPECT_NEAR(g.h, 42164.17 * km - constants::wgs84A, 1e-6);
    // GEO coverage to the horizon: acos(R / r) ~ 81.3 deg.
    EXPECT_NEAR(rad2deg(coverageHalfAngle(42164.17 * km - constants::wgs84A)), 81.30, 0.01);
    EXPECT_LT(coverageHalfAngle(500 * km, deg2rad(10.0)), coverageHalfAngle(500 * km));
    EXPECT_NEAR(maxSlantRange(500 * km) / km,
                std::sqrt(std::pow(constants::wgs84A + 500 * km, 2) - std::pow(constants::wgs84A, 2)) / km,
                1e-6);
    // Footprint points are all at the requested central angle.
    const double lam = deg2rad(20.0);
    for (const Geodetic& p : footprint(deg2rad(45.0), deg2rad(10.0), lam, 36)) {
        const double c = std::acos(std::sin(deg2rad(45.0)) * std::sin(p.lat)
                                   + std::cos(deg2rad(45.0)) * std::cos(p.lat) * std::cos(p.lon - deg2rad(10.0)));
        EXPECT_NEAR(c, lam, 1e-12);
    }
    // Line of sight.
    EXPECT_FALSE(lineOfSight({8000 * km, 0, 0}, {-8000 * km, 0, 0}));
    EXPECT_TRUE(lineOfSight({8000 * km, 0, 0}, {8000 * km, 1000 * km, 0}));
    EXPECT_TRUE(lineOfSight({0, 0, 7000 * km}, {7000 * km, 0, 7000 * km}));
    EXPECT_FALSE(lineOfSight({7000 * km, 0, 0}, {0, 7000 * km, 0})); // chord dips to 4950 km
    EXPECT_FALSE(lineOfSight({6400 * km, 0, 0}, {0, 6400 * km, 0}, Ellipsoid::wgs84(), 100 * km));
    // Ground track sampling.
    auto geo = [](const Time&) { return Vec3{42164.17 * km, 0, 0}; };
    EXPECT_EQ(groundTrack(geo, Time(), Time() + 3600.0, 600.0).size(), 7u);
}

TEST(OrbitalFrames, DefinitionsAndRelativeMotion) {
    const StateVector s = elementsToState({7000 * km, 0.1, 0.8, 0.3, 0.5, 1.0});
    for (const Mat3& m : {ricMatrix(s), lvlhMatrix(s), vncMatrix(s), ntwMatrix(s)}) {
        expectOrthonormal(m, 1e-14);
    }
    expectVecNear(ricMatrix(s) * s.r, {norm(s.r), 0, 0}, 1e-6);
    expectVecNear(lvlhMatrix(s) * s.r, {0, 0, -norm(s.r)}, 1e-6);
    expectVecNear(vncMatrix(s) * s.v, {norm(s.v), 0, 0}, 1e-9);
    expectVecNear(ntwMatrix(s) * s.v, {0, norm(s.v), 0}, 1e-9);
    // RIC and LVLH differ by a fixed permutation/sign change.
    const Mat3 P = lvlhMatrix(s) * ricMatrix(s).transpose();
    expectMatNear(P, Mat3{0, 1, 0, 0, 0, -1, -1, 0, 0}, 1e-14);
    // A deputy slightly ahead in the same circular orbit: pure in-track offset, ~zero relative velocity.
    const double      r      = 7000 * km;
    const StateVector chief  = elementsToState({r, 0.0, 0.5, 0.0, 0.0, 0.0});
    const StateVector deputy = elementsToState({r, 0.0, 0.5, 0.0, 0.0, 1e-4});
    const StateVector rel    = relativeStateRic(chief, deputy);
    EXPECT_NEAR(rel.r.y, r * 1e-4, 1.0);
    EXPECT_NEAR(rel.r.x, 0.0, 0.1);
    EXPECT_NEAR(norm(rel.v), 0.0, 1e-3);
    const StateVector back = deputyFromRelativeRic(chief, rel);
    expectVecNear(back.r, deputy.r, 1e-6);
    expectVecNear(back.v, deputy.v, 1e-9);
}

TEST(OrbitalFrames, CovarianceRotation) {
    Mat6         P{};
    const double sig[6] = {100.0, 50.0, 10.0, 0.1, 0.05, 0.01};
    for (std::size_t i = 0; i < 6; ++i) {
        P[i][i] = sig[i] * sig[i];
    }
    P[0][1] = P[1][0]     = 0.3 * 100.0 * 50.0;
    const StateVector s   = elementsToState({7000 * km, 0.01, 1.0, 0.2, 0.3, 0.4});
    const Mat3        R   = ricMatrix(s);
    const Mat6        Q   = rotateCovariance(P, R);
    double            trP = 0, trQ = 0, trPv = 0, trQv = 0;
    for (std::size_t i = 0; i < 3; ++i) {
        trP  += P[i][i];
        trQ  += Q[i][i];
        trPv += P[i + 3][i + 3];
        trQv += Q[i + 3][i + 3];
    }
    EXPECT_NEAR(trQ, trP, 1e-8);
    EXPECT_NEAR(trQv, trPv, 1e-12);
    for (std::size_t i = 0; i < 6; ++i) {
        for (std::size_t j = 0; j < 6; ++j) {
            EXPECT_NEAR(Q[i][j], Q[j][i], 1e-9);
        }
    }
    const Mat6 back = rotateCovariance(Q, R.transpose());
    for (std::size_t i = 0; i < 6; ++i) {
        for (std::size_t j = 0; j < 6; ++j) {
            EXPECT_NEAR(back[i][j], P[i][j], 1e-8);
        }
    }
}
