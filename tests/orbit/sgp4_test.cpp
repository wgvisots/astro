/**********************************************************************
 * sgp4_test.cpp
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

#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>

#include "astro/core/angle.h"
#include "astro/orbit/sgp4.h"
#include "astro/orbit/tle.h"
#include "test_util.h"

using namespace astro;
using astro::test::dataPath;

TEST(Tle, ParseIssExample) {
    const Tle t = Tle::parse("1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927",
                             "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537", "ISS (ZARYA)");
    EXPECT_EQ(t.satnum, 25544);
    EXPECT_EQ(t.name, "ISS (ZARYA)");
    EXPECT_EQ(t.intlDesignator, "98067A");
    EXPECT_EQ(t.epoch.iso8601(3), "2008-09-20T12:25:40.104Z");
    EXPECT_NEAR(t.meanMotionDot, -0.00002182, 1e-15);
    EXPECT_NEAR(t.bstar, -0.11606e-4, 1e-15);
    EXPECT_NEAR(rad2deg(t.inclination), 51.6416, 1e-12);
    EXPECT_NEAR(t.eccentricity, 0.0006703, 1e-15);
    EXPECT_NEAR(t.meanMotion, 15.72125391, 1e-12);
    EXPECT_EQ(t.revolutionNumber, 56353);
    EXPECT_EQ(t.elementSetNumber, 292);
    // Formatting reproduces the lines exactly.
    EXPECT_EQ(t.line1(), "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927");
    EXPECT_EQ(t.line2(), "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537");
}

TEST(Tle, ChecksumAlpha5AndErrors) {
    EXPECT_EQ(tleChecksum("1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  292"), 7);
    EXPECT_EQ(decodeAlpha5("A0001"), 100001);
    EXPECT_EQ(decodeAlpha5("Z9999"), 339999);
    EXPECT_EQ(decodeAlpha5("J0000"), 180000);  // I is skipped
    EXPECT_EQ(encodeAlpha5(100001), "A0001");
    EXPECT_EQ(encodeAlpha5(339999), "Z9999");
    EXPECT_EQ(encodeAlpha5(25544), "25544");
    EXPECT_THROW((void)Tle::parse("1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2928",
                            "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537"),
                 Error);
    EXPECT_NO_THROW((void)Tle::parse("1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2928",
                               "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537", "", false));
    const auto many = Tle::parseMany(
        "ISS (ZARYA)\n"
        "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927\n"
        "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537\n"
        "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927\n"
        "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537\n");
    ASSERT_EQ(many.size(), 2u);
    EXPECT_EQ(many[0].name, "ISS (ZARYA)");
    EXPECT_EQ(many[1].name, "");
}

namespace {
struct VerCase {
    std::string line1, line2;
    double start{0}, stop{0}, step{0};
};
struct VerRow {
    double t;
    Vec3 r, v;  // km, km/s
};
} // namespace

TEST(Sgp4, VerificationAgainstVallado) {
    // SGP4-VER.TLE and tcppver.out from Vallado et al. (2006), as distributed with python-sgp4.
    std::ifstream tf(dataPath("SGP4-VER.TLE"));
    ASSERT_TRUE(tf);
    std::vector<VerCase> cases;
    std::string line, pending;
    while (std::getline(tf, line)) {
        if (line.empty() || line[0] == '#') continue;
        if (line[0] == '1') {
            pending = line;
        } else if (line[0] == '2') {
            VerCase c{pending, line.substr(0, 69)};
            std::istringstream extra(line.size() > 69 ? line.substr(69) : "");
            extra >> c.start >> c.stop >> c.step;
            cases.push_back(c);
        }
    }
    ASSERT_GE(cases.size(), 30u);

    std::ifstream of(dataPath("tcppver.out"));
    ASSERT_TRUE(of);
    std::vector<std::pair<int, std::vector<VerRow>>> expected;
    while (std::getline(of, line)) {
        if (line.find("xx") != std::string::npos) {
            expected.push_back({std::stoi(line), {}});
            continue;
        }
        std::istringstream ls(line);
        VerRow r;
        if (ls >> r.t >> r.r.x >> r.r.y >> r.r.z >> r.v.x >> r.v.y >> r.v.z) expected.back().second.push_back(r);
    }
    ASSERT_EQ(expected.size(), cases.size());

    int compared = 0;
    double maxPosErr = 0.0, maxVelErr = 0.0;
    for (std::size_t k = 0; k < cases.size(); ++k) {
        const VerCase& c = cases[k];
        const Tle tle = Tle::parse(c.line1, c.line2, "", false);
        ASSERT_EQ(tle.satnum, expected[k].first);
        const auto& rows = expected[k].second;
        // 33334 is Vallado's deliberate error case: SGP4 fails at epoch (perturbed eccentricity out of range).
        // Its single t=0 row in tcppver.out is an artifact of the test driver printing that line unconditionally;
        // the reference implementation itself returns error 3 there.
        if (rows.empty() || tle.satnum == 33334) {
            try {
                Sgp4 bad{tle};
                ADD_FAILURE() << tle.satnum << " should fail at epoch";
            } catch (const Sgp4Error& e) {
                EXPECT_EQ(e.code(), 3);
            }
            continue;
        }
        const Sgp4 sgp4(tle);
        for (const VerRow& row : rows) {
            StateVector s;
            const int code = sgp4.propagate(row.t, s);
            EXPECT_EQ(code, 0) << tle.satnum << " t=" << row.t;
            const double dr = norm(s.r / 1000.0 - row.r);
            const double dv = norm(s.v / 1000.0 - row.v);
            maxPosErr = std::max(maxPosErr, dr);
            maxVelErr = std::max(maxVelErr, dv);
            EXPECT_LT(dr, 1e-6) << tle.satnum << " t=" << row.t << " (km)";
            EXPECT_LT(dv, 1e-9) << tle.satnum << " t=" << row.t << " (km/s)";
            ++compared;
        }
        // If the reference stops early, the next step must fail.
        const double last = rows.back().t;
        if (c.step != 0.0 && last + 1e-6 < c.stop && rows.size() > 1) {
            const double next = std::min(last + c.step, c.stop);
            StateVector s;
            EXPECT_NE(sgp4.propagate(next, s), 0) << tle.satnum << " expected error at t=" << next;
        }
    }
    EXPECT_GT(compared, 500);
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.3e", maxPosErr);
    RecordProperty("max_position_error_km", buf);
    std::snprintf(buf, sizeof(buf), "%.3e", maxVelErr);
    RecordProperty("max_velocity_error_kms", buf);
}

TEST(Sgp4, TimeInterfaceAndDeepSpaceFlag) {
    const Tle t = Tle::parse("1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927",
                             "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537");
    const Sgp4 s(t);
    EXPECT_FALSE(s.isDeepSpace());
    const StateVector a = s.propagate(90.0);
    const StateVector b = s.propagate(t.epoch + 90.0 * 60.0);
    test::expectVecNear(a.r, b.r, 1e-6);
    EXPECT_NEAR(norm(a.r) / 1000.0, 6720.0, 30.0);
    const Sgp4 copy = s;
    test::expectVecNear(copy.propagate(10.0).r, s.propagate(10.0).r, 0.0);
}
