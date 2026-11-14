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
#include <limits>
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
                             "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537",
                             "ISS (ZARYA)");
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
    EXPECT_EQ(decodeAlpha5("J0000"), 180000); // I is skipped
    EXPECT_EQ(encodeAlpha5(100001), "A0001");
    EXPECT_EQ(encodeAlpha5(339999), "Z9999");
    EXPECT_EQ(encodeAlpha5(25544), "25544");
    EXPECT_THROW((void)Tle::parse("1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2928",
                                  "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537"),
                 Error);
    EXPECT_NO_THROW((void)Tle::parse("1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2928",
                                     "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537",
                                     "",
                                     false));
    const auto many = Tle::parseMany("ISS (ZARYA)\n"
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
    double      start{0}, stop{0}, step{0};
};
struct VerRow {
    double t;
    Vec3   r, v; // km, km/s
};
} // namespace

TEST(Sgp4, VerificationAgainstVallado) {
    // SGP4-VER.TLE and tcppver.out from Vallado et al. (2006), as distributed with python-sgp4.
    std::ifstream tf(dataPath("SGP4-VER.TLE"));
    ASSERT_TRUE(tf) << "cannot open " << dataPath("SGP4-VER.TLE");
    std::vector<VerCase> cases;
    std::string          line, pending;
    while (std::getline(tf, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }
        if (line[0] == '1') {
            pending = line;
        } else if (line[0] == '2') {
            VerCase            c{pending, line.substr(0, 69)};
            std::istringstream extra(line.size() > 69 ? line.substr(69) : "");
            extra >> c.start >> c.stop >> c.step;
            cases.push_back(c);
        }
    }
    ASSERT_GE(cases.size(), 30u);

    std::ifstream of(dataPath("tcppver.out"));
    ASSERT_TRUE(of) << "cannot open " << dataPath("tcppver.out");
    std::vector<std::pair<int, std::vector<VerRow>>> expected;
    while (std::getline(of, line)) {
        if (line.find("xx") != std::string::npos) {
            expected.push_back({std::stoi(line), {}});
            continue;
        }
        std::istringstream ls(line);
        VerRow             r;
        if (ls >> r.t >> r.r.x >> r.r.y >> r.r.z >> r.v.x >> r.v.y >> r.v.z) {
            expected.back().second.push_back(r);
        }
    }
    ASSERT_EQ(expected.size(), cases.size());

    // The cases where the reference stops before its end time, and the error code Vallado's
    // reference implementation (sgp4unit) returns at the next step: 1, the mean eccentricity has left
    // [0, 1) or the mean semi-major axis is below 0.95 Earth radii; 4, the semi-latus rectum is
    // negative ("check error code 4" in SGP4-VER.TLE); 6, the satellite has decayed (28872 "lost in
    // 50 minutes", 29141 "lost in under 420 min"). 33334, which fails at the epoch with code 3, is
    // handled below.
    const std::map<int, int> earlyStopCodes{
        {22312, 1},
        {28350, 1},
        {28872, 6},
        {29141, 6},
        {33333, 4},
        {20413, 6},
    };
    int    earlyStops = 0;
    int    compared   = 0;
    double maxPosErr = 0.0, maxVelErr = 0.0;
    for (std::size_t k = 0; k < cases.size(); ++k) {
        const VerCase& c   = cases[k];
        const Tle      tle = Tle::parse(c.line1, c.line2, "", false);
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
            const int   code = sgp4.propagate(row.t, s);
            EXPECT_EQ(code, 0) << tle.satnum << " t=" << row.t;
            const double dr = norm(s.r / 1000.0 - row.r);
            const double dv = norm(s.v / 1000.0 - row.v);
            maxPosErr       = std::max(maxPosErr, dr);
            maxVelErr       = std::max(maxVelErr, dv);
            EXPECT_LT(dr, 1e-6) << tle.satnum << " t=" << row.t << " (km)";
            EXPECT_LT(dv, 1e-9) << tle.satnum << " t=" << row.t << " (km/s)";
            ++compared;
        }
        // If the reference stops early, the next step must fail with the reference's error code.
        const double last = rows.back().t;
        if (c.step != 0.0 && last + 1e-6 < c.stop && rows.size() > 1) {
            const double next = std::min(last + c.step, c.stop);
            StateVector  s;
            const auto   expectedCode = earlyStopCodes.find(tle.satnum);
            ASSERT_NE(expectedCode, earlyStopCodes.end()) << tle.satnum << " stops early in the reference";
            EXPECT_EQ(sgp4.propagate(next, s), expectedCode->second) << tle.satnum << " at t=" << next;
            ++earlyStops;
        }
    }
    EXPECT_GT(compared, 500);
    EXPECT_EQ(earlyStops, static_cast<int>(earlyStopCodes.size()));
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.3e", maxPosErr);
    RecordProperty("max_position_error_km", buf);
    std::snprintf(buf, sizeof(buf), "%.3e", maxVelErr);
    RecordProperty("max_velocity_error_kms", buf);
}

TEST(Sgp4, TimeInterfaceAndDeepSpaceFlag) {
    const Tle  t = Tle::parse("1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927",
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

TEST(Sgp4, TimeThatIsNotFiniteIsAnError) {
    // A NaN time used to come back as a NaN state with code 0: success. An infinite one on a resonant
    // deep-space orbit would never return, because the resonance terms are integrated step by step
    // towards the time asked for. (So only NaN is tried on such an orbit before the fix.)
    const Tle iss = Tle::parse("1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927",
                               "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537");
    // Molniya 2-14, 12 h resonant, from the verification set.
    const Tle    molniya = Tle::parse("1 08195U 75081A   06176.33215444  .00000099  00000-0  11873-3 0   813",
                                   "2 08195  64.1586 279.0717 6877146 264.7651  20.2257  2.00491383225656");
    const double nan     = std::numeric_limits<double>::quiet_NaN();
    const double inf     = std::numeric_limits<double>::infinity();
    EXPECT_FALSE(Sgp4(iss).isDeepSpace());
    EXPECT_TRUE(Sgp4(molniya).isDeepSpace());
    for (const Tle& tle : {iss, molniya}) {
        const Sgp4 sgp4(tle);
        for (double minutes : {nan, inf, -inf}) {
            StateVector s{
                {1.0, 2.0, 3.0},
                {4.0, 5.0, 6.0}
            };
            EXPECT_EQ(sgp4.propagate(minutes, s), 7) << minutes;
            EXPECT_TRUE(std::isnan(s.r.x) && std::isnan(s.r.z) && std::isnan(s.v.y)) << minutes;
            try {
                (void)sgp4.propagate(minutes);
                ADD_FAILURE() << "no exception for " << minutes;
            } catch (const Sgp4Error& e) {
                EXPECT_EQ(e.code(), 7) << minutes;
            }
        }
        // An ordinary time still works after those.
        StateVector s;
        EXPECT_EQ(sgp4.propagate(10.0, s), 0);
        EXPECT_TRUE(std::isfinite(norm(s.r)));
    }
}

TEST(Sgp4, RejectsNonFiniteElements) {
    // SGP4's own range tests (e >= 1, n <= 0, ...) are all false for NaN, so a NaN element used to go through to a
    // NaN state with error code 0. Each element the theory uses must be refused when the propagator is built.
    const Tle good = Tle::parse("1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927",
                                "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537");
    struct Element {
        const char* name;
        double Tle::*member;
    };
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    for (const Element& e : {
             Element{"B*",                                    &Tle::bstar       },
             Element{"inclination",                           &Tle::inclination },
             Element{"right ascension of the ascending node", &Tle::raan        },
             Element{"eccentricity",                          &Tle::eccentricity},
             Element{"argument of perigee",                   &Tle::argPerigee  },
             Element{"mean anomaly",                          &Tle::meanAnomaly },
             Element{"mean motion",                           &Tle::meanMotion  }
    }) {
        for (const double bad : {nan, inf, -inf}) {
            Tle tle       = good;
            tle.*e.member = bad;
            try {
                const Sgp4 sgp4{tle};
                // What used to happen: the propagator was built and gave NaN with code 0.
                StateVector state;
                const int   code = sgp4.propagate(10.0, state);
                ADD_FAILURE() << e.name << " = " << bad << " was accepted; propagate gave code " << code
                              << " and r.x = " << state.r.x;
            } catch (const Error& error) {
                EXPECT_NE(std::string(error.what()).find(e.name), std::string::npos)
                    << e.name << " = " << bad << " gave: " << error.what();
            }
        }
    }
    Tle tle   = good;
    tle.epoch = Time(JulianDate{nan, 0.0}, TimeScale::UTC);
    EXPECT_THROW(Sgp4{tle}, Error) << "epoch";
    EXPECT_NO_THROW(Sgp4{good});
}

TEST(Sgp4, RejectsFiniteElementsThatGiveNoState) {
    const Tle good = Tle::parse("1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927",
                                "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537");
    // A negative mean motion is SGP4's error 2. The reference test for it (n <= 0, made on the mean motion after
    // the Kozai correction) never sees it, because the correction takes a fractional power of the negative
    // number first and gives NaN.
    Tle tle        = good;
    tle.meanMotion = -2.0;
    try {
        const Sgp4 sgp4{tle};
        ADD_FAILURE() << "a mean motion of -2 rev/day was accepted";
    } catch (const Sgp4Error& e) {
        EXPECT_EQ(e.code(), 2);
    }
    // A drag term so large that the drag coefficients overflow: there is no error code for that, and no state
    // either.
    tle       = good;
    tle.bstar = 1e300;
    EXPECT_THROW(Sgp4{tle}, Error);
}

namespace {
// ITALSAT 2 of the verification set (24 h resonant, deep space), to be given other epochs. The sidereal time and
// the Sun and Moon terms of a deep-space orbit depend on the epoch's day count, so its states show how the epoch
// field was read; a near-Earth orbit's do not.
const char* const kItalsatLine1 = "1 24208U 96044A   06177.04061740 -.00000094  00000-0  10000-3 0  1600";
const char* const kItalsatLine2 = "2 24208   3.8536  80.0121 0026640 311.0977  48.3000  1.00778054 36119";

Tle italsatWithEpoch(const char* epochField) {
    std::string line1 = kItalsatLine1;
    line1.replace(18, 14, epochField);
    return Tle::parse(line1, kItalsatLine2, "", false); // the checksum no longer fits
}
} // namespace

TEST(Sgp4, EpochOnALeapSecondDayMatchesTheReferenceCode) {
    // The reference code reads the epoch field as a day count of 86400 s days (twoline2rv: days2mdhms, then jday),
    // also on 2016-12-31, which was 86401 s long. These rows are the output of Vallado's reference SGP4.cpp
    // (version 2020-07-13, WGS-72, improved mode: the code and settings that produced tcppver.out) for the element
    // set above with its epoch field changed to 16366.50000000. km and km/s, as in tcppver.out.
    // clang-format off
    const VerRow reference[] = {
        {    0.0, {  7545.707119046, 41264.404476331,   -5.035794759},
                  {    -3.026983938287,   0.559715911174,    0.208189999904}},
        { 1440.0, {  5513.450880877, 41588.726691865,  133.462716926},
                  {    -3.050551758799,   0.410130497634,    0.208086913514}},
        {43200.0, {-40413.070053964, 11424.877684897, 2875.559644864},
                  {    -0.840794413663,  -2.957214324862,    0.020695747422}},
    };
    // clang-format on
    const Tle  tle = italsatWithEpoch("16366.50000000");
    const Sgp4 sgp4(tle);
    ASSERT_TRUE(sgp4.isDeepSpace());
    for (const VerRow& row : reference) {
        const StateVector s = sgp4.propagate(row.t);
        // The tolerances of the verification set above: 1 mm and 1 micrometre per second. A day count that is
        // half a second off (the epoch at 12:00:00.5 instead of 12:00:00) moves this orbit by 8 mm at the epoch
        // and by 33 m at 30 days.
        EXPECT_LT(norm(s.r / 1000.0 - row.r), 1e-6) << "t = " << row.t << " min (km)";
        EXPECT_LT(norm(s.v / 1000.0 - row.v), 1e-9) << "t = " << row.t << " min (km/s)";
    }
    // And the epoch is the instant the field names: 0.5 x 86400 s after 0h UTC.
    EXPECT_NEAR(tle.epoch - Time::fromCalendar(TimeScale::UTC, 2016, 12, 31, 12, 0, 0.0), 0.0, 1e-7);
}

TEST(Sgp4, EpochDayCountRunsSmoothlyAcrossALeapSecond) {
    // In the reference code's day count 2016 day 366.99999999 and 2017 day 1.00000000 are 1e-8 day apart, although
    // a leap second lies between them. Two element sets that differ only by that step in the epoch must therefore
    // give almost the same state at the same time since epoch: the reference code's own two states at 30 days are
    // 6 cm apart. With a day count that jumps by the leap second they are 68 m apart.
    const Sgp4   before(italsatWithEpoch("16366.99999999")), after(italsatWithEpoch("17001.00000000"));
    const double thirtyDays = 30.0 * 1440.0; // minutes
    EXPECT_LT(norm(before.propagate(thirtyDays).r - after.propagate(thirtyDays).r), 0.5);
    // The instants themselves are a leap second and 0.864 ms (1e-8 day) apart.
    EXPECT_NEAR(after.epoch() - before.epoch(), 1.000864, 1e-6);
}
