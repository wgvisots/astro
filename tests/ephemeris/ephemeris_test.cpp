/**********************************************************************
 * ephemeris_test.cpp
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
#include "astro/core/error.h"
#include "astro/ephemeris/low_precision.h"
#include "astro/ephemeris/spk.h"
#include "test_util.h"

using namespace astro;
using astro::test::dataPath;
using astro::test::readCsv;

namespace {
Time tdb(double jd) { return Time::fromJd(jd, 0.0, TimeScale::TDB); }
}

class SpkTest : public ::testing::TestWithParam<const char*> {};

TEST_P(SpkTest, MatchesJplephem) {
    const SpkEphemeris spk(dataPath(GetParam()));
    EXPECT_EQ(spk.segments().size(), 5u);
    EXPECT_EQ(spk.bigEndian(), std::string(GetParam()).find("_be") != std::string::npos);
    const auto csv = readCsv(dataPath("spk_reference.csv"));
    ASSERT_EQ(csv.rows.size(), 125u);
    for (const auto& r : csv.rows) {
        const double jd = std::stod(r[0]);
        const int center = std::stoi(r[1]), target = std::stoi(r[2]);
        const Vec3 p{std::stod(r[3]) * 1e3, std::stod(r[4]) * 1e3, std::stod(r[5]) * 1e3};
        const Vec3 v{std::stod(r[6]) * 1e3, std::stod(r[7]) * 1e3, std::stod(r[8]) * 1e3};
        const StateVector s = spk.stateById(target, center, (jd - constants::jdJ2000) * 86400.0);
        EXPECT_LT(norm(s.r - p), 1e-12 * norm(p) + 1e-3) << center << "->" << target << " at " << jd;
        EXPECT_LT(norm(s.v - v), 1e-6) << center << "->" << target << " at " << jd;
    }
}

INSTANTIATE_TEST_SUITE_P(BothEndians, SpkTest,
                         ::testing::Values("de440s_2024q1_le.bsp", "de440s_2024q1_be.bsp"));

TEST(Spk, ChainingAndErrors) {
    const SpkEphemeris spk(dataPath("de440s_2024q1_le.bsp"));
    const Time t = tdb(2460330.5);
    const double et = ephemerisTime(t);
    // Moon relative to Earth = (EMB->Moon) - (EMB->Earth).
    const StateVector me = spk.state(Body::Moon, Body::Earth, t);
    const StateVector a = spk.stateById(301, 3, et), b = spk.stateById(399, 3, et);
    test::expectVecNear(me.r, a.r - b.r, 1e-6);
    EXPECT_GT(norm(me.r), 3.5e8);
    EXPECT_LT(norm(me.r), 4.1e8);
    // Sun relative to Earth goes through the SSB; reversing the pair negates the state.
    const StateVector se = spk.state(Body::Sun, Body::Earth, t);
    const StateVector es = spk.state(Body::Earth, Body::Sun, t);
    test::expectVecNear(se.r, -es.r, 1e-3);
    test::expectVecNear(se.v, -es.v, 1e-9);
    EXPECT_NEAR(norm(se.r) / constants::au, 0.984, 0.01);  // January: near perihelion
    // Type-3 segment (Mars barycenter).
    EXPECT_GT(norm(spk.state(Body::MarsBarycenter, Body::Sun, t).r), 1.3 * constants::au);
    EXPECT_THROW((void)spk.state(Body::Venus, Body::Earth, t), Error);   // not in file
    EXPECT_THROW((void)spk.state(Body::Moon, Body::Earth, tdb(2460500.5)), Error);  // out of range
    EXPECT_THROW(SpkEphemeris(dataPath("spk_reference.csv")), Error);
}

TEST(LowPrecision, SunAndMoonAgainstDe440) {
    const SpkEphemeris spk(dataPath("de440s_2024q1_le.bsp"));
    const LowPrecisionEphemeris lp;
    for (double jd = 2460311.0; jd < 2460370.0; jd += 2.3) {
        const Time t = tdb(jd);
        const Vec3 sunRef = spk.position(Body::Sun, Body::Earth, t);
        const Vec3 sun = sunPositionLowPrecision(t);
        EXPECT_LT(rad2deg(angleBetween(sun, sunRef)), 0.01) << jd;
        EXPECT_NEAR(norm(sun) / norm(sunRef), 1.0, 6e-5) << jd;
        const Vec3 moonRef = spk.position(Body::Moon, Body::Earth, t);
        const Vec3 moon = moonPositionLowPrecision(t);
        EXPECT_LT(rad2deg(angleBetween(moon, moonRef)), 0.1) << jd;
        EXPECT_NEAR(norm(moon) / norm(moonRef), 1.0, 2e-3) << jd;
        // Ephemeris interface + numerical velocity.
        const StateVector m = lp.state(Body::Moon, Body::Earth, t);
        EXPECT_LT(norm(m.v - spk.state(Body::Moon, Body::Earth, t).v), 30.0);
        test::expectVecNear(lp.position(Body::Earth, Body::Sun, t), -sun, 1e-2);
    }
    EXPECT_THROW((void)lp.state(Body::Mars, Body::Earth, tdb(2460311.0)), Error);
}
