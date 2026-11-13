/**********************************************************************
 * astrometry_test.cpp
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

#include "astro/astrometry/apparent.h"
#include "astro/astrometry/refraction.h"
#include "astro/astrometry/rise_set.h"
#include "astro/astrometry/spherical.h"
#include "astro/core/angle.h"
#include "astro/ephemeris/low_precision.h"
#include "astro/ephemeris/spk.h"
#include "test_util.h"

using namespace astro;
using astro::test::dataPath;
using astro::test::readCsv;

TEST(Spherical, RaDecAndSeparation) {
    const RaDec rd = raDecFromVector({0.0, -2.0, 0.0});
    EXPECT_NEAR(rd.ra, 1.5 * constants::pi, 1e-15);
    EXPECT_NEAR(rd.dec, 0.0, 1e-15);
    EXPECT_NEAR(rd.distance, 2.0, 1e-15);
    test::expectVecNear(vectorFromRaDec(rd), {0.0, -2.0, 0.0}, 1e-15);
    EXPECT_NEAR(angularSeparation(0.0, 0.0, constants::pi, 0.0), constants::pi, 1e-15);
    EXPECT_NEAR(angularSeparation(1.0, 0.5, 1.0 + 1e-9, 0.5), 1e-9 * std::cos(0.5), 1e-15);
    EXPECT_NEAR(positionAngle(0.0, 0.0, 0.0, 0.1), 0.0, 1e-15);               // due north
    EXPECT_NEAR(positionAngle(0.0, 0.0, 0.1, 0.0), constants::halfPi, 1e-15); // due east
    // Ecliptic: the equinox is on both equators; the pole of the ecliptic is at dec = 90 - eps.
    const double eps  = deg2rad(23.4392911);
    const RaDec  pole = equatorialToEcliptic(unitVectorFromRaDec(1.5 * constants::pi, constants::halfPi - eps), eps);
    EXPECT_NEAR(pole.dec, constants::halfPi, 1e-12);
    test::expectVecNear(eclipticToEquatorial(equatorialToEcliptic({1, 2, 3}, eps), eps), {1, 2, 3}, 1e-14);
}

TEST(Spherical, GalacticMatchesErfa) {
    const auto csv = readCsv(dataPath("galactic.csv"));
    for (const auto& r : csv.rows) {
        const RaDec g = icrsToGalactic(std::stod(r[0]), std::stod(r[1]));
        EXPECT_NEAR(wrapPi(g.ra - std::stod(r[2])), 0.0, 1e-14);
        EXPECT_NEAR(g.dec, std::stod(r[3]), 1e-14);
        const RaDec back = galacticToIcrs(g.ra, g.dec);
        EXPECT_NEAR(wrapPi(back.ra - std::stod(r[0])), 0.0, 1e-13);
    }
    // Galactic center: l = 0, b = 0 at RA 17h45m37.2s, Dec -28d56m10s (approximately).
    const RaDec gc = galacticToIcrs(0.0, 0.0);
    EXPECT_NEAR(rad2deg(gc.ra), 266.405, 0.001);
    EXPECT_NEAR(rad2deg(gc.dec), -28.936, 0.001);
}

TEST(Refraction, ConstantsMatchErfa) {
    const auto csv = readCsv(dataPath("refco.csv"));
    for (const auto& r : csv.rows) {
        const RefractionConstants k =
            refractionConstants({std::stod(r[0]), std::stod(r[1]), std::stod(r[2]), std::stod(r[3])});
        EXPECT_NEAR(k.A, std::stod(r[4]), 1e-15);
        EXPECT_NEAR(k.B, std::stod(r[5]), 1e-15);
    }
}

TEST(Refraction, ModelsAgree) {
    const RefractionConstants k = refractionConstants({1010.0, 10.0, 0.0, 0.574});
    for (double elDeg : {15.0, 30.0, 45.0, 60.0, 85.0}) {
        const double el  = deg2rad(elDeg);
        const double obs = refractTrueToObserved(el, k);
        EXPECT_GT(obs, el);
        EXPECT_NEAR(refractObservedToTrue(obs, k), el, arcsec2rad(0.05)) << elDeg;
        // The simple formulae are rules of thumb; they agree with the rigorous model to a few percent.
        const double tol = 0.05 * (obs - el) + arcsec2rad(1.0);
        EXPECT_NEAR(obs - el, refractionSaemundsson(el), tol) << elDeg;
        EXPECT_NEAR(obs - el, refractionBennett(obs), tol) << elDeg;
    }
    // ~34 arcmin at the horizon (Bennett).
    EXPECT_NEAR(rad2arcsec(refractionBennett(0.0)) / 60.0, 34.5, 0.6);
}

TEST(Refraction, BennettAndSaemundssonMatchMeeus) {
    // Meeus (1998), Astronomical Algorithms, 2nd ed., ch. 16, with the altitude h in degrees and the
    // refraction R in arcminutes, at 1010 hPa and 10 C:
    //     R = cot(h + 7.31 / (h + 4.4)) + 0.0013515          h apparent   (16.3, Bennett)
    //     R = 1.02 cot(h + 10.3 / (h + 5.11)) + 0.0019279    h true       (16.4, Saemundsson)
    // The last term of each is the one Meeus gives to make R zero at the zenith. He prints it to
    // 1e-7 arcmin, hence the tolerance.
    const auto arcmin = [](double rad) { return rad2arcsec(rad) / 60.0; };
    for (double h = 0.0; h <= 90.0; h += 2.5) {
        const double bennett     = 1.0 / std::tan(deg2rad(h + 7.31 / (h + 4.4))) + 0.0013515;
        const double saemundsson = 1.02 / std::tan(deg2rad(h + 10.3 / (h + 5.11))) + 0.0019279;
        EXPECT_NEAR(arcmin(refractionBennett(deg2rad(h))), bennett, 1e-7) << h;
        EXPECT_NEAR(arcmin(refractionSaemundsson(deg2rad(h))), saemundsson, 1e-7) << h;
    }
    // Meeus, example 16.a: the Sun's lower limb at an apparent altitude of 30 arcmin has R = 28.754
    // arcmin; its upper limb, at a true altitude of 33.246 arcmin, has R = 24.618 arcmin. (Both are
    // printed without the zenith terms, which add 0.0014 and 0.0019 arcmin.)
    EXPECT_NEAR(arcmin(refractionBennett(deg2rad(30.0 / 60.0))), 28.754, 0.002);
    EXPECT_NEAR(arcmin(refractionSaemundsson(deg2rad(33.246 / 60.0))), 24.618, 0.0025);
    // Proportional to pressure, and to 283 / (273 + T).
    EXPECT_DOUBLE_EQ(refractionBennett(0.2, 505.0, 10.0), 0.5 * refractionBennett(0.2));
    EXPECT_DOUBLE_EQ(refractionSaemundsson(0.2, 1010.0, 293.0), 0.5 * refractionSaemundsson(0.2));
}

TEST(Refraction, BennettAndSaemundssonAreDefinedFromNadirToZenith) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    // Smallest cotangent argument of each formula, h + p / (h + q) at h = sqrt(p) - q: there the
    // refraction is largest, and below it the formulas turn back and then run away.
    struct Formula {
        const char* name;
        double (*refraction)(double, double, double) noexcept;
        double largestArcmin;
    };
    const Formula formulas[] = {
        {"Bennett",     refractionBennett,     1.0 / std::tan(deg2rad(2.0 * std::sqrt(7.31) - 4.4)) + 0.0013515  },
        {"Saemundsson", refractionSaemundsson, 1.02 / std::tan(deg2rad(2.0 * std::sqrt(10.3) - 5.11)) + 0.0019279},
    };
    for (const Formula& f : formulas) {
        SCOPED_TRACE(f.name);
        const auto refraction = [&](double elevationDeg) { return f.refraction(deg2rad(elevationDeg), 1010.0, 10.0); };
        // Every millidegree from -90 to +90 deg: finite, never negative, never rising with elevation,
        // and with no jump. The steepest slope of either formula is 0.27 rad per rad (Bennett, near
        // -0.9 deg), so no step may exceed 0.3 of the spacing.
        constexpr double stepDeg  = 0.001;
        double           previous = refraction(-90.0), lowest = previous, largestRise = 0.0, largestDrop = 0.0;
        double           riseAt = -90.0, dropAt = -90.0;
        bool             finite = std::isfinite(previous);
        for (int i = -90000 + 1; i <= 90000; ++i) {
            const double r = refraction(i * stepDeg);
            finite         = finite && std::isfinite(r);
            lowest         = std::min(lowest, r);
            if (r - previous > largestRise) {
                largestRise = r - previous;
                riseAt      = i * stepDeg;
            }
            if (previous - r > largestDrop) {
                largestDrop = previous - r;
                dropAt      = i * stepDeg;
            }
            previous = r;
        }
        EXPECT_TRUE(finite);
        EXPECT_GE(lowest, 0.0);
        EXPECT_EQ(largestRise, 0.0) << "rises at " << riseAt << " deg";
        EXPECT_LE(largestDrop, 0.3 * deg2rad(stepDeg)) << "jumps at " << dropAt << " deg";
        // Exactly zero at the zenith; the largest value everywhere below the turning point.
        EXPECT_EQ(f.refraction(constants::halfPi, 1010.0, 10.0), 0.0);
        EXPECT_NEAR(rad2arcsec(refraction(-90.0)) / 60.0, f.largestArcmin, 1e-7);
        EXPECT_EQ(refraction(-4.4), refraction(-90.0));
        EXPECT_EQ(refraction(-5.11), refraction(-90.0));
        // NaN in, NaN out.
        EXPECT_TRUE(std::isnan(f.refraction(nan, 1010.0, 10.0)));
    }
}

TEST(Apparent, ObservedPlaceMatchesErfaAtco13) {
    const SpkEphemeris eph(dataPath("de440s_2024q1_le.bsp"));
    const auto         csv = readCsv(dataPath("atco13.csv"));
    const Atmosphere   atm{1013.25, 15.0, 0.5, 0.574};
    for (const auto& r : csv.rows) {
        auto        v = [&](const char* n) { return std::stod(r[csv.col(n)]); };
        const Time  t = Time::fromJd(v("utc1"), v("utc2"), TimeScale::UTC);
        CatalogStar star;
        star.ra             = v("rc");
        star.dec            = v("dc");
        star.pmRaCosDec     = v("pr") * std::cos(star.dec);
        star.pmDec          = v("pd");
        star.parallax       = arcsec2rad(v("px"));
        star.radialVelocity = v("rv") * 1e3;
        const ConstantEop eop(EopValues{.xp = v("xp"), .yp = v("yp"), .dut1 = v("dut1")});
        const Geodetic    site{v("phi"), v("elong"), v("hm")};
        const bool        refr  = v("refraction") != 0.0;
        const Aer         a     = starObserved(star, t, site, eop, eph, refr ? &atm : nullptr);
        const double      elRef = constants::halfPi - v("zob");
        // Unrefracted: sub-mas agreement. Refracted: same model, but very low elevations are sensitive.
        const double tol = refr ? (elRef > deg2rad(10.0) ? mas2rad(1.0) : arcsec2rad(1.0)) : mas2rad(1.0);
        EXPECT_NEAR(a.el, elRef, tol) << t.iso8601() << " el=" << rad2deg(elRef);
        if (elRef < deg2rad(89.0)) {
            EXPECT_NEAR(wrapPi(a.az - v("aob")) * std::cos(elRef), 0.0, tol) << t.iso8601();
        }
    }
}

TEST(Apparent, SunAndMoonObserved) {
    const SpkEphemeris          spk(dataPath("de440s_2024q1_le.bsp"));
    const LowPrecisionEphemeris lp;
    const ConstantEop           eop;
    const Geodetic              site{deg2rad(40.0), deg2rad(-105.0), 1600.0};
    for (double day = 5.0; day < 55.0; day += 7.3) {
        const Time t  = Time::fromIso8601("2024-01-01T00:00:00Z") + day * 86400.0;
        const Aer  s1 = bodyObserved(Body::Sun, t, site, eop, spk);
        const Aer  s2 = bodyObserved(Body::Sun, t, site, eop, lp);
        EXPECT_NEAR(s1.el, s2.el, deg2rad(0.01));
        EXPECT_NEAR(s1.range / constants::au, 0.986, 0.02);
        const Aer m1 = bodyObserved(Body::Moon, t, site, eop, spk);
        const Aer m2 = bodyObserved(Body::Moon, t, site, eop, lp);
        EXPECT_NEAR(m1.el, m2.el, deg2rad(0.12));
        // Geometric vs apparent: aberration ~20" for the Sun.
        const FrameTransformer ft(t, eop);
        const Aer              g = topocentricFromGcrs(spk.position(Body::Sun, Body::Earth, t), ft, site);
        const double           d = angleBetween(aerToEnu({g.az, g.el, 1.0}), aerToEnu({s1.az, s1.el, 1.0}));
        EXPECT_NEAR(rad2arcsec(d), 20.5, 1.0);
    }
}

TEST(RiseSet, LondonSolstice) {
    // Greenwich, 2024-06-21: sunrise 03:43 UTC, sunset 20:21 UTC (published almanac values, +/- 1 min).
    const LowPrecisionEphemeris eph;
    const ConstantEop           eop;
    const Geodetic              site{deg2rad(51.4769), deg2rad(-0.0005), 0.0};
    const Time                  start = Time::fromIso8601("2024-06-21T00:00:00Z");
    const auto                  ev    = sunRiseSet(site, start, start + 86400.0, eop, eph);
    ASSERT_EQ(ev.size(), 3u);
    EXPECT_EQ(ev[0].type, SkyEventType::Rise);
    EXPECT_EQ(ev[1].type, SkyEventType::Transit);
    EXPECT_EQ(ev[2].type, SkyEventType::Set);
    EXPECT_NEAR(ev[0].time - Time::fromIso8601("2024-06-21T03:43:00Z"), 0.0, 90.0);
    EXPECT_NEAR(ev[2].time - Time::fromIso8601("2024-06-21T20:21:00Z"), 0.0, 90.0);
    EXPECT_NEAR(ev[1].time - Time::fromIso8601("2024-06-21T12:02:00Z"), 0.0, 90.0); // local noon
    // The refined events sit on the threshold.
    EXPECT_NEAR(bodyElevation(Body::Sun, ev[0].time, site, eop, eph), kSunriseElevation, 1e-6);
    // Astronomical twilight does not end at 51.5 N in midsummer.
    EXPECT_TRUE(twilight(Twilight::Astronomical, site, start, start + 86400.0, eop, eph).empty());
    EXPECT_EQ(twilight(Twilight::Civil, site, start, start + 86400.0, eop, eph).size(), 2u);
    // Polar day: no sunrise at 80 N.
    const Geodetic arctic{deg2rad(80.0), 0.0, 0.0};
    const auto     polar = sunRiseSet(arctic, start, start + 86400.0, eop, eph);
    EXPECT_TRUE(
        std::none_of(polar.begin(), polar.end(), [](const SkyEvent& e) { return e.type != SkyEventType::Transit; }));
}

TEST(RiseSet, MoonRiseSetAndPhase) {
    const SpkEphemeris eph(dataPath("de440s_2024q1_le.bsp"));
    const ConstantEop  eop;
    const Geodetic     site{deg2rad(40.0), deg2rad(-105.0), 1600.0};
    const Time         start = Time::fromIso8601("2024-02-01T00:00:00Z");
    const auto         ev    = moonRiseSet(site, start, start + 3.0 * 86400.0, eop, eph);
    int                rises = 0, sets = 0;
    for (const auto& e : ev) {
        rises += e.type == SkyEventType::Rise;
        sets  += e.type == SkyEventType::Set;
    }
    EXPECT_GE(rises, 2);
    EXPECT_GE(sets, 2);
    // Full moon 2024-01-25 17:54 UTC; new moon 2024-02-09 22:59 UTC.
    const MoonPhase full = moonPhase(Time::fromIso8601("2024-01-25T17:54:00Z"), eph);
    EXPECT_GT(full.illuminatedFraction, 0.99);
    EXPECT_NEAR(full.age, constants::pi, deg2rad(0.5));
    const MoonPhase nm = moonPhase(Time::fromIso8601("2024-02-09T22:59:00Z"), eph);
    EXPECT_LT(nm.illuminatedFraction, 0.01);
    EXPECT_NEAR(wrapPi(nm.age), 0.0, deg2rad(0.5));
    const MoonPhase wax = moonPhase(Time::fromIso8601("2024-02-14T00:00:00Z"), eph);
    EXPECT_TRUE(wax.waxing);
}

TEST(Sidereal, LocalSiderealTimeAndHourAngle) {
    const Time        t = Time::fromIso8601("2024-01-10T00:00:00Z");
    const ConstantEop eop;
    const double      lmst0 = localMeanSiderealTime(t, eop, 0.0);
    // GMST at 0h UT1 on 2024-01-10 is 7h 16m 05.626s (ERFA gmst06).
    EXPECT_NEAR(rad2deg(lmst0) / 15.0, 7.0 + 16.0 / 60.0 + 5.626 / 3600.0, 0.001 / 3600.0);
    EXPECT_NEAR(wrapPi(localMeanSiderealTime(t, eop, deg2rad(15.0)) - lmst0), deg2rad(15.0), 1e-12);
    EXPECT_NEAR(hourAngle(deg2rad(10.0), deg2rad(350.0)), deg2rad(20.0), 1e-12);
    EXPECT_LT(std::abs(localApparentSiderealTime(t, eop, 0.0) - lmst0), arcsec2rad(1.2 * 15.0));
}
