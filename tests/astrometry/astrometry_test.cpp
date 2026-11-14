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
#include <functional>
#include <limits>
#include <string>
#include <vector>

#include "astro/astrometry/apparent.h"
#include "astro/astrometry/refraction.h"
#include "astro/astrometry/rise_set.h"
#include "astro/astrometry/spherical.h"
#include "astro/core/angle.h"
#include "astro/core/error.h"
#include "astro/core/rootfind.h"
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
        // The same model as ERFA's, refracted or not and at every elevation in the file: measured
        // 0.032 mas at most (0.027 mas below 10 deg, refracted). Tolerance 0.1 mas.
        const double tol = mas2rad(0.1);
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
    // Full moon 2024-01-25 17:54 UTC; new moon 2024-02-09 22:59 UTC (the US Naval Observatory's
    // phases of the Moon, to the minute). The age (Moon minus Sun in ecliptic longitude) moves 0.0085
    // deg per minute, so the rounding to the minute is 0.0042 deg; the almanac's apparent longitudes
    // differ from these geometric ones by the Sun's aberration, 20.5" (0.0057 deg). Measured 0.0058
    // and 0.0075 deg; tolerance 0.02 deg.
    const MoonPhase full = moonPhase(Time::fromIso8601("2024-01-25T17:54:00Z"), eph);
    EXPECT_GT(full.illuminatedFraction, 0.99);
    EXPECT_NEAR(full.age, constants::pi, deg2rad(0.02));
    const MoonPhase nm = moonPhase(Time::fromIso8601("2024-02-09T22:59:00Z"), eph);
    EXPECT_LT(nm.illuminatedFraction, 0.01);
    EXPECT_NEAR(wrapPi(nm.age), 0.0, deg2rad(0.02));
    const MoonPhase wax = moonPhase(Time::fromIso8601("2024-02-14T00:00:00Z"), eph);
    EXPECT_TRUE(wax.waxing);
}

namespace {

/// A day-long cosine "elevation" of 1 rad amplitude, greatest at `peak` and once every `period`
/// seconds from it.
std::function<double(const Time&)> cosineElevation(const Time& peak, double period) {
    return [=](const Time& t) { return std::cos(constants::twoPi * (t - peak) / period); };
}

/// The transits among `events`.
std::vector<SkyEvent> transitsOf(const std::vector<SkyEvent>& events) {
    std::vector<SkyEvent> transits;
    for (const SkyEvent& e : events) {
        if (e.type == SkyEventType::Transit) {
            transits.push_back(e);
        }
    }
    return transits;
}

// A cosine of period one day is exactly 1 in double precision within 0.15 ms of its peak, so no
// search can place the peak better than that; twice that is allowed.
constexpr double kTransitTolerance = 3e-4;

} // namespace

TEST(RiseSet, TransitInTheFirstOrLastSampleIntervalIsFound) {
    const Time   peak = Time::fromIso8601("2024-03-01T00:00:00Z");
    const auto   el   = cosineElevation(peak, 86400.0);
    const double step = 600.0, span = 36000.0; // 60 sample intervals
    // The peak 0.01 s to 599 s after the start of the window, and as far before its end. Up to half
    // a step from the edge, the edge sample is higher than the one next to it.
    for (double offset : {0.01, 1.0, 100.0, 299.0, 301.0, 450.0, 599.0}) {
        SCOPED_TRACE(offset);
        const auto atStart = transitsOf(findElevationEvents(el, peak - offset, peak - offset + span, 2.0, step));
        ASSERT_EQ(atStart.size(), 1u);
        EXPECT_NEAR(atStart[0].time - peak, 0.0, kTransitTolerance);
        EXPECT_NEAR(atStart[0].elevation, 1.0, 3e-16); // 0.3 ms from the peak the cosine is 1 - 2.4e-16
        const auto atEnd = transitsOf(findElevationEvents(el, peak + offset - span, peak + offset, 2.0, step));
        ASSERT_EQ(atEnd.size(), 1u);
        EXPECT_NEAR(atEnd[0].time - peak, 0.0, kTransitTolerance);
        EXPECT_NEAR(atEnd[0].elevation, 1.0, 3e-16);
    }
}

TEST(RiseSet, WindowEdgeIsNotATransit) {
    const Time   peak = Time::fromIso8601("2024-03-01T00:00:00Z");
    const auto   el   = cosineElevation(peak, 86400.0);
    const double step = 600.0;
    // Falling all through the window, and rising all through it: the highest point of the window is
    // its edge, which is not a maximum of the elevation.
    EXPECT_TRUE(transitsOf(findElevationEvents(el, peak + 1000.0, peak + 31000.0, 2.0, step)).empty());
    EXPECT_TRUE(transitsOf(findElevationEvents(el, peak - 31000.0, peak - 1000.0, 2.0, step)).empty());
    // The maximum exactly at the start, or exactly at the end, of the window is not reported either:
    // from inside the window it looks the same.
    EXPECT_TRUE(transitsOf(findElevationEvents(el, peak, peak + 30000.0, 2.0, step)).empty());
    EXPECT_TRUE(transitsOf(findElevationEvents(el, peak - 30000.0, peak, 2.0, step)).empty());
}

TEST(RiseSet, EveryTransitIsFoundOnceByWindowsLaidEndToEnd) {
    // Day-long windows, with a peak that comes 37 s later each day: it starts 500 s before the first
    // window, crosses a window boundary on the way, and is within half a step of one for 16 days.
    const Time            t0  = Time::fromIso8601("2024-03-01T00:00:00Z");
    const double          day = 86400.0, period = day + 37.0;
    const auto            el   = cosineElevation(t0 - 500.0, period);
    const int             days = 40;
    std::vector<SkyEvent> found;
    for (int d = 0; d < days; ++d) {
        const auto transits = transitsOf(findElevationEvents(el, t0 + d * day, t0 + (d + 1) * day, 2.0, 600.0));
        found.insert(found.end(), transits.begin(), transits.end());
    }
    // Peak k is at t0 - 500 + k * period; those of k = 1 .. 39 are inside the 40 days. Peak 14 is 18 s
    // after a boundary, peak 13 is 19 s before one.
    std::vector<double> expected;
    for (int k = 0; k <= days; ++k) {
        const double s = -500.0 + k * period;
        if (s > 0.0 && s < days * day) {
            expected.push_back(s);
        }
    }
    ASSERT_EQ(expected.size(), 39u);
    ASSERT_EQ(found.size(), expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
        EXPECT_NEAR(found[i].time - t0, expected[i], kTransitTolerance) << i;
    }
}

TEST(RiseSet, SolarTransitNextToMidnightUtcIsFoundInDailyWindows) {
    // At 17 S 179 E the Sun transits within 20 minutes of 00:00 UTC all year, so day-long windows
    // from 00:00 UTC have the transit next to an edge.
    const LowPrecisionEphemeris eph;
    const ConstantEop           eop;
    const Geodetic              site{deg2rad(-17.0), deg2rad(179.0), 0.0};
    const Time                  t0   = Time::fromIso8601("2024-01-01T00:00:00Z");
    const double                day  = 86400.0;
    const int                   days = 366;
    std::vector<double>         found; // seconds after t0
    for (int d = 0; d < days; ++d) {
        for (const SkyEvent& e : transitsOf(sunRiseSet(site, t0 + d * day, t0 + (d + 1) * day, eop, eph))) {
            found.push_back(e.time - t0);
        }
    }
    // Truth, without the search under test: near each midnight the elevation has its maximum where the
    // elevations 30 s before and 30 s after are equal (exact for a parabola, which the top of the
    // Sun's path is over a minute).
    std::vector<double> expected;
    for (int d = 0; d <= days; ++d) {
        const auto slope = [&](double s) {
            return bodyElevation(Body::Sun, t0 + (s + 30.0), site, eop, eph)
                 - bodyElevation(Body::Sun, t0 + (s - 30.0), site, eop, eph);
        };
        const double s = brentRoot(slope, d * day - 1800.0, d * day + 1800.0, 1e-6);
        if (s > 0.0 && s < days * day) {
            expected.push_back(s);
        }
    }
    ASSERT_EQ(expected.size(), 366u);
    // Each expected transit is found, and nothing else is. The two times agree to 1.4 ms at worst
    // over the year; 10 ms is allowed.
    std::size_t lost = 0;
    for (double s : expected) {
        lost += std::none_of(found.begin(), found.end(), [&](double f) { return std::abs(f - s) < 0.01; });
    }
    EXPECT_EQ(lost, 0u) << "of " << expected.size() << " transits";
    EXPECT_EQ(found.size(), expected.size());
}

TEST(RiseSet, BadWindowOrStepThrows) {
    const Time t0   = Time::fromIso8601("2024-03-01T00:00:00Z");
    const auto el   = cosineElevation(t0 + 43200.0, 86400.0);
    const auto find = [&](const Time& start, const Time& end, double threshold, double step) {
        return findElevationEvents(el, start, end, threshold, step);
    };
    // The message of the astro::Error thrown, or "no error".
    const auto message = [&](const Time& start, const Time& end, double threshold, double step) -> std::string {
        try {
            (void)find(start, end, threshold, step);
        } catch (const Error& e) {
            return e.what();
        }
        return "no error";
    };
    const auto mentions = [](const std::string& text, const char* word) {
        return text.find(word) != std::string::npos;
    };
    const double nan = std::numeric_limits<double>::quiet_NaN(), inf = std::numeric_limits<double>::infinity();
    // The window as given is searched normally: two crossings of 0.5 and the transit between them.
    EXPECT_EQ(find(t0, t0 + 86400.0, 0.5, 600.0).size(), 3u);
    // End before start.
    EXPECT_TRUE(mentions(message(t0 + 86400.0, t0, 0.5, 600.0), "end")) << message(t0 + 86400.0, t0, 0.5, 600.0);
    // A step that is zero, negative, NaN or infinite.
    for (double step : {0.0, -1.0, nan, inf}) {
        EXPECT_TRUE(mentions(message(t0, t0 + 86400.0, 0.5, step), "step")) << step;
    }
    // More steps than can be counted (1e12).
    EXPECT_TRUE(mentions(message(t0, t0 + 1e9, 0.5, 1e-3), "step")) << message(t0, t0 + 1e9, 0.5, 1e-3);
    // A threshold that is NaN.
    EXPECT_TRUE(mentions(message(t0, t0 + 86400.0, nan, 600.0), "threshold"));
    // An empty window holds no event.
    EXPECT_TRUE(find(t0, t0, 0.5, 600.0).empty());
    // A step longer than the window is not an error: the window is still sampled at its ends and its
    // middle. Here the middle is the peak.
    const auto coarse = find(t0, t0 + 86400.0, 0.5, 1e9);
    ASSERT_EQ(coarse.size(), 3u);
    EXPECT_EQ(coarse[1].type, SkyEventType::Transit);
    // The wrappers pass the window on.
    const LowPrecisionEphemeris eph;
    const ConstantEop           eop;
    const Geodetic              site{deg2rad(51.4769), deg2rad(-0.0005), 0.0};
    EXPECT_THROW((void)sunRiseSet(site, t0 + 86400.0, t0, eop, eph), Error);
    EXPECT_THROW((void)twilight(Twilight::Civil, site, t0 + 86400.0, t0, eop, eph), Error);
    EXPECT_THROW((void)moonRiseSet(site, t0 + 86400.0, t0, eop, eph), Error);
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

TEST(Spherical, GalacticMatrixIsTheHipparcosOne) {
    // The Hipparcos and Tycho Catalogues (ESA 1997), vol. 1, sec. 1.5.3, eq. 1.5.11: the matrix from
    // ICRS to galactic coordinates, printed to 10 decimals.
    // clang-format off
    const Mat3 published{-0.0548755604, -0.8734370902, -0.4838350155,
                         +0.4941094279, -0.4448296300, +0.7469822445,
                         -0.8676661490, -0.1980763734, +0.4559837762};
    // clang-format on
    test::expectMatNear(icrsToGalacticMatrix(), published, 0.5e-10);
}

TEST(Spherical, EclipticMatrixIsARotationAboutTheEquinox) {
    // The ecliptic and the equator share the x axis (the equinox); the matrix turns y and z by eps
    // about it: rows (1, 0, 0), (0, cos eps, sin eps), (0, -sin eps, cos eps).
    const double eps = 84381.406 * constants::arcsecToRad; // IAU 2006 obliquity at J2000
    const double c = std::cos(eps), s = std::sin(eps);
    test::expectMatNear(equatorialToEclipticMatrix(eps), Mat3{1, 0, 0, 0, c, s, 0, -s, c}, 1e-16);
}

TEST(Apparent, LightDeflectionByTheSunIsTwiceTheNewtonianAngle) {
    // A ray from a distant source passing a mass M, seen from a distance d of it at an angle theta from
    // it, is bent away from it by (2 G M / (c^2 d)) cot(theta / 2) (to first post-Newtonian order;
    // 4 G M / (c^2 b) for impact parameter b = d theta when theta is small: 1.75" at the Sun's limb).
    const double au = constants::au;
    const double k  = 2.0 * constants::gmSun / (constants::c * constants::c * au); // rad, at 1 au
    const Vec3   e{1.0, 0.0, 0.0};                                                 // from the Sun to the observer
    for (double thetaDeg : {0.2666, 1.0, 10.0, 45.0, 90.0, 150.0}) {
        SCOPED_TRACE(testing::Message() << "theta = " << thetaDeg << " deg");
        const double theta = deg2rad(thetaDeg);
        // The source direction, theta from the direction to the Sun (-e).
        const Vec3   p{-std::cos(theta), std::sin(theta), 0.0};
        const Vec3   bent     = lightDeflection(p, p, e, 1.0);
        const double expected = k / std::tan(0.5 * theta);
        // Measured 3e-11 of the deflection (2.6e-16 rad at the limb); tolerance 1e-10 of it, plus rounding.
        EXPECT_NEAR(angleBetween(p, bent), expected, 1e-10 * expected + 1e-15);
        // Away from the Sun. The correction is added at right angles to p, not normalised (as in ERFA's
        // ld): the length is that of the hypotenuse.
        EXPECT_GT(dot(normalized(bent), e), dot(p, e));
        EXPECT_NEAR(norm(bent), std::hypot(1.0, expected), 1e-15);
    }
    // At the limb (959.63" from the centre) the deflection is 1.75".
    EXPECT_NEAR(rad2arcsec(k / std::tan(0.5 * arcsec2rad(959.63))), 1.751, 0.001);
    // Twice the mass, twice the bending; twice the distance, half.
    const Vec3 p{-std::cos(0.1), std::sin(0.1), 0.0};
    EXPECT_NEAR(angleBetween(p, lightDeflection(p, p, e, 1.0, 2.0)),
                2.0 * k / std::tan(0.05),
                1e-9 * k / std::tan(0.05));
    EXPECT_NEAR(angleBetween(p, lightDeflection(p, p, e, 2.0)), 0.5 * k / std::tan(0.05), 1e-9 * k / std::tan(0.05));
}

namespace {

/// Direction seen by an observer moving at velocity beta (units of c) of light arriving from direction
/// p (unit vector), by the Lorentz transformation:
///     p' = (p / gamma + (1 + p.beta / (1 + 1 / gamma)) beta) / (1 + p.beta).
Vec3 lorentzAberration(const Vec3& p, const Vec3& beta) {
    const double invGamma = std::sqrt(1.0 - dot(beta, beta));
    const double pb       = dot(p, beta);
    return normalized((invGamma * p + (1.0 + pb / (1.0 + invGamma)) * beta) / (1.0 + pb));
}

} // namespace

TEST(Apparent, AberrationIsTheLorentzTransformationOfDirection) {
    // ERFA's ab adds to the Lorentz transformation a term (Schwarzschild radius of the Sun / distance)
    // times beta, 2e-8 x 1e-4 = 2e-12 rad at 1 au: the tolerance is 5e-12 rad.
    const Vec3 beta = Vec3{-2.9e4, 5.1e3, 2.2e3} / constants::c; // an Earth-like orbital velocity
    for (const Vec3& dir : {
             Vec3{1,   0,    0    },
             Vec3{0,   1,    0    },
             Vec3{0.3, -0.4, 0.866},
             Vec3{-1,  0.2,  0.1  }
    }) {
        const Vec3 p   = normalized(dir);
        const Vec3 got = aberration(p, beta, 1.0);
        SCOPED_TRACE(testing::Message() << "direction " << p.x << ", " << p.y << ", " << p.z);
        EXPECT_LT(angleBetween(got, lorentzAberration(p, beta)), 5e-12);
        // Toward the velocity, by about |beta| sin(angle): up to 20.5".
        EXPECT_LE(rad2arcsec(angleBetween(got, p)), 20.6);
    }
    // At right angles to the velocity the Lorentz shift is asin(|beta|) (|beta| = 1e-4), plus ERFA's
    // 2e-12 rad term (measured 2.1e-12).
    const Vec3 perpendicular = normalized(cross(beta, Vec3{0, 0, 1}));
    EXPECT_NEAR(angleBetween(aberration(perpendicular, beta, 1.0), perpendicular), norm(beta), 1e-11);
}

TEST(Apparent, StarSpaceMotionIsLinearInMetres) {
    // A star at distance au / parallax moves in a straight line at its space velocity: the proper
    // motion times the distance across the line of sight, the radial velocity along it. After the
    // interval from the catalog epoch (here J2016.0, the epoch of Gaia DR3) to t, in Julian years.
    CatalogStar star;
    star.ra             = deg2rad(269.452);
    star.dec            = deg2rad(4.693);
    star.pmRaCosDec     = arcsec2rad(-0.802);
    star.pmDec          = arcsec2rad(10.362);
    star.parallax       = arcsec2rad(0.547);
    star.radialVelocity = -110.5e3;
    star.epochJdTdb     = 2457389.0;                                    // J2016.0
    const Time   t      = Time::fromJd(2460676.5, 0.0, TimeScale::TDB); // 2025-01-01
    const double years  = (2460676.5 - star.epochJdTdb) / 365.25;
    const double dist   = constants::au / star.parallax;
    const Vec3   u      = unitVectorFromRaDec(star.ra, star.dec);
    const Vec3   east{-std::sin(star.ra), std::cos(star.ra), 0.0};
    const Vec3   north  = cross(u, east);
    const double year   = 365.25 * 86400.0;
    const Vec3   v      = dist * (star.pmRaCosDec * east + star.pmDec * north) / year + star.radialVelocity * u;
    const Vec3   expect = dist * u + v * years * year;
    const Vec3   got    = starBarycentricPosition(star, t);
    // Rounding in coordinates of 5.6e16 m (measured: no difference at all); tolerance 1e-14 of the
    // distance.
    EXPECT_LT(norm(got - expect), 1e-14 * dist);
    EXPECT_GT(norm(got - dist * u), 1e-5 * dist); // it has moved

    // Without a parallax the result is a unit vector: the same direction with the proper motion alone.
    star.parallax       = 0.0;
    const Vec3 noPx     = starBarycentricPosition(star, t);
    const Vec3 expected = normalized(u + (star.pmRaCosDec * east + star.pmDec * north) * years);
    EXPECT_NEAR(norm(noPx), 1.0, 1e-15);
    EXPECT_LT(angleBetween(noPx, expected), 1e-14); // rounding (measured: none)
}

TEST(Apparent, ObserverEphemerisIsTheGeocentreOffsetByTheObserver) {
    // From an SPK the barycentric state is the Earth's relative to the solar-system barycenter and the
    // heliocentric position the Earth's relative to the Sun, each plus the offset. Without the
    // barycenter (LowPrecisionEphemeris) the heliocentric state stands in for it.
    const SpkEphemeris          spk(dataPath("de440s_2024q1_le.bsp"));
    const LowPrecisionEphemeris lp;
    const Time                  t = Time::fromIso8601("2024-02-10T06:00:00Z");
    const StateVector           offset{
                  {4.0e6, -3.0e6, 2.0e6 },
                  {200.0, 300.0,  -100.0}
    };
    const ObserverEphemeris o     = observerEphemeris(t, spk, offset);
    const StateVector       earth = spk.state(Body::Earth, Body::SolarSystemBarycenter, t);
    test::expectVecNear(o.barycentric.r, earth.r + offset.r, 1e-4, "barycentric r"); // 1.5e11 m
    test::expectVecNear(o.barycentric.v, earth.v + offset.v, 1e-9, "barycentric v");
    test::expectVecNear(o.heliocentric, spk.position(Body::Earth, Body::Sun, t) + offset.r, 1e-4, "heliocentric");
    const ObserverEphemeris l = observerEphemeris(t, lp);
    const StateVector       h = lp.state(Body::Earth, Body::Sun, t);
    test::expectVecNear(l.barycentric.r, h.r, 1e-4, "low precision r");
    test::expectVecNear(l.barycentric.v, h.v, 1e-9, "low precision v");
}

TEST(Apparent, StarApparentGcrsIsSpaceMotionThenDeflectionThenAberration) {
    // The pipeline of apparent.h written out with the truths of the tests above: the star's
    // barycentric position at t, seen from the geocentre (parallax), bent by the Sun (deflection
    // formula), shifted by the Earth's barycentric velocity (Lorentz transformation), as RA, Dec.
    // Barnard's star from its J2016.0 Gaia epoch, 57 deg from the Sun.
    const SpkEphemeris spk(dataPath("de440s_2024q1_le.bsp"));
    CatalogStar        star;
    star.ra             = deg2rad(269.452);
    star.dec            = deg2rad(4.693);
    star.pmRaCosDec     = arcsec2rad(-0.802);
    star.pmDec          = arcsec2rad(10.362);
    star.parallax       = arcsec2rad(0.547);
    star.radialVelocity = -110.5e3;
    star.epochJdTdb     = 2457389.0;
    const Time t        = Time::fromIso8601("2024-02-10T06:00:00Z");

    const StateVector earth = spk.state(Body::Earth, Body::SolarSystemBarycenter, t);
    const Vec3        sun   = spk.position(Body::Sun, Body::SolarSystemBarycenter, t);
    const Vec3        p     = normalized(starBarycentricPosition(star, t) - earth.r); // natural direction
    const Vec3        e     = normalized(earth.r - sun);                              // Sun to observer
    const double      theta = std::acos(dot(p, -e));
    const double      k     = 2.0 * constants::gmSun / (constants::c * constants::c * norm(earth.r - sun));
    // Bent by k cot(theta / 2) away from the Sun, in the plane of p and the Sun.
    const Vec3  away     = normalized(e - dot(e, p) * p);
    const Vec3  deflect  = normalized(p + k / std::tan(0.5 * theta) * away);
    const Vec3  apparent = lorentzAberration(deflect, earth.v / constants::c);
    const RaDec expected = raDecFromVector(apparent);
    const RaDec got      = starApparentGcrs(star, t, spk);
    // Left out here: ERFA's gravitational term in the aberration (2e-12 rad; see
    // AberrationIsTheLorentzTransformationOfDirection), and the finite distance of the star in the
    // deflection (its parallax times the deflection, 3e-6 x 9 mas). Measured 1.3e-12 rad; tolerance
    // 5e-12 rad (1 uas).
    EXPECT_NEAR(wrapPi(got.ra - expected.ra) * std::cos(expected.dec), 0.0, 5e-12);
    EXPECT_NEAR(got.dec - expected.dec, 0.0, 5e-12);
    EXPECT_NEAR(rad2deg(theta), 56.5, 0.1);
}

TEST(RiseSet, TwilightEventsSitOnTheConventionalSunElevations) {
    // Civil, nautical and astronomical twilight begin and end when the Sun's centre is 6, 12 and 18
    // deg below the horizon (geometric). At Greenwich on 2024-03-20 each happens once in the morning
    // and once in the evening, in order around sunrise and sunset.
    EXPECT_EQ(twilightElevation(Twilight::Civil), deg2rad(-6.0));
    EXPECT_EQ(twilightElevation(Twilight::Nautical), deg2rad(-12.0));
    EXPECT_EQ(twilightElevation(Twilight::Astronomical), deg2rad(-18.0));
    const LowPrecisionEphemeris eph;
    const ConstantEop           eop;
    const Geodetic              site{deg2rad(51.4769), deg2rad(-0.0005), 0.0};
    const Time                  start = Time::fromIso8601("2024-03-20T00:00:00Z");
    const auto                  sun   = sunRiseSet(site, start, start + 86400.0, eop, eph);
    ASSERT_EQ(sun.size(), 3u);
    Time dawnBefore = sun[0].time, duskAfter = sun[2].time;
    for (Twilight tw : {Twilight::Civil, Twilight::Nautical, Twilight::Astronomical}) {
        SCOPED_TRACE(testing::Message() << tw);
        const auto ev = twilight(tw, site, start, start + 86400.0, eop, eph);
        ASSERT_EQ(ev.size(), 2u);
        EXPECT_EQ(ev[0].type, SkyEventType::Rise);
        EXPECT_EQ(ev[1].type, SkyEventType::Set);
        for (const SkyEvent& e : ev) {
            // Refined to 1 ms, in which the Sun's elevation changes by 3e-7 rad at most.
            EXPECT_NEAR(bodyElevation(Body::Sun, e.time, site, eop, eph), twilightElevation(tw), 1e-6);
        }
        // Each deeper twilight begins earlier and ends later.
        EXPECT_LT(ev[0].time, dawnBefore);
        EXPECT_GT(ev[1].time, duskAfter);
        dawnBefore = ev[0].time;
        duskAfter  = ev[1].time;
    }
}

TEST(RiseSet, MoonRiseAndSetAreWhenTheRefractedUpperLimbTouchesTheHorizon) {
    // Documented: the upper limb with 34' of refraction, topocentric. At each rise and set the geometric
    // elevation of the centre plus the Moon's angular semidiameter (from its mean radius, 1737.4 km,
    // and its topocentric distance) plus 34' is zero. Transits are local maxima of the elevation.
    const SpkEphemeris eph(dataPath("de440s_2024q1_le.bsp"));
    const ConstantEop  eop;
    const Geodetic     site{deg2rad(40.0), deg2rad(-105.0), 1600.0};
    const Time         start = Time::fromIso8601("2024-02-01T00:00:00Z");
    const auto         ev    = moonRiseSet(site, start, start + 3.0 * 86400.0, eop, eph);
    int                rises = 0, sets = 0, transits = 0;
    for (const SkyEvent& e : ev) {
        SCOPED_TRACE(testing::Message() << e.type << " at " << e.time.iso8601());
        const double el = bodyElevation(Body::Moon, e.time, site, eop, eph);
        if (e.type == SkyEventType::Transit) {
            ++transits;
            EXPECT_GT(el, bodyElevation(Body::Moon, e.time - 60.0, site, eop, eph));
            EXPECT_GT(el, bodyElevation(Body::Moon, e.time + 60.0, site, eop, eph));
            continue;
        }
        rises += e.type == SkyEventType::Rise;
        sets  += e.type == SkyEventType::Set;
        const double range =
            topocentricFromGcrs(eph.position(Body::Moon, Body::Earth, e.time), FrameTransformer(e.time, eop), site)
                .range;
        const double limb = el + std::asin(1737.4e3 / range) + deg2rad(34.0 / 60.0);
        // Refined to 1 ms, in which the Moon's elevation changes by 2e-7 rad at most.
        EXPECT_NEAR(limb, 0.0, 1e-6);
    }
    // Three days hold three rises, three sets and three transits, give or take one at the ends.
    EXPECT_GE(rises, 2);
    EXPECT_GE(sets, 2);
    EXPECT_GE(transits, 2);
}

TEST(Apparent, RefractedBodyElevationIsTheUnrefractedOneRefracted) {
    // bodyObserved with an atmosphere applies the refraction of refraction.h (checked against ERFA's
    // refco in Refraction.ConstantsMatchErfa) to the unrefracted observed elevation, and leaves the
    // azimuth and the range alone.
    const SpkEphemeris eph(dataPath("de440s_2024q1_le.bsp"));
    const ConstantEop  eop;
    const Geodetic     site{deg2rad(40.0), deg2rad(-105.0), 1600.0};
    const Atmosphere   atm{830.0, 5.0, 0.3, 0.574};
    for (double hours : {16.0, 19.0, 22.5}) { // the Sun at several elevations above the horizon
        const Time t   = Time::fromIso8601("2024-02-10T00:00:00Z") + hours * 3600.0;
        const Aer  raw = bodyObserved(Body::Sun, t, site, eop, eph);
        const Aer  ref = bodyObserved(Body::Sun, t, site, eop, eph, &atm);
        SCOPED_TRACE(testing::Message() << "elevation " << rad2deg(raw.el) << " deg");
        ASSERT_GT(raw.el, deg2rad(3.0));
        EXPECT_NEAR(ref.el, refractTrueToObserved(raw.el, refractionConstants(atm)), 1e-12);
        EXPECT_GT(ref.el - raw.el, arcsec2rad(10.0)); // refraction raises the body
        EXPECT_NEAR(wrapPi(ref.az - raw.az), 0.0, 1e-12);
        EXPECT_NEAR(ref.range, raw.range, 1e-3);
    }
}

TEST(Spherical, RaDecOnTheZAxisHasRaZero) {
    // Documented: RA is 0 on the z axis; declination +-90 deg there.
    const RaDec north = raDecFromVector({0.0, 0.0, 3.0});
    EXPECT_EQ(north.ra, 0.0);
    EXPECT_EQ(north.dec, constants::halfPi);
    EXPECT_EQ(north.distance, 3.0);
    const RaDec south = raDecFromVector({0.0, 0.0, -1.0});
    EXPECT_EQ(south.ra, 0.0);
    EXPECT_EQ(south.dec, -constants::halfPi);
    // RA just short of 2 pi stays in [0, 2 pi).
    const RaDec r = raDecFromVector({1.0, -1e-17, 0.0});
    EXPECT_GE(r.ra, 0.0);
    EXPECT_LT(r.ra, constants::twoPi);
}

TEST(RiseSet, MoonIsWaningBetweenFullAndNew) {
    // Between the full moon of 2024-01-25 and the new moon of 2024-02-09 the age runs from 180 to
    // 360 deg: waning, and the illuminated fraction falls.
    const SpkEphemeris eph(dataPath("de440s_2024q1_le.bsp"));
    const MoonPhase    a = moonPhase(Time::fromIso8601("2024-01-29T00:00:00Z"), eph);
    const MoonPhase    b = moonPhase(Time::fromIso8601("2024-02-04T00:00:00Z"), eph);
    EXPECT_FALSE(a.waxing);
    EXPECT_FALSE(b.waxing);
    EXPECT_GT(a.age, constants::pi);
    EXPECT_GT(b.age, a.age);
    EXPECT_GT(a.illuminatedFraction, b.illuminatedFraction);
    // The illuminated fraction is (1 + cos(phase angle)) / 2.
    EXPECT_NEAR(a.illuminatedFraction, 0.5 * (1.0 + std::cos(a.phaseAngle)), 1e-15);
}
