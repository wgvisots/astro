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

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

#include "astro/core/angle.h"
#include "astro/core/error.h"
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
    const Geodetic site{degToRad(30.0), degToRad(-90.0), 0.0};
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
    const Geodetic    site{degToRad(40.0), degToRad(-105.0), 1600.0};
    const Time        start  = kIss.epoch;
    const Time        end    = start + 86400.0;
    const double      mask   = degToRad(10.0);
    const auto        passes = findPasses(sat, site, start, end, mask, 30.0);
    ASSERT_GE(passes.size(), 2u);
    for (const Pass& p : passes) {
        // A pass cut by the window starts (ends) at its edge, with the satellite above the mask there;
        // otherwise rise and set are where the elevation crosses the mask.
        EXPECT_LE(p.aos, p.tca);
        EXPECT_LE(p.tca, p.los);
        EXPECT_LT(p.los - p.aos, 15 * 60.0); // ISS passes last < 15 minutes
        if (p.startsBeforeWindow) {
            EXPECT_EQ(p.aos, start);
            EXPECT_GT(lookAngles(sat(start), site).el, mask);
        } else {
            EXPECT_NEAR(lookAngles(sat(p.aos), site).el, mask, 1e-5);
        }
        if (p.endsAfterWindow) {
            EXPECT_EQ(p.los, end);
            EXPECT_GT(lookAngles(sat(end), site).el, mask);
        } else {
            EXPECT_NEAR(lookAngles(sat(p.los), site).el, mask, 1e-5);
        }
        // Brute-force maximum elevation within the pass (and the window): 1 s scan, then 10 ms around the
        // best sample (high passes peak sharply).
        double best = -1.0, bestS = 0.0;
        for (double s = 0.0; s <= p.los - p.aos; s += 1.0) {
            const double el = lookAngles(sat(p.aos + s), site).el;
            if (el > best) {
                best  = el;
                bestS = s;
            }
        }
        for (double s = std::max(0.0, bestS - 1.0); s <= std::min(p.los - p.aos, bestS + 1.0); s += 0.01) {
            best = std::max(best, lookAngles(sat(p.aos + s), site).el);
        }
        EXPECT_NEAR(p.maxElevation, best, degToRad(0.001));
        EXPECT_GE(p.maxElevation, best - 1e-9);
    }
}

namespace {

/// A satellite flying due east at a constant height above a site's horizontal plane, passing
/// `offset` metres north of the zenith at time `peak`. `s` seconds after the peak its east, north, up
/// position is (speed s, offset, height), so every quantity of the pass has a closed form.
struct StraightPass {
    Geodetic site{degToRad(40.0), degToRad(-105.0), 1600.0};
    Time     peak   = Time::fromIso8601("2024-03-01T12:00:00Z");
    double   height = 400 * km, offset = 50 * km, speed = 7 * km;

    [[nodiscard]] ItrsStateFunction state() const {
        return [*this](const Time& t) {
            const Vec3 velocity = ecefToEnuMatrix(site.lat, site.lon).transpose() * Vec3{speed, 0.0, 0.0};
            return StateVector{enuToEcef({speed * (t - peak), offset, height}, site), velocity};
        };
    }
    /// Elevation (rad) `s` seconds after the peak; the greatest is atan(height / offset), 82.9 deg.
    [[nodiscard]] double elevation(double s) const { return std::atan2(height, std::hypot(speed * s, offset)); }
    /// Azimuth (rad, [0, 2 pi), clockwise from north) `s` seconds after the peak.
    [[nodiscard]] double azimuth(double s) const { return wrapTwoPi(std::atan2(speed * s, offset)); }
    /// Seconds from the peak to where the elevation is `mask`, on either side.
    [[nodiscard]] double crossing(double mask) const {
        return std::sqrt(std::pow(height / std::tan(mask), 2) - offset * offset) / speed;
    }
};

// Rise and set are refined to 1 ms (the tolerance given to the root finder).
constexpr double kCrossingTolerance = 1e-3;
// Angles of the closed form against those from the ITRS position: 1e-9 m of rounding in coordinates
// of 6e6 m, seen from 400 km, is 3e-15 rad.
constexpr double kAngleTolerance = 1e-13;

} // namespace

TEST(Passes, WholePassMatchesClosedForm) {
    const StraightPass p;
    const double       mask = degToRad(10.0), half = p.crossing(mask); // 324 s
    const auto         passes = findPasses(p.state(), p.site, p.peak - 900.0, p.peak + 900.0, mask, 30.0);
    ASSERT_EQ(passes.size(), 1u);
    EXPECT_FALSE(passes[0].startsBeforeWindow);
    EXPECT_FALSE(passes[0].endsAfterWindow);
    EXPECT_NEAR(passes[0].aos - p.peak, -half, kCrossingTolerance);
    EXPECT_NEAR(passes[0].los - p.peak, half, kCrossingTolerance);
    // The elevation falls off as 1.2e-3 rad/s^2 times the square of the time from the peak, so the
    // 1 ms the search promises for the time of the maximum is 1.2e-9 rad of elevation.
    EXPECT_NEAR(passes[0].tca - p.peak, 0.0, 1e-3);
    EXPECT_NEAR(passes[0].maxElevation, p.elevation(0.0), 1.2e-9);
    // The azimuth turns by 6.8e-5 rad/s at rise and set: 1 ms is 7e-8 rad.
    EXPECT_NEAR(passes[0].aosAzimuth, p.azimuth(-half), 7e-8);
    EXPECT_NEAR(passes[0].losAzimuth, p.azimuth(half), 7e-8);
}

TEST(Passes, PassAlreadyDescendingAtWindowStartPeaksAtTheStart) {
    // The window opens 100 s after the peak: within the window the elevation only falls, so its
    // maximum is at the window's start.
    const StraightPass p;
    const double       mask   = degToRad(10.0);
    const Time         start  = p.peak + 100.0;
    const auto         passes = findPasses(p.state(), p.site, start, start + 900.0, mask, 30.0);
    ASSERT_EQ(passes.size(), 1u);
    EXPECT_TRUE(passes[0].startsBeforeWindow);
    EXPECT_FALSE(passes[0].endsAfterWindow);
    EXPECT_EQ(passes[0].aos, start);
    EXPECT_EQ(passes[0].tca, start);
    EXPECT_NEAR(passes[0].maxElevation, p.elevation(100.0), kAngleTolerance); // 29.6 deg
    EXPECT_NEAR(passes[0].aosAzimuth, p.azimuth(100.0), kAngleTolerance);
    EXPECT_NEAR(passes[0].los - p.peak, p.crossing(mask), kCrossingTolerance);
}

TEST(Passes, PassStillRisingAtWindowEndPeaksAtTheEnd) {
    // The window closes 100 s before the peak.
    const StraightPass p;
    const double       mask   = degToRad(10.0);
    const Time         end    = p.peak - 100.0;
    const auto         passes = findPasses(p.state(), p.site, end - 900.0, end, mask, 30.0);
    ASSERT_EQ(passes.size(), 1u);
    EXPECT_FALSE(passes[0].startsBeforeWindow);
    EXPECT_TRUE(passes[0].endsAfterWindow);
    EXPECT_NEAR(passes[0].aos - p.peak, -p.crossing(mask), kCrossingTolerance);
    EXPECT_EQ(passes[0].los, end);
    EXPECT_EQ(passes[0].tca, end);
    EXPECT_NEAR(passes[0].maxElevation, p.elevation(-100.0), kAngleTolerance);
    EXPECT_NEAR(passes[0].losAzimuth, p.azimuth(-100.0), kAngleTolerance);
}

namespace {

/// A satellite on a circular equatorial orbit of radius r over a spherical Earth of radius R, seen from
/// a site at latitude `lat` on the prime meridian. At angle th = w (t - peak) east of the site's
/// meridian the satellite is at r (cos th, sin th, 0) in the ITRS. With the site's east (0, 1, 0),
/// north (-sin lat, 0, cos lat) and up (cos lat, 0, sin lat), the line of sight has components
///     east r sin th,   north -r cos th sin lat,   up r cos th cos lat - R,
/// its length is rho = sqrt(r^2 + R^2 - 2 r R cos th cos lat), and the range rate is
/// r R w cos lat sin th / rho. The elevation is `mask` where cos th cos lat = k with
///     k = (R cos^2 mask + sin mask sqrt(r^2 - R^2 cos^2 mask)) / r,
/// the root of the quadratic (r k - R)^2 = sin^2 mask (r^2 + R^2 - 2 r R k) on the visible side.
struct EquatorialPass {
    double   R = 6371 * km, r = R + 1000 * km, lat = degToRad(10.0), w = constants::twoPi / 6000.0;
    Time     peak = Time::fromIso8601("2024-03-01T12:00:00Z");
    Geodetic site{lat, 0.0, 0.0};

    [[nodiscard]] Ellipsoid   sphere() const { return Ellipsoid::sphere(R); }
    [[nodiscard]] double      angle(const Time& t) const { return w * (t - peak); }
    [[nodiscard]] StateVector at(double th) const {
        return {
            {r * std::cos(th),      r * std::sin(th),     0.0},
            {-r * w * std::sin(th), r * w * std::cos(th), 0.0}
        };
    }
    [[nodiscard]] ItrsStateFunction state() const {
        return [*this](const Time& t) { return at(angle(t)); };
    }
    [[nodiscard]] double east(double th) const { return r * std::sin(th); }
    [[nodiscard]] double north(double th) const { return -r * std::cos(th) * std::sin(lat); }
    [[nodiscard]] double up(double th) const { return r * std::cos(th) * std::cos(lat) - R; }
    [[nodiscard]] double range(double th) const {
        return std::sqrt(r * r + R * R - 2.0 * r * R * std::cos(th) * std::cos(lat));
    }
    [[nodiscard]] double azimuth(double th) const { return wrapTwoPi(std::atan2(east(th), north(th))); }
    [[nodiscard]] double elevation(double th) const { return std::asin(up(th) / range(th)); }
    [[nodiscard]] double rangeRate(double th) const { return r * R * w * std::cos(lat) * std::sin(th) / range(th); }
    /// cos th cos lat where the elevation is `mask`.
    [[nodiscard]] double k(double mask) const {
        const double c = std::cos(mask);
        return (R * c * c + std::sin(mask) * std::sqrt(r * r - R * R * c * c)) / r;
    }
    /// Angle th > 0 at which the satellite sets below `mask`.
    [[nodiscard]] double setAngle(double mask) const { return std::acos(k(mask) / std::cos(lat)); }
};

} // namespace

TEST(Look, CircularEquatorialOrbitOverASphereMatchesClosedForm) {
    // Look angles of EquatorialPass at angles before, at and after the meridian, and below the horizon.
    const EquatorialPass p;
    for (double th : {-0.6, -0.2, -1e-3, 0.0, 0.15, 0.4, 1.2}) {
        SCOPED_TRACE(testing::Message() << "angle " << th);
        const LookAngles la = lookAngles(p.at(th), p.site, p.sphere());
        // Coordinates of 7e6 m carry 1e-9 m of rounding: 1e-15 rad seen from 1000 km or more (measured
        // 8.9e-16 rad, 1.9e-9 m, 3.6e-12 m/s).
        EXPECT_NEAR(wrapPi(la.az - p.azimuth(th)), 0.0, 5e-15);
        EXPECT_NEAR(la.el, p.elevation(th), 5e-15);
        EXPECT_NEAR(la.range, p.range(th), 1e-8);
        EXPECT_NEAR(la.rangeRate, p.rangeRate(th), 2e-11);
    }
    // On the meridian the satellite is due south (the site is north of the orbit) and the range rate is 0.
    EXPECT_NEAR(lookAngles(p.at(0.0), p.site, p.sphere()).az, constants::pi, 1e-15);
}

TEST(Passes, CircularEquatorialOrbitOverASphereMatchesClosedForm) {
    // Rise, set, culmination and the azimuths at rise and set of EquatorialPass, from the closed form.
    const EquatorialPass p;
    const double         mask   = degToRad(10.0);
    const double         th     = p.setAngle(mask); // 0.228 rad: 218 s from the meridian
    const auto           passes = findPasses(p.state(), p.site, p.peak - 900.0, p.peak + 900.0, mask, 30.0, p.sphere());
    ASSERT_EQ(passes.size(), 1u);
    const Pass& pass = passes[0];
    EXPECT_FALSE(pass.startsBeforeWindow);
    EXPECT_FALSE(pass.endsAfterWindow);
    EXPECT_NEAR(pass.aos - p.peak, -th / p.w, kCrossingTolerance);
    EXPECT_NEAR(pass.los - p.peak, th / p.w, kCrossingTolerance);
    EXPECT_NEAR(pass.tca - p.peak, 0.0, 1e-3);
    // Highest on the meridian: tan(el) = (r cos lat - R) / (r sin lat), 34.7 deg. The elevation is
    // flat there (its second derivative is about w^2 = 1e-6 rad/s^2), so even 1 ms of timing is 5e-13
    // rad; measured 1.1e-16 rad (the time came within 1.3e-6 s).
    EXPECT_NEAR(pass.maxElevation, std::atan2(p.r * std::cos(p.lat) - p.R, p.r * std::sin(p.lat)), 1e-12);
    // The azimuth turns by 4e-3 rad/s at rise and set: 1 ms is 4e-6 rad.
    EXPECT_NEAR(pass.aosAzimuth, p.azimuth(-th), 4e-6); // south-west
    EXPECT_NEAR(pass.losAzimuth, p.azimuth(th), 4e-6);  // south-east
    // The same circle seen as coverage: the central angle from the site to the sub-satellite point at
    // rise is acos(k), and the slant range then is rho at that angle.
    const double central = std::acos(p.k(mask));
    EXPECT_NEAR(coverageHalfAngle(p.r - p.R, mask, p.R), central, 2e-15); // measured 3.9e-16
    EXPECT_NEAR(maxSlantRange(p.r - p.R, mask, p.R), p.range(th), 1e-8);  // measured 2.8e-9 m
}

TEST(Look, SatelliteAtTheSiteHasNoDirectionAndNaNRangeRate) {
    // Documented: at zero range the range rate does not exist and is NaN.
    const Geodetic   site{degToRad(30.0), degToRad(-90.0), 0.0};
    const LookAngles la = lookAngles(
        {
            geodeticToEcef(site),
            {1.0, 2.0, 3.0}
    },
        site);
    EXPECT_EQ(la.range, 0.0);
    EXPECT_TRUE(std::isnan(la.rangeRate));
}

TEST(Ground, GroundTrackOfACircularEquatorialOrbitRunsAlongTheEquator) {
    // The sub-satellite point of EquatorialPass is on the equator at longitude th, at height r - R, at
    // start, start + step, ... and at `end`.
    const EquatorialPass p;
    const auto           sat   = [&](const Time& t) { return p.state()(t).r; };
    const Time           start = p.peak - 1000.0;
    const auto           track = groundTrack(sat, start, p.peak + 1100.0, 600.0, p.sphere());
    ASSERT_EQ(track.size(), 5u); // -1000, -400, 200, 800 s and the end at 1100 s
    const double times[] = {-1000.0, -400.0, 200.0, 800.0, 1100.0};
    for (std::size_t i = 0; i < track.size(); ++i) {
        SCOPED_TRACE(testing::Message() << "point " << i);
        EXPECT_NEAR(track[i].lat, 0.0, 1e-15);
        // The times differ from `peak` by 1e-11 s of Julian-date rounding, 1e-14 rad at w = 1e-3 rad/s
        // (measured 1.4e-14 rad).
        EXPECT_NEAR(wrapPi(track[i].lon - p.w * times[i]), 0.0, 5e-14);
        EXPECT_NEAR(track[i].h, p.r - p.R, 1e-8);
    }
}

TEST(Passes, MaximumInsideAWindowThatCutsThePassIsFound) {
    // The window opens 7 s before the peak, or closes 7 s after it: the pass is cut, and its maximum
    // is inside the window, in the first or the last sample interval.
    const StraightPass p;
    const double       mask = degToRad(10.0);
    const auto         late = findPasses(p.state(), p.site, p.peak - 7.0, p.peak + 893.0, mask, 30.0);
    ASSERT_EQ(late.size(), 1u);
    EXPECT_TRUE(late[0].startsBeforeWindow);
    EXPECT_EQ(late[0].aos, p.peak - 7.0);
    EXPECT_NEAR(late[0].tca - p.peak, 0.0, 1e-3);
    EXPECT_NEAR(late[0].maxElevation, p.elevation(0.0), 1.2e-9); // see WholePassMatchesClosedForm
    const auto early = findPasses(p.state(), p.site, p.peak - 893.0, p.peak + 7.0, mask, 30.0);
    ASSERT_EQ(early.size(), 1u);
    EXPECT_TRUE(early[0].endsAfterWindow);
    EXPECT_EQ(early[0].los, p.peak + 7.0);
    EXPECT_NEAR(early[0].tca - p.peak, 0.0, 1e-3);
    EXPECT_NEAR(early[0].maxElevation, p.elevation(0.0), 1.2e-9);
}

namespace {

/// Apparent radii of the Sun (a) and the Earth (b), and the angle between their centres (c), seen from
/// `sat`, written out from the geometry (Montenbruck & Gill 2000, sec. 3.4.2).
struct Disks {
    double a, b, c;
};

Disks disksSeenFrom(const Vec3& sat, const Vec3& sun) {
    const Vec3 toSun = sun - sat;
    return {std::asin(constants::sunRadius / norm(toSun)),
            std::asin(constants::wgs84A / norm(sat)),
            std::acos(dot(-sat, toSun) / (norm(sat) * norm(toSun)))};
}

/// Fraction of a disk of radius `a` left uncovered by a disk of radius `b` whose centre is `c` away,
/// by integrating across the line of centres the overlap of the two disks' vertical chords (midpoint
/// rule, `n` strips). A route independent of the closed-form lens area the library uses.
double uncoveredFraction(double a, double b, double c, int n = 1000000) {
    const double lo = std::max(-a, c - b), hi = std::min(a, c + b);
    if (hi <= lo) {
        return 1.0;
    }
    const double h       = (hi - lo) / n;
    double       overlap = 0.0;
    for (int i = 0; i < n; ++i) {
        const double x   = lo + (i + 0.5) * h;
        const double ha  = std::sqrt(std::max(0.0, a * a - x * x));
        const double hb  = std::sqrt(std::max(0.0, b * b - (x - c) * (x - c)));
        overlap         += 2.0 * std::min(ha, hb) * h;
    }
    return 1.0 - overlap / (constants::pi * a * a);
}

/// The root of `f` in [lo, hi] (f changes sign there), by bisection to the resolution of a double.
template <class F>
double bisect(F&& f, double lo, double hi) {
    const bool loNegative = f(lo) < 0.0;
    for (int i = 0; i < 200; ++i) {
        const double mid                         = 0.5 * (lo + hi);
        ((f(mid) < 0.0) == loNegative ? lo : hi) = mid;
    }
    return 0.5 * (lo + hi);
}

} // namespace

TEST(Eclipse, ShadowStatesAndBetaAngle) {
    const Vec3   sun{constants::au, 0, 0};
    const double r = 7000 * km;
    // Sunward: lit. Directly behind: umbra.
    EXPECT_EQ(earthShadow({r, 0, 0}, sun).state, ShadowState::Sunlit);
    EXPECT_EQ(earthShadow({r, 0, 0}, sun).illumination, 1.0);
    EXPECT_EQ(earthShadow({-r, 0, 0}, sun).state, ShadowState::Umbra);
    EXPECT_EQ(earthShadow({-r, 0, 0}, sun).illumination, 0.0);
    EXPECT_EQ(earthShadow({-r, 0, 0}, sun, ShadowModel::Cylindrical).state, ShadowState::Umbra);
    // The cylinder's edge is the Earth's radius from the Sun line, whatever the distance behind.
    EXPECT_EQ(earthShadow({-r, constants::wgs84A - 1.0, 0}, sun, ShadowModel::Cylindrical).state, ShadowState::Umbra);
    EXPECT_EQ(earthShadow({-r, constants::wgs84A + 1.0, 0}, sun, ShadowModel::Cylindrical).state, ShadowState::Sunlit);
    // Far beyond the umbra cone (~1.4e6 km), the Earth only partially covers the Sun: annular, with the
    // Earth's whole disk inside the Sun's, so the visible fraction is 1 - b^2 / a^2.
    const Vec3   far{-3.0e9, 0, 0};
    const Disks  d = disksSeenFrom(far, sun);
    ShadowResult s = earthShadow(far, sun);
    EXPECT_EQ(s.state, ShadowState::Annular);
    EXPECT_NEAR(s.illumination, 1.0 - d.b * d.b / (d.a * d.a), 1e-12);
    // Beta angle.
    const StateVector polar{
        {r, 0, 0       },
        {0, 0, 7.5 * km}
    }; // orbit normal along -y
    EXPECT_NEAR(betaAngle(polar, {0, -constants::au, 0}), constants::halfPi, 1e-12);
    EXPECT_NEAR(betaAngle(polar, sun), 0.0, 1e-12);
    EXPECT_NEAR(betaAngle(polar, {constants::au, -constants::au, 0}), constants::pi / 4.0, 1e-12);
}

TEST(Eclipse, PenumbraBoundariesAndIlluminationMatchTheDiskGeometry) {
    // Across the shadow 7000 km behind the Earth, along y: the penumbra begins where the disks touch
    // (c = a + b) and the umbra where the Sun's disk is just inside the Earth's (c = b - a). Both edges
    // are found here by bisection on the geometry; inside the penumbra the illumination is the
    // uncovered fraction of the solar disk, integrated numerically.
    const Vec3   sun{constants::au, 0, 0};
    const double r     = 7000 * km;
    const auto   point = [&](double y) { return Vec3{-r, y, 0}; };
    const auto   outer = [&](double y) {
        const Disks d = disksSeenFrom(point(y), sun);
        return d.c - (d.a + d.b);
    };
    const auto inner = [&](double y) {
        const Disks d = disksSeenFrom(point(y), sun);
        return d.c - (d.b - d.a);
    };
    const double yPenumbra = bisect(outer, 6000 * km, 7000 * km);
    const double yUmbra    = bisect(inner, 6000 * km, 7000 * km);
    ASSERT_LT(yUmbra, yPenumbra); // the penumbra is about 61 km wide here
    // 1 m on either side of each edge (the edges move 1 m per 1.4e-7 rad of c, far above rounding).
    EXPECT_EQ(earthShadow(point(yPenumbra + 1.0), sun).state, ShadowState::Sunlit);
    EXPECT_EQ(earthShadow(point(yPenumbra - 1.0), sun).state, ShadowState::Penumbra);
    EXPECT_EQ(earthShadow(point(yUmbra + 1.0), sun).state, ShadowState::Penumbra);
    EXPECT_EQ(earthShadow(point(yUmbra - 1.0), sun).state, ShadowState::Umbra);
    // Illumination at a quarter, half and three quarters of the way across, increasing outward. The
    // chords have square-root ends, so the midpoint rule converges as n^-1.5: measured differences
    // 1.5e-8 with 2e5 strips, 2e-10 with 4e6. With 1e6 strips (a few ms), tolerance 5e-9.
    double previous = 0.0;
    for (double f : {0.25, 0.5, 0.75}) {
        const double y = yUmbra + f * (yPenumbra - yUmbra);
        const Disks  d = disksSeenFrom(point(y), sun);
        SCOPED_TRACE(testing::Message() << "fraction " << f << " of the way out");
        const ShadowResult s = earthShadow(point(y), sun);
        EXPECT_EQ(s.state, ShadowState::Penumbra);
        EXPECT_NEAR(s.illumination, uncoveredFraction(d.a, d.b, d.c), 5e-9);
        EXPECT_GT(s.illumination, previous);
        previous = s.illumination;
    }
}

TEST(Eclipse, ConicalAndCylindricalEclipseTimesMatchTheGeometry) {
    // A circular equatorial orbit 500 km up with the Sun fixed along x: the satellite is at angle
    // u = n t from x. Cylindrical shadow: entry where r sin u = R (u past 90 deg). Conical: penumbra
    // entry where c = a + b, umbra entry where c = b - a, both solved for u by bisection here. Exits are
    // the mirror images about u = 180 deg.
    const double      r      = constants::wgs84A + 500 * km;
    const StateVector s0     = elementsToState({r, 0.0, 0.0, 0.0, 0.0, 0.0});
    const double      period = orbitalPeriod(r);
    const double      n      = constants::twoPi / period;
    const Time        t0;
    auto              sat = [&](const Time& t) { return propagateKepler(s0, t - t0).r; };
    const Vec3        sunAt{constants::au, 0, 0};
    auto              sun       = [&](const Time&) { return sunAt; };
    const auto        position  = [&](double u) { return Vec3{r * std::cos(u), r * std::sin(u), 0.0}; };
    const double      uCylinder = constants::pi - std::asin(constants::wgs84A / r);
    const double      uPenumbra = bisect(
        [&](double u) {
            const Disks d = disksSeenFrom(position(u), sunAt);
            return d.c - (d.a + d.b);
        },
        constants::halfPi,
        constants::pi);
    const double uUmbra = bisect(
        [&](double u) {
            const Disks d = disksSeenFrom(position(u), sunAt);
            return d.c - (d.b - d.a);
        },
        constants::halfPi,
        constants::pi);
    // Brent is asked for 1 ms and, converging faster than linearly, gives each crossing to 4.3e-5 s
    // (measured; 2.3e-6 s for the conical ones). Tolerance 2e-4 s.
    constexpr double tol = 2e-4;

    const auto cyl = findEclipses(sat, sun, t0, t0 + 2.0 * period, 20.0, ShadowModel::Cylindrical);
    ASSERT_EQ(cyl.size(), 2u);
    for (const EclipseInterval& e : cyl) {
        EXPECT_EQ(e.umbraEntry, e.penumbraEntry); // one boundary only, as documented
        EXPECT_EQ(e.umbraExit, e.penumbraExit);
    }
    EXPECT_NEAR(cyl[0].penumbraEntry - t0, uCylinder / n, tol);
    EXPECT_NEAR(cyl[0].penumbraExit - t0, (constants::twoPi - uCylinder) / n, tol);
    EXPECT_NEAR(cyl[1].penumbraEntry - t0, period + uCylinder / n, tol);

    const auto con = findEclipses(sat, sun, t0, t0 + 2.0 * period, 20.0, ShadowModel::Conical);
    ASSERT_EQ(con.size(), 2u);
    EXPECT_TRUE(con[0].umbral);
    EXPECT_NEAR(con[0].penumbraEntry - t0, uPenumbra / n, tol);
    EXPECT_NEAR(con[0].umbraEntry - t0, uUmbra / n, tol); // 8.4 s after the penumbra entry
    EXPECT_NEAR(con[0].umbraExit - t0, (constants::twoPi - uUmbra) / n, tol);
    EXPECT_NEAR(con[0].penumbraExit - t0, (constants::twoPi - uPenumbra) / n, tol);
}

TEST(Eclipse, EclipseInProgressAtTheWindowEdgesIsClippedToThem) {
    // The window opens and closes in mid-eclipse (cylindrical shadow, so entry and exit are known in
    // closed form): the first interval starts at the window's start and the last ends at its end; the
    // crossing inside the window is the analytic one.
    const double      r      = constants::wgs84A + 500 * km;
    const StateVector s0     = elementsToState({r, 0.0, 0.0, 0.0, 0.0, 0.0});
    const double      period = orbitalPeriod(r);
    const double      n      = constants::twoPi / period;
    const Time        t0;
    auto              sat = [&](const Time& t) { return propagateKepler(s0, t - t0).r; };
    auto              sun = [](const Time&) { return Vec3{constants::au, 0, 0}; };
    const double      uIn = constants::pi - std::asin(constants::wgs84A / r), uOut = constants::twoPi - uIn;
    const Time        start = t0 + constants::pi / n;          // behind the Earth
    const Time        end   = t0 + period + constants::pi / n; // behind the Earth again, one orbit later
    for (ShadowModel model : {ShadowModel::Cylindrical, ShadowModel::Conical}) {
        SCOPED_TRACE(testing::Message() << model);
        const auto ecl = findEclipses(sat, sun, start, end, 20.0, model);
        ASSERT_EQ(ecl.size(), 2u);
        EXPECT_EQ(ecl[0].penumbraEntry, start);
        EXPECT_EQ(ecl[0].umbraEntry, start);
        EXPECT_TRUE(ecl[0].umbral);
        EXPECT_EQ(ecl[1].penumbraExit, end);
        EXPECT_EQ(ecl[1].umbraExit, end);
        if (model == ShadowModel::Cylindrical) {
            EXPECT_NEAR(ecl[0].penumbraExit - t0,
                        uOut / n,
                        2e-4); // as in ConicalAndCylindricalEclipseTimesMatchTheGeometry
            EXPECT_NEAR(ecl[1].penumbraEntry - t0,
                        period + uIn / n,
                        2e-4); // as in ConicalAndCylindricalEclipseTimesMatchTheGeometry
        }
    }
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
    const Time end   = start + 86400.0;
    const auto state = [&](const Time& t) { return earthShadow(sat(t), sun(t)).state; };
    for (const auto& e : ecl) {
        // An eclipse cut by the window starts (or ends) at the window's edge, where the satellite is
        // then in shadow, and is no longer than a whole one. One that is not cut has a sunlit moment
        // just outside each end, and lasts 20 to 40 minutes.
        const bool   cutAtStart = e.penumbraEntry == start, cutAtEnd = e.penumbraExit == end;
        const double duration = e.penumbraExit - e.penumbraEntry;
        EXPECT_LT(duration, 40 * 60.0);
        if (cutAtStart) {
            EXPECT_NE(state(start), ShadowState::Sunlit);
        } else {
            EXPECT_EQ(state(e.penumbraEntry - 1.0), ShadowState::Sunlit);
            EXPECT_NE(state(e.penumbraEntry + 1.0), ShadowState::Sunlit);
        }
        if (cutAtEnd) {
            EXPECT_NE(state(end), ShadowState::Sunlit);
        } else {
            EXPECT_EQ(state(e.penumbraExit + 1.0), ShadowState::Sunlit);
            EXPECT_NE(state(e.penumbraExit - 1.0), ShadowState::Sunlit);
        }
        if (!cutAtStart && !cutAtEnd) {
            EXPECT_GT(duration, 20 * 60.0);
        }
        EXPECT_EQ(state(e.umbraEntry + 0.5 * (e.umbraExit - e.umbraEntry)), ShadowState::Umbra);
    }
}

TEST(Ground, SubpointFootprintAndLineOfSight) {
    // GEO over lon 0.
    const Geodetic g = subSatellitePoint({42164.17 * km, 0, 0});
    EXPECT_NEAR(g.lat, 0.0, 1e-12);
    EXPECT_NEAR(g.lon, 0.0, 1e-12);
    EXPECT_NEAR(g.h, 42164.17 * km - constants::wgs84A, 1e-6);
    // GEO coverage to the horizon: acos(R / r) ~ 81.3 deg.
    EXPECT_NEAR(radToDeg(coverageHalfAngle(42164.17 * km - constants::wgs84A)), 81.30, 0.01);
    EXPECT_LT(coverageHalfAngle(500 * km, degToRad(10.0)), coverageHalfAngle(500 * km));
    EXPECT_NEAR(maxSlantRange(500 * km) / km,
                std::sqrt(std::pow(constants::wgs84A + 500 * km, 2) - std::pow(constants::wgs84A, 2)) / km,
                1e-6);
    // Footprint points are all at the requested central angle.
    const double lam = degToRad(20.0);
    for (const Geodetic& p : coverageFootprint(degToRad(45.0), degToRad(10.0), lam, 36)) {
        const double c = std::acos(std::sin(degToRad(45.0)) * std::sin(p.lat)
                                   + std::cos(degToRad(45.0)) * std::cos(p.lat) * std::cos(p.lon - degToRad(10.0)));
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

TEST(Ground, FootprintAboutAPoleIsAnEvenlySpacedParallel) {
    // About a pole, the circle of angular radius `lam` is the parallel of latitude 90 deg - lam, and
    // its n points are 360/n deg of longitude apart. Every azimuth at a pole points along a meridian:
    // the one reached is the limit for a center just short of the pole on the meridian `lon0`. So from
    // the north pole azimuth 180 deg (south) runs down `lon0` itself and azimuth 0 down the opposite
    // meridian, and from the south pole azimuth 0 (north) runs up `lon0`.
    const double lam = degToRad(20.0), lon0 = degToRad(10.0);
    const int    n = 12;
    for (double pole : {constants::halfPi, -constants::halfPi}) {
        const auto ring = coverageFootprint(pole, lon0, lam, n);
        // The same circle about a center 1e-9 rad (6 mm) short of the pole.
        const auto nearby = coverageFootprint(pole > 0.0 ? pole - 1e-9 : pole + 1e-9, lon0, lam, n);
        ASSERT_EQ(ring.size(), 12u);
        ASSERT_EQ(nearby.size(), 12u);
        for (std::size_t i = 0; i < ring.size(); ++i) {
            SCOPED_TRACE(testing::Message() << (pole > 0.0 ? "north" : "south") << " pole, point " << i);
            const double az  = constants::twoPi * static_cast<double>(i) / n;
            const double lon = pole > 0.0 ? lon0 + constants::pi - az : lon0 + az;
            // 1e-14 rad: asin near 1 magnifies rounding by 1 / sin(lam) = 3.
            EXPECT_NEAR(ring[i].lat, pole > 0.0 ? constants::halfPi - lam : lam - constants::halfPi, 1e-14);
            EXPECT_NEAR(wrapPi(ring[i].lon - lon), 0.0, 1e-14);
            // Moving the center by 1e-9 rad moves no point of the circle by more than 1e-9 rad of
            // latitude, or 1e-9 / sin(lam) = 3e-9 rad of longitude.
            EXPECT_NEAR(nearby[i].lat, ring[i].lat, 1e-8);
            EXPECT_NEAR(wrapPi(nearby[i].lon - ring[i].lon), 0.0, 1e-8);
        }
    }
}

namespace {

/// The message of the astro::Error that `call` throws, or "no error" if it throws nothing.
template <class Call>
std::string errorMessage(Call&& call) {
    try {
        call();
    } catch (const Error& e) {
        return e.what();
    }
    return "no error";
}

/// Whether `call` throws an astro::Error whose message names both `function` and `what`.
template <class Call>
testing::AssertionResult throwsNaming(Call&& call, const std::string& function, const std::string& what) {
    const std::string message = errorMessage(call);
    if (message.find(function) != std::string::npos && message.find(what) != std::string::npos) {
        return testing::AssertionSuccess();
    }
    return testing::AssertionFailure() << "got \"" << message << "\"";
}

constexpr double kNan = std::numeric_limits<double>::quiet_NaN();
constexpr double kInf = std::numeric_limits<double>::infinity();

} // namespace

TEST(Passes, BadWindowOrStepThrows) {
    const StraightPass p;
    const auto         sat = p.state();
    const Time         a = p.peak - 600.0, b = p.peak + 600.0;
    // The satellite is up at `b`: the old code returned a "pass" from b back to a.
    EXPECT_TRUE(throwsNaming([&] { (void)findPasses(sat, p.site, b, a); }, "findPasses", "end"));
    for (double step : {0.0, -60.0, kNan, kInf}) {
        EXPECT_TRUE(throwsNaming([&] { (void)findPasses(sat, p.site, a, b, 0.0, step); }, "findPasses", "step"))
            << step;
    }
    EXPECT_TRUE(throwsNaming([&] { (void)findPasses(sat, p.site, a, b, kNan); }, "findPasses", "minElevation"));
    // More steps than can be counted.
    EXPECT_THROW((void)findPasses(sat, p.site, a, a + 1e9, 0.0, 1e-3), Error);
    // An empty window holds no pass, even at a moment when the satellite is up.
    EXPECT_TRUE(findPasses(sat, p.site, p.peak, p.peak).empty());
}

TEST(Eclipse, BadWindowOrStepThrows) {
    const double      r  = constants::wgs84A + 500 * km;
    const StateVector s0 = elementsToState({r, 0.0, 0.0, 0.0, 0.0, 0.0});
    const Time        t0;
    auto              sat = [&](const Time& t) { return propagateKepler(s0, t - t0).r; };
    auto              sun = [](const Time&) { return Vec3{constants::au, 0, 0}; };
    const Time        mid = t0 + 0.5 * orbitalPeriod(r); // behind the Earth: in eclipse
    ASSERT_EQ(earthShadow(sat(mid), sun(mid)).state, ShadowState::Umbra);
    EXPECT_TRUE(throwsNaming([&] { (void)findEclipses(sat, sun, mid, t0); }, "findEclipses", "end"));
    for (double step : {0.0, -30.0, kNan, kInf}) {
        EXPECT_TRUE(throwsNaming([&] { (void)findEclipses(sat, sun, t0, mid, step); }, "findEclipses", "step")) << step;
    }
    EXPECT_TRUE(throwsNaming([&] { (void)findEclipses(sat, sun, t0, t0 + 1e9, 1e-3); }, "findEclipses", "step"));
    // An empty window holds no eclipse, even at a moment in the umbra.
    EXPECT_TRUE(findEclipses(sat, sun, mid, mid).empty());
}

TEST(Ground, GroundTrackBadWindowOrStepThrows) {
    auto       geo = [](const Time&) { return Vec3{42164.17 * km, 0, 0}; };
    const Time t0;
    EXPECT_TRUE(throwsNaming([&] { (void)groundTrack(geo, t0 + 3600.0, t0, 600.0); }, "groundTrack", "end"));
    for (double step : {0.0, -600.0, kNan, kInf}) {
        EXPECT_TRUE(throwsNaming([&] { (void)groundTrack(geo, t0, t0 + 3600.0, step); }, "groundTrack", "step"))
            << step;
    }
    EXPECT_TRUE(throwsNaming([&] { (void)groundTrack(geo, t0, t0 + 1e9, 1e-3); }, "groundTrack", "step"));
    // start, start + step, ..., and `end` itself when it is not a whole number of steps from start.
    EXPECT_EQ(groundTrack(geo, t0, t0 + 3600.0, 600.0).size(), 7u);
    EXPECT_EQ(groundTrack(geo, t0, t0 + 3700.0, 600.0).size(), 8u);
    // An empty window is the one point at `start`.
    EXPECT_EQ(groundTrack(geo, t0, t0, 600.0).size(), 1u);
}

namespace {

/// A state whose orbital frames are written out by hand: at (7000 km, 0, 0), moving 1 km/s outward
/// and 7 km/s along the horizontal direction (0, cos i, sin i), inclination i = 30 deg. Its angular
/// momentum is along (0, -sin i, cos i).
struct HandState {
    double      c = std::sqrt(3.0) / 2.0, s = 0.5, k = 1.0 / std::sqrt(50.0); // cos i, sin i, 1 / |v| in km/s
    StateVector state{
        {7000 * km, 0.0,       0.0      },
        {1 * km,    7 * km* c, 7 * km* s}
    };
    Vec3 radial{1.0, 0.0, 0.0};
    Vec3 horizontal{0.0, c, s};           // in the orbit plane, perpendicular to r, prograde
    Vec3 normal{0.0, -s, c};              // along r x v
    Vec3 velocity{k, 7 * k* c, 7 * k* s}; // v / |v|
    Vec3 coNormal{7 * k, -k* c, -k* s};   // velocity x normal
};

} // namespace

TEST(OrbitalFrames, MatricesOfAnInclinedStateAreKnownByHand) {
    // Each row of each matrix written out from the definitions in orbital_frames.h for HandState.
    const HandState h;
    // 1e-15: rounding in unit vectors.
    expectMatNear(ricMatrix(h.state), Mat3::fromRows(h.radial, h.horizontal, h.normal), 1e-15, "RIC");
    expectMatNear(lvlhMatrix(h.state), Mat3::fromRows(h.horizontal, -h.normal, -h.radial), 1e-15, "LVLH");
    expectMatNear(vncMatrix(h.state), Mat3::fromRows(h.velocity, h.normal, h.coNormal), 1e-15, "VNC");
    expectMatNear(ntwMatrix(h.state), Mat3::fromRows(h.coNormal, h.velocity, h.normal), 1e-15, "NTW");
}

TEST(OrbitalFrames, MatricesAreRightHandedRotations) {
    const StateVector s = elementsToState({7000 * km, 0.1, 0.8, 0.3, 0.5, 1.0});
    for (const Mat3& m : {ricMatrix(s), lvlhMatrix(s), vncMatrix(s), ntwMatrix(s)}) {
        expectOrthonormal(m, 1e-14);
    }
}

TEST(OrbitalFrames, RelativeStateOfADeputyOnTheSameCircularOrbit) {
    // Chief and deputy on one circular orbit (radius r, inclination 30 deg), the deputy du ahead. In the
    // chief's rotating RIC frame the deputy sits still at (r (cos du - 1), r sin du, 0): its velocity
    // relative to the rotating frame is zero.
    const double      r = 7000 * km, v = std::sqrt(constants::gmEarth / r), du = 1e-4;
    const double      c = std::sqrt(3.0) / 2.0, s = 0.5;
    const StateVector chief{
        {r,   0.0,   0.0  },
        {0.0, v * c, v * s}
    };
    const StateVector deputy{
        {r * std::cos(du),  r * std::sin(du) * c, r * std::sin(du) * s},
        {-v * std::sin(du), v * std::cos(du) * c, v * std::cos(du) * s}
    };
    const StateVector rel = relativeStateRic(chief, deputy);
    // Positions of 7e6 m carry 1e-9 m of rounding (measured 3.6e-10 m); velocities of 7.5e3 m/s carry
    // 1e-12 m/s (measured 9e-13 m/s). Tolerances 2e-9 m and 3e-12 m/s.
    expectVecNear(rel.r, {r * (std::cos(du) - 1.0), r * std::sin(du), 0.0}, 2e-9, "position"); // x is -0.035 m
    expectVecNear(rel.v, {0.0, 0.0, 0.0}, 3e-12, "velocity");
    expectVecNear(relativeRicToDeputy(chief, rel).r, deputy.r, 2e-9, "back r");
    expectVecNear(relativeRicToDeputy(chief, rel).v, deputy.v, 3e-12, "back v");

    // A deputy given an extra 0.1 m/s along the orbit normal moves at 0.1 m/s along cross-track.
    StateVector pushed = deputy;
    pushed.v           = deputy.v + 0.1 * Vec3{0.0, -s, c};
    expectVecNear(relativeStateRic(chief, pushed).v, {0.0, 0.0, 0.1}, 3e-12, "cross-track velocity");
    expectVecNear(relativeRicToDeputy(chief,
                                      {
                                          rel.r,
                                          {0.0, 0.0, 0.1}
    })
                      .v,
                  pushed.v,
                  3e-12,
                  "back pushed v");
}

TEST(OrbitalFrames, CovarianceRotatedIntoRicIsKnownByHand) {
    // The RIC matrix of HandState is a rotation by i about x: rows (1, 0, 0), (0, c, s), (0, -s, c).
    // A covariance with a position block [[A, D, 0], [D, B, 0], [0, 0, C]], a velocity block
    // diag(a, b, cc) and two position-velocity terms k = cov(y, vy) and m = cov(x, vz) becomes, by
    // R P R^T worked out by hand, the matrix below. A transposed rotation flips the sign of every term
    // in s; the m terms make the position-velocity block unsymmetric, so that its transpose differs.
    const HandState h;
    const double    c = h.c, s = h.s;
    const double    A = 1.0e4, B = 2.5e3, C = 1.0e2, D = 1.5e3, a = 1e-2, b = 2.5e-3, cc = 1e-4, k = 0.3, m = -0.2;
    Mat6            P{};
    P[0][0] = A;
    P[1][1] = B;
    P[2][2] = C;
    P[0][1] = P[1][0] = D;
    P[3][3]           = a;
    P[4][4]           = b;
    P[5][5]           = cc;
    P[1][4] = P[4][1] = k;
    P[0][5] = P[5][0] = m;
    Mat6 expected{};
    // clang-format off
    const double rr[3][3] = {{A,      D * c,                  -D * s                },
                             {D * c,  c * c * B + s * s * C,  c * s * (C - B)       },
                             {-D * s, c * s * (C - B),        s * s * B + c * c * C }};
    const double vv[3][3] = {{a,      0.0,                    0.0                   },
                             {0.0,    c * c * b + s * s * cc, c * s * (cc - b)      },
                             {0.0,    c * s * (cc - b),       s * s * b + c * c * cc}};
    const double rv[3][3] = {{0.0,    m * s,                  m * c                 },
                             {0.0,    k * c * c,              -k * c * s            },
                             {0.0,    -k * c * s,             k * s * s             }};
    // clang-format on
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            expected[i][j]         = rr[i][j];
            expected[i + 3][j + 3] = vv[i][j];
            expected[i][j + 3]     = rv[i][j];
            expected[i + 3][j]     = rv[j][i];
        }
    }
    const Mat6 Q = rotateCovariance(P, ricMatrix(h.state));
    for (std::size_t i = 0; i < 6; ++i) {
        for (std::size_t j = 0; j < 6; ++j) {
            // Rounding: a few units in the last place of the largest term, 1e4 m^2.
            EXPECT_NEAR(Q[i][j], expected[i][j], 1e-11) << i << ", " << j;
        }
    }
}

TEST(Eclipse, GrazingPassThroughThePenumbraOnlyIsAPenumbralEclipse) {
    // A circular orbit with the Sun raised out of its plane by beta: the satellite comes nearest the
    // shadow axis at u = 180 deg. With beta between the angle at which it would just touch the umbra
    // (c = b - a there) and the angle at which it would just touch the penumbra (c = a + b), both
    // found by bisection, the eclipse is penumbral only: documented as umbral = false, with the umbra
    // times equal to the penumbra times. (From 500 km the Earth's disk has a radius of 67 deg, so beta
    // is about 1.2 rad.)
    const double      r      = constants::wgs84A + 500 * km;
    const StateVector s0     = elementsToState({r, 0.0, 0.0, 0.0, 0.0, 0.0});
    const double      period = orbitalPeriod(r);
    const Time        t0;
    auto              sat   = [&](const Time& t) { return propagateKepler(s0, t - t0).r; };
    const auto        sunAt = [](double beta) { return constants::au * Vec3{std::cos(beta), 0.0, std::sin(beta)}; };
    const Vec3        behind{-r, 0.0, 0.0};
    const double      betaUmbra = bisect(
        [&](double beta) {
            const Disks d = disksSeenFrom(behind, sunAt(beta));
            return d.c - (d.b - d.a);
        },
        0.0,
        constants::halfPi);
    const double betaPenumbra = bisect(
        [&](double beta) {
            const Disks d = disksSeenFrom(behind, sunAt(beta));
            return d.c - (d.a + d.b);
        },
        0.0,
        constants::halfPi);
    ASSERT_LT(betaUmbra, betaPenumbra);
    const double beta = 0.5 * (betaUmbra + betaPenumbra);
    auto         sun  = [&](const Time&) { return sunAt(beta); };
    ASSERT_EQ(earthShadow(behind, sunAt(beta)).state, ShadowState::Penumbra);
    const auto ecl = findEclipses(sat, sun, t0, t0 + period, 20.0);
    ASSERT_EQ(ecl.size(), 1u);
    EXPECT_FALSE(ecl[0].umbral);
    EXPECT_EQ(ecl[0].umbraEntry, ecl[0].penumbraEntry);
    EXPECT_EQ(ecl[0].umbraExit, ecl[0].penumbraExit);
    // Centred on u = 180 deg (Brent's 1 ms on each end).
    EXPECT_NEAR(0.5 * ((ecl[0].penumbraEntry - t0) + (ecl[0].penumbraExit - t0)), 0.5 * period, 2e-3);
}

TEST(Ground, LineOfSightIsBlockedWhenAnEndIsInsideTheEllipsoid) {
    // The segment meets the solid ellipsoid as soon as one end is inside it (or both are).
    const double a = constants::wgs84A;
    EXPECT_FALSE(lineOfSight({0.5 * a, 0, 0}, {2.0 * a, 0, 0}));
    EXPECT_FALSE(lineOfSight({2.0 * a, 0, 0}, {0.5 * a, 0, 0}));
    EXPECT_FALSE(lineOfSight({0.1 * a, 0, 0}, {0, 0.1 * a, 0}));
    // Just above the pole, along the polar radius b < a: clear 1 m above it, blocked 1 m below.
    const double b = Ellipsoid::wgs84().b();
    EXPECT_TRUE(lineOfSight({-1000.0, 0, b + 1.0}, {1000.0, 0, b + 1.0}));
    EXPECT_FALSE(lineOfSight({-1000.0, 0, b - 1.0}, {1000.0, 0, b - 1.0}));
}
