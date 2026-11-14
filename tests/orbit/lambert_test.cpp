/**********************************************************************
 * lambert_test.cpp
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

// Lambert's problem (lambertIzzo, lambert): textbook examples, classical closed forms that do not
// use the solver's variables (Euler's parabolic time, Lagrange's minimum-energy time), closure by
// propagation, symmetry, and the cases an earlier implementation got wrong.

#include <gtest/gtest.h>

#include <cmath>
#include <random>
#include <vector>

#include "astro/core/angle.h"
#include "astro/core/error.h"
#include "astro/orbit/elements.h"
#include "astro/orbit/two_body.h"
#include "test_util.h"

using namespace astro;
using astro::test::expectVecNear;

namespace {

constexpr double km = 1000.0;

/// r turned by `angle` about +z.
Vec3 turned(const Vec3& r, double angle) {
    return {std::cos(angle) * r.x - std::sin(angle) * r.y, std::sin(angle) * r.x + std::cos(angle) * r.y, r.z};
}

bool isFinite(const LambertSolution& s) {
    return std::isfinite(s.v1.x + s.v1.y + s.v1.z + s.v2.x + s.v2.y + s.v2.z);
}

/// Semi-major axis of the orbit through r with velocity v (negative for a hyperbola).
double semiMajorAxis(const Vec3& r, const Vec3& v, double mu) {
    return -0.5 * mu / (0.5 * norm2(v) - mu / norm(r));
}

/// Complete revolutions of an elliptic transfer, from Kepler's equation at both ends.
int revolutionsOf(const Vec3& r1, const Vec3& r2, const LambertSolution& s, double tof, double mu) {
    const KeplerianElements a = stateToElements({r1, s.v1}, mu), b = stateToElements({r2, s.v2}, mu);
    const double            n  = std::sqrt(mu / (a.a * a.a * a.a));
    const double            dM = wrapTwoPi(b.meanAnomaly() - a.meanAnomaly());
    return static_cast<int>(std::lround((n * tof - dM) / constants::twoPi));
}

} // namespace

TEST(Lambert, CurtisExample5_2) {
    // Curtis, Orbital Mechanics for Engineering Students, Example 5.2 (mu = 398600 km^3/s^2).
    // Printed to five significant figures: 0.05 m/s.
    const double          mu = 398600e9;
    const Vec3            r1{5000 * km, 10000 * km, 2100 * km}, r2{-14600 * km, 2500 * km, 7000 * km};
    const LambertSolution s = lambert(r1, r2, 3600.0, mu);
    expectVecNear(s.v1, {-5.9925 * km, 1.9254 * km, 3.2456 * km}, 0.06);
    expectVecNear(s.v2, {-3.3125 * km, -4.1966 * km, -0.38529 * km}, 0.06);
}

TEST(Lambert, ValladoExample7_5) {
    // Vallado, Example 7-5 (mu = 398600.4418 km^3/s^2). Printed to 1 mm/s.
    const double          mu = 398600.4418e9;
    const Vec3            r1{15945.34 * km, 0, 0}, r2{12214.83899 * km, 10249.46731 * km, 0};
    const LambertSolution s = lambert(r1, r2, 76.0 * 60.0, mu);
    expectVecNear(s.v1, {2.058913 * km, 2.915965 * km, 0}, 2e-3);
    expectVecNear(s.v2, {-3.451565 * km, 0.910315 * km, 0}, 2e-3);
}

TEST(Lambert, ParabolicTimeGivesAParabola) {
    // Euler's equation: the transfer is a parabola when
    //   6 sqrt(mu) t = (r1 + r2 + c)^(3/2) -+ (r1 + r2 - c)^(3/2)   (- the short way, + the long way).
    const double mu = constants::gmEarth;
    for (double angle : {0.3, 1.2, 2.5, 3.0}) {
        for (bool prograde : {true, false}) { // r1 x r2 points up: prograde is the short way
            const Vec3   r1{9000 * km, 0, 2000 * km};
            const Vec3   r2 = 1.7 * turned(r1, angle);
            const double R1 = norm(r1), R2 = norm(r2), c = norm(r2 - r1);
            const double sign = prograde ? -1.0 : 1.0;
            const double t = (std::pow(R1 + R2 + c, 1.5) + sign * std::pow(R1 + R2 - c, 1.5)) / (6.0 * std::sqrt(mu));
            const LambertSolution s = lambert(r1, r2, t, mu, prograde);
            EXPECT_NEAR(0.5 * norm2(s.v1) - mu / R1, 0.0, 1e-12 * mu / R1) << angle << " " << prograde;
            EXPECT_NEAR(0.5 * norm2(s.v2) - mu / R2, 0.0, 1e-12 * mu / R2) << angle << " " << prograde;
            EXPECT_EQ(cross(r1, s.v1).z > 0.0, prograde);
        }
    }
}

TEST(Lambert, MinimumEnergyTimeGivesTheMinimumEnergyEllipse) {
    // Lagrange's form of the time of flight at the smallest possible semi-major axis, a = s / 2
    // (alpha = pi): sqrt(mu) t = a^(3/2) (pi - (beta - sin beta)), sin(beta / 2) = sqrt((s - c) / s),
    // beta negative the long way round.
    const double mu = constants::gmEarth;
    for (double angle : {0.4, 1.5, 2.8}) {
        for (bool prograde : {true, false}) {
            const Vec3            r1{7500 * km, 1000 * km, -500 * km};
            const Vec3            r2 = 2.2 * turned(r1, angle);
            const double          R1 = norm(r1), R2 = norm(r2), c = norm(r2 - r1), s = 0.5 * (R1 + R2 + c);
            const double          a    = 0.5 * s;
            const double          beta = (prograde ? 1.0 : -1.0) * 2.0 * std::asin(std::sqrt((s - c) / s));
            const double          t    = std::pow(a, 1.5) * (constants::pi - (beta - std::sin(beta))) / std::sqrt(mu);
            const LambertSolution sol  = lambert(r1, r2, t, mu, prograde);
            EXPECT_NEAR(semiMajorAxis(r1, sol.v1, mu), a, 1e-11 * a) << angle << " " << prograde;
            EXPECT_NEAR(semiMajorAxis(r2, sol.v2, mu), a, 1e-11 * a) << angle << " " << prograde;
        }
    }
}

TEST(Lambert, ClosesOnTheArrivalPointIncludingMultipleRevolutions) {
    std::mt19937                           rng(23);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    int                                    multi = 0, total = 0;
    for (int k = 0; k < 300; ++k) {
        const Vec3 r1{(7000 + 20000 * u(rng)) * km * (u(rng) - 0.5),
                      (7000 + 20000 * u(rng)) * km * (u(rng) - 0.5),
                      3000 * km * (u(rng) - 0.5)};
        const Vec3 r2{(7000 + 20000 * u(rng)) * km * (u(rng) - 0.5),
                      (7000 + 20000 * u(rng)) * km * (u(rng) - 0.5),
                      3000 * km * (u(rng) - 0.5)};
        if (norm(r1) < 6600 * km || norm(r2) < 6600 * km) {
            continue;
        }
        const double                       tof      = 1800.0 + 80000.0 * u(rng);
        const bool                         prograde = u(rng) > 0.3;
        const std::vector<LambertSolution> all      = lambertIzzo(r1, r2, tof, constants::gmEarth, prograde, 3);
        // The single-revolution solution, then pairs: left then right branch of 1, 2, ... revolutions.
        ASSERT_EQ(all.size() % 2, 1u);
        for (std::size_t i = 0; i < all.size(); ++i) {
            const LambertSolution& s = all[i];
            EXPECT_EQ(s.revolutions, static_cast<int>((i + 1) / 2));
            EXPECT_EQ(s.rightBranch, i > 0 && i % 2 == 0);
            const StateVector end = propagateKepler({r1, s.v1}, tof);
            expectVecNear(end.r, r2, 1e-9 * norm(r2));
            expectVecNear(end.v, s.v2, 1e-9 * norm(s.v2));
            const double hz = cross(r1, s.v1).z;
            if (std::abs(hz) > 1e3) {
                EXPECT_EQ(hz > 0.0, prograde);
            }
            const double a = semiMajorAxis(r1, s.v1, constants::gmEarth);
            if (a > 0.0) {
                EXPECT_EQ(revolutionsOf(r1, r2, s, tof, constants::gmEarth), s.revolutions);
            }
            if (s.rightBranch) { // the right branch is the longer-period one of its pair
                EXPECT_GT(a, semiMajorAxis(r1, all[i - 1].v1, constants::gmEarth));
            }
            multi += s.revolutions > 0;
            ++total;
        }
    }
    EXPECT_GT(multi, 100);
    EXPECT_GT(total, 400);
}

TEST(Lambert, ReversedTransferIsTheSameOrbitBackwards) {
    // From r2 back to r1 in the same time, the other way round: the velocities swap and change sign.
    const double mu = constants::gmEarth;
    const Vec3   r1{8000 * km, 500 * km, 1500 * km}, r2{-3000 * km, 12000 * km, -4000 * km};
    for (double tof : {2000.0, 9000.0, 40000.0}) {
        for (bool prograde : {true, false}) {
            const LambertSolution out = lambert(r1, r2, tof, mu, prograde), back = lambert(r2, r1, tof, mu, !prograde);
            expectVecNear(back.v1, -1.0 * out.v2, 1e-12 * norm(out.v2));
            expectVecNear(back.v2, -1.0 * out.v1, 1e-12 * norm(out.v1));
        }
    }
}

TEST(Lambert, VeryShortFlightsApproachTheStraightLine) {
    // The short way round in a time far below any orbital period: v ~ (r2 - r1) / t, bent by gravity
    // by a relative mu t^2 / r^3 or so.
    const double mu = 1.0;
    const Vec3   r1{1.0, 0.2, -0.1}, r2{0.3, 1.4, 0.5};
    for (double tof : {1e-3, 1e-5, 1e-7}) {
        const LambertSolution s    = lambert(r1, r2, tof, mu, cross(r1, r2).z > 0.0);
        const Vec3            line = (r2 - r1) / tof;
        EXPECT_LT(norm(s.v1 - line), 3.0 * tof * tof * norm(line)) << tof;
        EXPECT_LT(norm(s.v2 - line), 3.0 * tof * tof * norm(line)) << tof;
        EXPECT_GT(0.5 * norm2(s.v1) - mu / norm(r1), 0.0); // a hyperbola
    }
}

TEST(Lambert, TransferAnglesNearZeroAnd360Degrees) {
    // lambda^2 close to 1: the time-of-flight curves bend sharply near x = 0, and for lambda near -1
    // the iteration for the minimum of a multi-revolution curve leaves the domain. An earlier
    // implementation returned NaN velocities in these cases (the first two taken from a search).
    struct Case {
        double      angle, R2, tof;
        bool        prograde;
        std::size_t solutions;
    };
    const Case cases[] = {
        {1.2432111579751993e-05, 0.98116209739840898, 3.2129880145328582, false, 1 }, // no 1-revolution solution exists
        {0.00047109649668282036, 0.99905045626675448, 12.175229007713314, false, 9 },
        {1e-4,                   1.0,                 20.0,               true,  11}, // time for 9 revolutions: all 5 asked for
        {1e-4,                   1.0,                 20.0,               false, 11},
    };

    for (const Case& c : cases) {
        const Vec3                         r1{1.0, 0.0, 0.0}, r2 = c.R2 * turned(r1, c.angle);
        const std::vector<LambertSolution> all = lambertIzzo(r1, r2, c.tof, 1.0, c.prograde, 5);
        EXPECT_EQ(all.size(), c.solutions) << c.angle << " " << c.prograde;
        for (const LambertSolution& s : all) {
            ASSERT_TRUE(isFinite(s)) << c.angle << " " << c.prograde;
            const StateVector end = propagateKepler({r1, s.v1}, c.tof, 1.0);
            EXPECT_LT(norm(end.r - r2), 1e-9) << c.angle << " " << c.prograde << " " << s.revolutions;
        }
    }
}

TEST(Lambert, MultiRevolutionSolutionsAppearAtTheMinimumTime) {
    // lambda = -0.9899: one time just above the minimum of the 1-revolution curve (both solutions
    // exist) and one below it (neither does). T = 5.7919 and 4.9570 in Izzo's units; the minimum is
    // 5.7442.
    const double sinHalf = (1.0 - 0.979911) / (1.0 + 0.979911); // lambda^2 = 0.979911 for |r1| = |r2|
    const Vec3   r1{1.0, 0.0, 0.0}, r2 = turned(r1, 2.0 * std::asin(sinHalf));
    const double s      = 1.0 + sinHalf;
    const double toTime = std::sqrt(s * s * s / 2.0);
    EXPECT_EQ(lambertIzzo(r1, r2, 5.79192468439986 * toTime, 1.0, false, 5).size(), 3u);
    EXPECT_EQ(lambertIzzo(r1, r2, 4.9569764862249 * toTime, 1.0, false, 5).size(), 1u);
    for (const LambertSolution& sol : lambertIzzo(r1, r2, 5.79192468439986 * toTime, 1.0, false, 5)) {
        EXPECT_LT(norm(propagateKepler({r1, sol.v1}, 5.79192468439986 * toTime, 1.0).r - r2), 1e-10);
    }
}

TEST(Lambert, VeryLongFlights) {
    // Flights of 1e5 to 1e6 time units between points at radius ~1: enormous ellipses, x within 1e-6
    // of -1. (An earlier implementation returned NaN.) The period is no longer resolved by the
    // velocity, so the check is on the orbit, not on closure by propagation.
    struct Case {
        double angle, R2, tof;
    };
    const Case cases[] = {
        {1.1697278571229852,  3.2069986527951637, 279821.12180283264},
        {0.94004290603109109, 1.0378602623602025, 1731537.5612165609},
    };
    for (const Case& c : cases) {
        const Vec3            r1{1.0, 0.0, 0.0}, r2 = c.R2 * turned(r1, c.angle);
        const LambertSolution s = lambert(r1, r2, c.tof, 1.0);
        ASSERT_TRUE(isFinite(s));
        const double a = semiMajorAxis(r1, s.v1, 1.0);
        ASSERT_GT(a, 0.0);
        const double period = constants::twoPi * std::pow(a, 1.5);
        EXPECT_GT(c.tof, 0.5 * period); // out to the far apsis and back
        EXPECT_LT(c.tof, period);
        EXPECT_NEAR(semiMajorAxis(r2, s.v2, 1.0), a, 1e-6 * a);
        expectVecNear(cross(r1, s.v1), cross(r2, s.v2), 1e-12 * norm(cross(r1, s.v1)));
    }
}

TEST(Lambert, PlaneThroughTheZAxisTakesTheShorterWay) {
    // In a plane that contains the z axis "prograde" has no meaning: r1 x r2 has no z component. In
    // floating point that component is zero or a rounding error of either sign, and its sign used to
    // pick the long or the short way at random. Now both settings give the shorter way.
    int noisy = 0;
    for (int k = 0; k < 400; ++k) {
        const double phi = 0.1 + 0.0157 * k; // the plane's azimuth
        const double a = 0.2, b = 1.9;       // angles from the equator, within the plane
        const Vec3   r1{std::cos(phi) * std::cos(a), std::sin(phi) * std::cos(a), std::sin(a)};
        const Vec3   r2{1.3 * std::cos(phi) * std::cos(b), 1.3 * std::sin(phi) * std::cos(b), 1.3 * std::sin(b)};
        noisy += cross(r1 / norm(r1), r2 / norm(r2)).z != 0.0 ? 1 : 0;
        for (bool prograde : {true, false}) {
            const LambertSolution s = lambert(r1, r2, 2.0, 1.0, prograde);
            // The shorter way: the orbit turns from r1 towards r2, so its normal is along r1 x r2.
            EXPECT_GT(dot(cross(r1, s.v1), cross(r1, r2)), 0.0) << "phi " << phi << " prograde " << prograde;
            EXPECT_LT(norm(propagateKepler({r1, s.v1}, 2.0, 1.0).r - r2), 1e-12) << phi;
        }
    }
    EXPECT_GT(noisy, 0) << "no case here has a rounding error in the z component: the test shows nothing";
}

TEST(Lambert, Errors) {
    const Vec3   r1{7000 * km, 0, 0}, r2{0, 9000 * km, 0};
    const double nan = std::numeric_limits<double>::quiet_NaN(), inf = std::numeric_limits<double>::infinity();
    EXPECT_THROW((void)lambert(r1, r2, 0.0), Error);
    EXPECT_THROW((void)lambert(r1, r2, -10.0), Error);
    EXPECT_THROW((void)lambert(r1, r2, nan), Error);
    EXPECT_THROW((void)lambert(r1, r2, inf), Error);
    EXPECT_THROW((void)lambert(r1, r2, 3600.0, 0.0), Error);
    EXPECT_THROW((void)lambert(r1, r2, 3600.0, nan), Error);
    EXPECT_THROW((void)lambert(Vec3{}, r2, 3600.0), Error);
    EXPECT_THROW((void)lambert(r1, Vec3{nan, 0, 0}, 3600.0), Error);
    EXPECT_THROW((void)lambert(r1, 2.0 * r1, 3600.0), Error);  // same direction
    EXPECT_THROW((void)lambert(r1, -2.0 * r1, 3600.0), Error); // opposite: the plane is undefined
    EXPECT_THROW((void)lambertIzzo(r1, r2, 3600.0, constants::gmEarth, true, 3).at(7), std::out_of_range);
    EXPECT_EQ(lambertIzzo(r1, r2, 1e6, constants::gmEarth, true, -2).size(), 1u); // negative: no revolutions
}
