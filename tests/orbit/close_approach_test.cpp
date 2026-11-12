/**********************************************************************
 * close_approach_test.cpp
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

#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <random>
#include <sstream>
#include <vector>

#include "astro/core/angle.h"
#include "astro/core/error.h"
#include "astro/core/rootfind.h"
#include "astro/orbit/close_approach.h"
#include "astro/orbit/elements.h"
#include "astro/orbit/gravity.h"
#include "astro/orbit/numerical.h"
#include "astro/orbit/two_body.h"
#include "test_util.h"

using namespace astro;

namespace {

constexpr double km = 1000.0;

const Time kEpoch = Time::fromIso8601("2026-11-11T00:00:00Z");

/// A minimum of the distance found by brute force: time (s after the window start) and distance.
struct Minimum {
    double s{0.0};
    double distance{0.0};
};

/// Every local minimum of the distance over [0, T], by sampling every `step` seconds and refining
/// each discrete minimum with Brent's minimizer (in a local variable, for full time precision).
std::vector<Minimum>
bruteForceMinima(const Trajectory& a, const Trajectory& b, const Time& start, double T, double step) {
    const auto          fa    = a.secondsAfter(start);
    const auto          fb    = b.secondsAfter(start);
    const auto          dist2 = [&](double s) { return norm2(fb(s).r - fa(s).r); };
    const int           n     = static_cast<int>(T / step);
    std::vector<double> f(static_cast<std::size_t>(n) + 1);
    for (int k = 0; k <= n; ++k) {
        f[static_cast<std::size_t>(k)] = dist2(k * step);
    }
    std::vector<Minimum> out;
    for (int k = 1; k < n; ++k) {
        const auto i = static_cast<std::size_t>(k);
        if (f[i] < f[i - 1] && f[i] <= f[i + 1]) {
            const double         s0 = (k - 1) * step;
            const MinimizeResult m  = brentMinimize([&](double u) { return dist2(s0 + u); }, 0.0, 2.0 * step, 1e-12);
            out.push_back({s0 + m.x, std::sqrt(m.fx)});
        }
    }
    return out;
}

/// Circular LEO at 51.6 deg, and a crossing orbit whose velocity is turned 60 deg from it, built so that at kTca the two
/// are exactly 100 m apart along the radius, perpendicular to both velocities: kTca is then an
/// exact closest approach (zero range rate) with a miss distance of exactly 100 m.
struct Conjunction {
    static constexpr double kTca  = 37000.0; // s after kEpoch
    static constexpr double kMiss = 100.0;   // m

    KeplerianElements elements{7000 * km, 0.0, deg2rad(51.6), 0.3, 0.0, 0.0};
    StateVector       atTca; // second object at kTca
    StateVector       firstAtTca;

    Conjunction() {
        firstAtTca         = propagateKepler(elementsToState(elements), kTca);
        const Vec3   r     = firstAtTca.r;
        const double speed = norm(firstAtTca.v);
        // A horizontal velocity of the same speed through the same point: turn the
        // first object's velocity 60 deg about the radius.
        const Vec3   up     = normalized(r);
        const Vec3   along  = normalized(firstAtTca.v);
        const Vec3   across = cross(up, along);
        const double turn   = deg2rad(60.0);
        const Vec3   v2     = speed * (std::cos(turn) * along + std::sin(turn) * across);
        const Vec3   n      = normalized(cross(firstAtTca.v, v2)); // radial, perpendicular to both velocities
        atTca               = {r + kMiss * n, v2};
    }
};

} // namespace

TEST(CloseApproach, FindsAnExactConstructedConjunction) {
    const Conjunction c;
    // One object as elements, the other as a state at a different epoch.
    const Trajectory first{kEpoch, c.elements};
    const Trajectory second{kEpoch + Conjunction::kTca, c.atTca};
    const Time       end = kEpoch + 86400.0;

    const CloseApproach ca = closestApproach(first, second, kEpoch, end);
    EXPECT_FALSE(ca.atBoundary);
    EXPECT_NEAR(ca.time - kEpoch, Conjunction::kTca, 1e-6);
    EXPECT_NEAR(ca.distance, Conjunction::kMiss, 1e-6);
    EXPECT_NEAR(ca.relativeSpeed, norm(c.atTca.v - c.firstAtTca.v), 1e-6);
    test::expectVecNear(ca.relativePosition(), ca.second.r - ca.first.r, 0.0);
    test::expectVecNear(ca.relativeVelocity(), ca.second.v - ca.first.v, 0.0);
    test::expectVecNear(ca.first.r, c.firstAtTca.r, 1e-2);
    test::expectVecNear(ca.second.r, c.atTca.r, 1e-2);
    EXPECT_NEAR(dot(ca.relativePosition(), ca.relativeVelocity()), 0.0, 1e-3); // zero range rate

    // The same answer through the function-of-time constructor, and with the arguments swapped.
    const Trajectory    viaFunction([&first](const Time& t) { return first.state(t); });
    const CloseApproach swapped = closestApproach(second, viaFunction, kEpoch, end);
    EXPECT_NEAR(swapped.time - ca.time, 0.0, 1e-6);
    EXPECT_NEAR(swapped.distance, ca.distance, 1e-9);
}

TEST(CloseApproach, ListsEveryMinimumAgainstBruteForce) {
    const Conjunction c;
    const Trajectory  first{kEpoch, c.elements};
    const Trajectory  second{kEpoch + Conjunction::kTca, c.atTca};
    const double      T = 86400.0;

    // The orbits cross at two antipodal points, so the objects meet every half revolution, a
    // little farther apart each time (the second orbit is slightly eccentric).
    const std::vector<CloseApproach> all   = closeApproaches(first, second, kEpoch, kEpoch + T);
    const std::vector<Minimum>       brute = bruteForceMinima(first, second, kEpoch, T, 1.0);
    ASSERT_EQ(all.size(), brute.size());
    for (std::size_t k = 0; k < all.size(); ++k) {
        EXPECT_NEAR(all[k].time - kEpoch, brute[k].s, 1e-4) << k;
        EXPECT_LE(all[k].distance, brute[k].distance + 1e-6) << k;
        EXPECT_NEAR(all[k].distance, brute[k].distance, 1e-3) << k;
        EXPECT_FALSE(all[k].atBoundary);
        if (k > 0) {
            EXPECT_GT(all[k].time - all[k - 1].time, 0.0);
        }
    }

    // A distance threshold keeps the conjunctions only.
    const std::vector<CloseApproach> close = closeApproaches(first, second, kEpoch, kEpoch + T, 10 * km);
    ASSERT_FALSE(close.empty());
    EXPECT_LT(close.size(), all.size());
    bool haveTca = false;
    for (const CloseApproach& x : close) {
        EXPECT_LE(x.distance, 10 * km);
        haveTca = haveTca || std::abs(x.time - kEpoch - Conjunction::kTca) < 1e-6;
    }
    EXPECT_TRUE(haveTca);
}

TEST(CloseApproach, RandomOrbitPairsMatchBruteForce) {
    // Elliptic orbits from LEO to GEO, any inclination, perigee above 6600 km: no minimum may be
    // missed at the default step, and the global minimum must match.
    std::mt19937                           rng(20261111);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    const auto                             randomOrbit = [&] {
        const double a  = (6700 + unit(rng) * 35500) * km;
        const double e  = unit(rng) * std::max(0.0, std::min(0.75, 1.0 - 6600 * km / a));
        const double i  = unit(rng) * constants::pi;
        const double nu = unit(rng) * constants::twoPi;
        return KeplerianElements{a, e, i, unit(rng) * constants::twoPi, unit(rng) * constants::twoPi, nu};
    };
    const double T = 86400.0;
    for (int trial = 0; trial < 12; ++trial) {
        const Trajectory                 first{kEpoch, randomOrbit()};
        const Trajectory                 second{kEpoch - 3600.0 * trial, elementsToState(randomOrbit())};
        const std::vector<CloseApproach> all   = closeApproaches(first, second, kEpoch, kEpoch + T);
        const std::vector<Minimum>       brute = bruteForceMinima(first, second, kEpoch, T, 2.0);
        ASSERT_EQ(all.size(), brute.size()) << "trial " << trial;
        for (std::size_t k = 0; k < all.size(); ++k) {
            EXPECT_NEAR(all[k].time - kEpoch, brute[k].s, 1e-3) << "trial " << trial << " minimum " << k;
            EXPECT_LE(all[k].distance, brute[k].distance + 1e-6) << "trial " << trial << " minimum " << k;
        }

        const CloseApproach best     = closestApproach(first, second, kEpoch, kEpoch + T);
        double              expected = std::min(norm(second.state(kEpoch).r - first.state(kEpoch).r),
                                   norm(second.state(kEpoch + T).r - first.state(kEpoch + T).r));
        for (const Minimum& m : brute) {
            expected = std::min(expected, m.distance);
        }
        EXPECT_LE(best.distance, expected + 1e-6) << "trial " << trial;
        EXPECT_NEAR(best.distance, expected, 1e-3) << "trial " << trial;
    }
}

TEST(CloseApproach, NumericalPropagatorThroughAFunction) {
    // A point-mass NumericalPropagator wrapped as a function gives the two-body answer.
    const Conjunction         c;
    const StateVector         s0 = elementsToState(c.elements);
    const NumericalPropagator prop{kEpoch, s0, {std::make_shared<GravityField>(std::vector<double>{})}};
    const Trajectory          numerical([&prop](const Time& t) { return prop.propagate(t); });
    const Trajectory          second{kEpoch + Conjunction::kTca, c.atTca};
    const CloseApproach       ca =
        closestApproach(numerical, second, kEpoch + Conjunction::kTca - 3000.0, kEpoch + Conjunction::kTca + 3000.0);
    EXPECT_FALSE(ca.atBoundary);
    EXPECT_NEAR(ca.time - kEpoch, Conjunction::kTca, 1e-5);
    EXPECT_NEAR(ca.distance, Conjunction::kMiss, 1e-3);
}

TEST(CloseApproach, SecularModelsAndEquinoctialElements) {
    const KeplerianElements el{7000 * km, 0.01, deg2rad(51.6), 0.3, 0.2, 0.1};
    const SecularPropagator sp(el);
    const Trajectory        secular{kEpoch, el, TrajectoryModel::SecularJ2J4};
    test::expectVecNear(secular.state(kEpoch + 5000.0).r, sp.propagate(5000.0).r, 1e-6);
    const SecularPropagator spJ2(el, {}, SecularOrder::J2);
    const Trajectory        secularJ2{kEpoch, el, TrajectoryModel::SecularJ2};
    test::expectVecNear(secularJ2.state(kEpoch - 5000.0).r, spJ2.propagate(-5000.0).r, 1e-6);
    const StateVector s0 = elementsToState(el);
    const Trajectory  secularState{kEpoch, s0, TrajectoryModel::SecularJ2J4};
    test::expectVecNear(secularState.state(kEpoch + 5000.0).r, SecularPropagator(s0).propagate(5000.0).r, 1e-6);

    // The secular models are SecularPropagatorCorrected (velocity = d(position)/dt).
    const StateVector corrected = SecularPropagatorCorrected(el).propagate(7000.0);
    EXPECT_EQ(secular.state(kEpoch + 7000.0).r, corrected.r);
    EXPECT_EQ(secular.state(kEpoch + 7000.0).v, corrected.v);

    // Equinoctial elements give the same trajectory as the Keplerian ones.
    const Trajectory kepler{kEpoch, el};
    const Trajectory equinoctial{kEpoch, keplerianToEquinoctial(el)};
    test::expectVecNear(equinoctial.state(kEpoch + 5000.0).r, kepler.state(kEpoch + 5000.0).r, 1e-4);

    // The node regression separates a J2/J4 orbit from its two-body twin; the search follows it.
    // The twins drift apart slowly, so each minimum is nearly flat: minimizing the distance pins
    // its time only to seconds, while the range rate root is sharp. Compare distances, and check
    // the range rate at each reported time. The twins start at the same point, so the window
    // starts a minute later, past that degenerate minimum.
    const Time                       start = kEpoch + 60.0;
    const std::vector<CloseApproach> all   = closeApproaches(kepler, secular, start, start + 86400.0);
    const std::vector<Minimum>       brute = bruteForceMinima(kepler, secular, start, 86400.0, 2.0);
    ASSERT_EQ(all.size(), brute.size());
    for (std::size_t k = 0; k < all.size(); ++k) {
        EXPECT_LE(all[k].distance, brute[k].distance + 1e-6) << k;
        EXPECT_NEAR(all[k].distance, brute[k].distance, 1e-6 * brute[k].distance) << k;
        const double rangeRate = dot(all[k].relativePosition(), all[k].relativeVelocity()) / all[k].distance;
        EXPECT_NEAR(rangeRate, 0.0, 1e-5 * all[k].relativeSpeed) << k;
    }
}

TEST(CloseApproach, WindowEnds) {
    const Conjunction c;
    const Trajectory  first{kEpoch, c.elements};
    const Trajectory  second{kEpoch + Conjunction::kTca, c.atTca};
    const Time        tca = kEpoch + Conjunction::kTca;

    // Still closing at the end of the window: the minimum is the end.
    const CloseApproach closing = closestApproach(first, second, tca - 60.0, tca - 10.0);
    EXPECT_TRUE(closing.atBoundary);
    EXPECT_NEAR(closing.time - (tca - 10.0), 0.0, 1e-9);
    EXPECT_LT(dot(closing.relativePosition(), closing.relativeVelocity()), 0.0);
    EXPECT_TRUE(closeApproaches(first, second, tca - 60.0, tca - 10.0).empty());

    // Already separating at the start: the minimum is the start.
    const CloseApproach separating = closestApproach(first, second, tca + 10.0, tca + 60.0);
    EXPECT_TRUE(separating.atBoundary);
    EXPECT_NEAR(separating.time - (tca + 10.0), 0.0, 1e-9);

    // A short window around the approach finds it; a zero-length one is its start.
    const CloseApproach around = closestApproach(first, second, tca - 1.0, tca + 1.0);
    EXPECT_FALSE(around.atBoundary);
    EXPECT_NEAR(around.distance, Conjunction::kMiss, 1e-6);
    const CloseApproach point = closestApproach(first, second, tca, tca);
    EXPECT_NEAR(point.time - tca, 0.0, 1e-9);
    EXPECT_NEAR(point.distance, Conjunction::kMiss, 1e-6);

    // An object and itself: zero distance everywhere, no minimum strictly inside.
    const CloseApproach self = closestApproach(first, first, kEpoch, kEpoch + 3600.0);
    EXPECT_EQ(self.distance, 0.0);
    EXPECT_TRUE(closeApproaches(first, first, kEpoch, kEpoch + 3600.0).empty());
}

TEST(CloseApproach, Errors) {
    const Conjunction c;
    const Trajectory  first{kEpoch, c.elements};
    const Trajectory  second{kEpoch, c.atTca};
    const Time        end = kEpoch + 3600.0;

    EXPECT_THROW((void)closestApproach(first, second, end, kEpoch), Error);
    CloseApproachOptions bad;
    bad.stepFraction = 0.0;
    EXPECT_THROW((void)closestApproach(first, second, kEpoch, end, bad), Error);
    bad               = {};
    bad.timeTolerance = std::nan("");
    EXPECT_THROW((void)closestApproach(first, second, kEpoch, end, bad), Error);
    EXPECT_THROW((void)closeApproaches(first, second, kEpoch, end, std::nan("")), Error);

    EXPECT_THROW(Trajectory{Trajectory::StateFunction{}}, Error);
    const Trajectory nan([](const Time&) {
        return StateVector{
            {std::nan(""), 0.0,    0.0},
            {0.0,          7500.0, 0.0}
        };
    });
    EXPECT_THROW((void)closestApproach(first, nan, kEpoch, end), Error);

    // The secular models need an elliptic orbit.
    const StateVector hyperbolic{
        {7000 * km, 0.0,     0.0},
        {0.0,       12 * km, 0.0}
    };
    EXPECT_NO_THROW(Trajectory(kEpoch, hyperbolic));
    EXPECT_THROW(Trajectory(kEpoch, hyperbolic, TrajectoryModel::SecularJ2J4), Error);
    EXPECT_THROW(Trajectory(kEpoch, hyperbolic, static_cast<TrajectoryModel>(42)), Error);

    std::ostringstream os;
    os << TrajectoryModel::SecularJ2 << ' ' << static_cast<TrajectoryModel>(42);
    EXPECT_EQ(os.str(), "TrajectoryModel::SecularJ2 TrajectoryModel::INVALID_VALUE");
}
