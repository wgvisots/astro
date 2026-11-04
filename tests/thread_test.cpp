/**********************************************************************
 * thread_test.cpp
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

// Shared objects used from several threads at once must give the same results as sequential use.
// Run under ThreadSanitizer (-fsanitize=thread) to check for data races as well.

#include <gtest/gtest.h>

#include <memory>
#include <optional>
#include <thread>
#include <vector>

#include "astro/core/angle.h"
#include "astro/ephemeris/spk.h"
#include "astro/orbit/elements.h"
#include "astro/orbit/gravity.h"
#include "astro/orbit/impact.h"
#include "astro/orbit/numerical.h"
#include "astro/orbit/secular.h"
#include "astro/satellite/estimate_tracker.h"
#include "test_util.h"

using namespace astro;

namespace {

constexpr int kThreads = 6;

const Time kEpoch = Time::fromIso8601("2024-02-01T00:00:00Z"); // inside the DE440s test kernel

StateVector leo() {
    return elementsToState({6778e3, 0.0005, deg2rad(51.6), 0.3, 0.5, 0.7});
}

/// Runs work(i) for i = 0 .. kThreads-1 concurrently.
template <class F>
void concurrently(F&& work) {
    std::vector<std::thread> threads;
    for (int i = 0; i < kThreads; ++i) {
        threads.emplace_back([&work, i] { work(i); });
    }
    for (auto& t : threads) {
        t.join();
    }
}

} // namespace

TEST(ThreadSafety, SharedNumericalPropagator) {
    // J2/J4 plus the Moon from the SPK kernel, so the ephemeris is shared too.
    const auto                spk  = std::make_shared<SpkEphemeris>(astro::test::dataPath("de440s_2024q1_le.bsp"));
    const auto                moon = std::make_shared<ThirdBodyForce>(Body::Moon, constants::gmMoon, spk, kEpoch);
    const NumericalPropagator prop(
        kEpoch,
        leo(),
        {std::make_shared<GravityField>(std::vector<double>{0, 0, constants::J2, 0, constants::J4}), moon});
    Mat6 P0{};
    for (std::size_t i = 0; i < 3; ++i) {
        P0[i][i]         = 100.0;
        P0[i + 3][i + 3] = 1e-4;
    }
    const auto                                refGrid = prop.propagateGrid(0.0, 60.0, 121);
    const auto                                refCov  = prop.propagateGridWithCovariance(0.0, 60.0, 121, P0, 1e-10);
    std::vector<std::vector<StateVector>>     grids(kThreads);
    std::vector<std::vector<StateCovariance>> covs(kThreads);
    std::vector<AdaptiveStats>                stats(kThreads);
    concurrently([&](int i) {
        for (int rep = 0; rep < 3; ++rep) {
            grids[i] = prop.propagateGrid(0.0, 60.0, 121, &stats[i]);
            covs[i]  = prop.propagateGridWithCovariance(0.0, 60.0, 121, P0, 1e-10);
        }
    });
    for (int i = 0; i < kThreads; ++i) {
        EXPECT_EQ(grids[i], refGrid);
        ASSERT_EQ(covs[i].size(), refCov.size());
        for (std::size_t k = 0; k < refCov.size(); ++k) {
            EXPECT_EQ(covs[i][k].state, refCov[k].state);
            EXPECT_EQ(covs[i][k].covariance, refCov[k].covariance);
        }
        EXPECT_GT(stats[i].accepted, 0);
        EXPECT_EQ(stats[i].accepted, stats[0].accepted);
    }
}

TEST(ThreadSafety, SharedSecularPropagatorAndTracker) {
    const SecularPropagator sec(leo());
    const auto              refSec = sec.propagateGrid(0.0, 10.0, 2000);

    EstimateTracker tracker;
    Mat6            P{};
    for (std::size_t i = 0; i < 3; ++i) {
        P[i][i]         = 100.0;
        P[i + 3][i + 3] = 1e-4;
    }
    (void)tracker.update({kEpoch, leo(), P});
    const auto refTrack = tracker.predictGrid(0.0, 60.0, 61);

    std::vector<std::vector<StateVector>>     secs(kThreads);
    std::vector<std::vector<StateCovariance>> tracks(kThreads);
    concurrently([&](int i) {
        secs[i]   = sec.propagateGrid(0.0, 10.0, 2000);
        tracks[i] = tracker.predictGrid(0.0, 60.0, 61); // const use of a shared tracker
    });
    for (int i = 0; i < kThreads; ++i) {
        EXPECT_EQ(secs[i], refSec);
        ASSERT_EQ(tracks[i].size(), refTrack.size());
        for (std::size_t k = 0; k < refTrack.size(); ++k) {
            EXPECT_EQ(tracks[i][k].state, refTrack[k].state);
            EXPECT_EQ(tracks[i][k].covariance, refTrack[k].covariance);
        }
    }
}

TEST(ThreadSafety, SharedRotationSeriesForImpacts) {
    // Many impact predictions in one window, sharing one ItrsRotationSeries across threads.
    const ConstantEop        eop;
    const ItrsRotationSeries rotations(kEpoch, 3 * 3600.0, eop, ImpactOptions{}.frameOptions);
    ImpactOptions            options;
    options.rotations = &rotations;
    Mat6 P{};
    for (std::size_t i = 0; i < 3; ++i) {
        P[i][i]         = 1e4;
        P[i + 3][i + 3] = 1e-2;
    }
    // Suborbital: a LEO-like state slowed to 80 % of circular speed comes down within an orbit.
    const StateVector s0  = leo();
    const StateVector sub = {s0.r, s0.v * 0.8};
    const auto        ref = predictImpact(kEpoch, sub, P, Frame::GCRS, eop, options);
    ASSERT_TRUE(ref && ref->uncertainty);

    std::vector<std::optional<Impact>> got(kThreads);
    concurrently([&](int i) {
        for (int rep = 0; rep < 20; ++rep) {
            got[i] = predictImpact(kEpoch, sub, P, Frame::GCRS, eop, options);
        }
    });
    for (int i = 0; i < kThreads; ++i) {
        ASSERT_TRUE(got[i] && got[i]->uncertainty);
        EXPECT_EQ(got[i]->flightTime, ref->flightTime);
        EXPECT_EQ(got[i]->itrs, ref->itrs);
        EXPECT_EQ(got[i]->uncertainty->ellipse.semiMajor, ref->uncertainty->ellipse.semiMajor);
    }
}
