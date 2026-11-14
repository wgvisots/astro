/**********************************************************************
 * numerical_grid_test.cpp
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

// Output grids of the numerical propagator that do not simply lead away from the epoch: grids that
// span it and grids that point towards it. Every output must be as good as integrating straight
// from the epoch to its own time.

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

#include <gtest/gtest.h>

#include "astro/core/angle.h"
#include "astro/core/constants.h"
#include "astro/core/mat6.h"
#include "astro/orbit/elements.h"
#include "astro/orbit/gravity.h"
#include "astro/orbit/numerical.h"

using namespace astro;

namespace {

constexpr double km = 1000.0;

/// Zonal coefficients J0..J4 with J3 left out: the point mass with J2 and J4.
const CentralBody kJ2J4 = CentralBody::earthJ2J4();

StateVector leo() {
    return elementsToState({7000 * km, 0.001, degToRad(51.6), 0.3, 0.6, 0.9});
}

/// The reviewer's setting: relative tolerance 1e-9, at which a day of low Earth orbit is good to metres.
NumericalPropagator loosePropagator() {
    AdaptiveOptions opt;
    opt.relTol = 1e-9;
    opt.absTol = 1e-6;
    return NumericalPropagator(Time(), leo(), {std::make_shared<GravityField>(kJ2J4)}, opt);
}

/// The same orbit integrated a hundred thousand times more tightly: the truth for the tests below.
NumericalPropagator tightPropagator() {
    AdaptiveOptions opt;
    opt.relTol = 1e-14;
    opt.absTol = 1e-11;
    return NumericalPropagator(Time(), leo(), {std::make_shared<GravityField>(kJ2J4)}, opt);
}

std::vector<double> gridTimes(double start, double step, std::size_t count) {
    std::vector<double> t(count);
    for (std::size_t k = 0; k < count; ++k) {
        t[k] = start + step * static_cast<double>(k);
    }
    return t;
}

/// Largest position difference (m) between two lists of states.
double worstDifference(const std::vector<StateVector>& a, const std::vector<StateVector>& b) {
    double worst = 0.0;
    for (std::size_t k = 0; k < a.size(); ++k) {
        worst = std::max(worst, norm(a[k].r - b[k].r));
    }
    return worst;
}

/// True if the symmetric matrix has a Cholesky factor (it is positive definite).
bool isPositiveDefinite(const Mat6& m) {
    double L[6][6] = {};
    for (std::size_t i = 0; i < 6; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            double sum = m[i][j];
            for (std::size_t k = 0; k < j; ++k) {
                sum -= L[i][k] * L[j][k];
            }
            if (i == j) {
                if (!(sum > 0.0)) {
                    return false;
                }
                L[i][i] = std::sqrt(sum);
            } else {
                L[i][j] = sum / L[j][j];
            }
        }
    }
    return true;
}

} // namespace

TEST(NumericalGridAcrossEpoch, TheEpochOutputIsTheEpochState) {
    // One day either side of the epoch, every 10 minutes: output 144 is at dt = 0. It used to come
    // from integrating a day back and a day forward again, and was 8 m from the state given.
    const NumericalPropagator prop = loosePropagator();
    const auto                grid = prop.propagateGrid(-86400.0, 600.0, 289);
    EXPECT_EQ(grid[144], leo());
    const auto backwards = prop.propagateGrid(86400.0, -600.0, 289);
    EXPECT_EQ(backwards[144], leo());
}

TEST(NumericalGridAcrossEpoch, EveryOutputIsAsGoodAsPropagatingToItDirectly) {
    // Truth: a much tighter integration. propagate(dts) integrates outwards from the epoch in each
    // direction; a grid across the epoch must not be worse. It used to integrate epoch -> start ->
    // end, three times the path to the far end, and was 31 m out there against 5 m.
    const NumericalPropagator prop  = loosePropagator();
    const NumericalPropagator truth = tightPropagator();
    for (const auto& [start, step] : {
             std::pair{-86400.0, 600.0 },
             std::pair{86400.0,  -600.0},
             std::pair{-86399.0, 599.0 }  // no output falls on the epoch
    }) {
        const std::vector<double> t      = gridTimes(start, step, 289);
        const auto                exact  = truth.propagate(t);
        const double              direct = worstDifference(prop.propagate(t), exact);
        const double              grid   = worstDifference(prop.propagateGrid(start, step, 289), exact);
        // The grid and propagate(dts) take different steps, so their errors differ in detail
        // (measured: 4.7 m against 3.6 m). The old path gave 30 m.
        EXPECT_LT(grid, 2.0 * direct) << "start " << start << ": grid " << grid << " m, direct " << direct << " m";
        EXPECT_GT(direct, 0.1) << "the tolerance is too tight for this test to show anything";
    }
}

TEST(NumericalGridTowardsEpoch, EveryOutputIsAsGoodAsPropagatingToItDirectly) {
    // A grid that starts far from the epoch and steps towards it, without reaching it. It used to
    // integrate out to the start and back, so the output nearest the epoch carried the error of the
    // whole round trip.
    const NumericalPropagator prop  = loosePropagator();
    const NumericalPropagator truth = tightPropagator();
    for (const auto& [start, step] : {
             std::pair{86400.0,  -600.0},
             std::pair{-86400.0, 600.0 }
    }) {
        const std::size_t         count = 140; // ends 3000 s from the epoch
        const std::vector<double> t     = gridTimes(start, step, count);
        const auto                grid  = prop.propagateGrid(start, step, count);
        const auto                exact = truth.propagate(t);
        // The output nearest the epoch is the one integrated to exactly: 2 cm from the truth, where
        // the round trip left 7.7 m.
        EXPECT_EQ(grid.back(), prop.propagate(t.back())) << "start " << start;
        EXPECT_LT(norm(grid.back().r - exact.back().r), 0.1) << "start " << start;
        EXPECT_LT(worstDifference(grid, exact), 2.0 * worstDifference(prop.propagate(t), exact)) << "start " << start;
    }
}

TEST(NumericalGridAcrossEpoch, CovarianceAtTheEpochIsTheInitialCovariance) {
    Mat6 P0{};
    for (std::size_t i = 0; i < 3; ++i) {
        P0[i][i]         = 100.0;
        P0[i + 3][i + 3] = 1e-4;
        P0[i][i + 3] = P0[i + 3][i] = 0.03;
    }
    const NumericalPropagator prop = loosePropagator();
    for (double q : {0.0, 1e-8}) {
        const auto grid = prop.propagateGridWithCovariance(-3600.0, 60.0, 121, P0, q);
        EXPECT_EQ(grid[60].state, leo()) << q;
        EXPECT_EQ(grid[60].covariance, P0) << q;
    }
}

TEST(NumericalGridAcrossEpoch, ProcessNoiseOnlyEverAccumulates) {
    // With no initial uncertainty the covariance is the process noise alone, which is zero at the
    // epoch and positive definite everywhere else. Integrating through the epoch formed it there as
    // the difference of two larger matrices, and outputs near the epoch had negative variances.
    const NumericalPropagator prop = loosePropagator();
    const double              q    = 1e-8;
    const std::size_t         n    = 2058;
    const double              step = 7200.0 / static_cast<double>(n - 1);
    const auto                grid = prop.propagateGridWithCovariance(-3600.0, step, n, Mat6{}, q);
    for (std::size_t k = 0; k < n; ++k) {
        const double t = -3600.0 + step * static_cast<double>(k);
        ASSERT_TRUE(isPositiveDefinite(grid[k].covariance)) << "t = " << t;
        // Close to the epoch the noise is that of free motion: variance q |t|^3 / 3 in position and
        // q |t| in velocity (to 1 % within a minute, where gravity has hardly acted on it).
        if (std::abs(t) < 60.0) {
            for (std::size_t i = 0; i < 3; ++i) {
                ASSERT_NEAR(grid[k].covariance[i][i],
                            q * std::abs(t * t * t) / 3.0,
                            0.01 * q * std::abs(t * t * t) / 3.0)
                    << "t = " << t;
                ASSERT_NEAR(grid[k].covariance[i + 3][i + 3], q * std::abs(t), 0.01 * q * std::abs(t)) << "t = " << t;
            }
        }
    }
}
