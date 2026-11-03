/**********************************************************************
 * grid_test.cpp
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

#include <memory>
#include <tuple>
#include <vector>

#include "astro/core/angle.h"
#include "astro/core/error.h"
#include "astro/frames/transform.h"
#include "astro/orbit/elements.h"
#include "astro/orbit/gravity.h"
#include "astro/orbit/numerical.h"
#include "astro/orbit/secular.h"
#include "test_util.h"

using namespace astro;
using astro::test::expectVecNear;

namespace {

constexpr double km = 1000.0;

const std::vector<double> kJ2J4{0.0, 0.0, constants::J2, 0.0, constants::J4};

double maxPositionDifference(const std::vector<StateVector>& a, const std::vector<StateVector>& b) {
    double m = 0.0;
    for (std::size_t k = 0; k < a.size(); ++k) {
        m = std::max(m, norm(a[k].r - b[k].r));
    }
    return m;
}

} // namespace

TEST(Gravity, ZonalJerkIsDirectionalDerivative) {
    // The zonal acceleration depends on position only, so its jerk is the derivative along v.
    const std::vector<double> J{0.0, 0.0, constants::J2, constants::J3, constants::J4, -2.3e-7, 5.4e-7};
    for (const Vec3 r : {
             Vec3{7000 * km, 0,          0        },
             Vec3{3000 * km, -2000 * km, 6000 * km},
             Vec3{0,         1,          7100 * km}
    }) {
        const Vec3   v{1.2 * km, 6.8 * km, 2.1 * km};
        const double h  = 1e-2;
        const Vec3   fd = (zonalAcceleration(r + v * h, J) - zonalAcceleration(r - v * h, J)) / (2 * h);
        const Vec3   j  = zonalJerk(r, v, J);
        expectVecNear(j, fd, 1e-9 * norm(fd) + 1e-20);
        Vec3 a2, j2;
        zonalAccelerationAndJerk(r, v, J, a2, j2);
        expectVecNear(a2, zonalAcceleration(r, J), 1e-22);
        expectVecNear(j2, j, 1e-24);
    }
}

TEST(SecularPropagator, MatchesFreeFunctions) {
    const KeplerianElements el{7100 * km, 0.02, deg2rad(63.0), 1.0, 2.0, 3.0};
    const StateVector       s0 = elementsToState(el);
    const SecularPropagator prop(s0);
    for (double dt : {0.0, 100.0, -5000.0, 86400.0, 30 * 86400.0}) {
        const StateVector s = propagateSecular(s0, dt);
        EXPECT_EQ(prop.propagate(dt), s);
        // Same orbit through the element form (true anomaly and elementsToState), to rounding: the
        // two forms reduce the mean anomaly and build the rotation differently (1e-12 relative).
        const StateVector viaElements = elementsToState(propagateSecular(stateToElements(s0), dt));
        EXPECT_LT(norm(s.r - viaElements.r), 1e-12 * norm(s.r)) << dt;
        EXPECT_LT(norm(s.v - viaElements.v), 1e-12 * norm(s.v)) << dt;
    }
    EXPECT_THROW(SecularPropagator(KeplerianElements{7000 * km, 1.2, 0.5, 0, 0, 0}), Error);
}

TEST(SecularPropagator, GridMatchesDirectEvaluation) {
    // The grid steps its angles incrementally; it must agree with propagate(dt) at every point,
    // across mean-anomaly wraps, for near-circular to very eccentric orbits, both directions, and
    // a step large enough to disable the incremental path.
    for (double e : {0.0, 1e-4, 0.001, 0.05, 0.3, 0.7, 0.95}) {
        const KeplerianElements el{(e > 0.5 ? 26000 : 7200) * km, e, deg2rad(51.6), 0.4, 1.3, 0.2};
        const SecularPropagator prop(el);
        for (const auto& [start, step, count] : {
                 std::tuple{0.0,    10.0,   20000},
                 std::tuple{-300.0, 60.0,   3000 },
                 std::tuple{1000.0, -7.0,   5000 },
                 std::tuple{0.0,    3000.0, 200  }
        }) {
            const auto grid = prop.propagateGrid(start, step, count);
            ASSERT_EQ(grid.size(), static_cast<std::size_t>(count));
            // Relative to a (position) and |v| (velocity): near perigee of a very eccentric orbit
            // Kepler's equation amplifies rounding by 1 / (1 - e cos E) (20x at e = 0.95).
            double worstR = 0.0, worstV = 0.0;
            for (int k = 0; k < count; ++k) {
                const StateVector direct = prop.propagate(start + step * k);
                worstR                   = std::max(worstR, norm(grid[k].r - direct.r) / el.a);
                worstV                   = std::max(worstV, norm(grid[k].v - direct.v) / norm(direct.v));
            }
            EXPECT_LT(worstR, 1e-13) << "e = " << e << ", step " << step;
            EXPECT_LT(worstV, 1e-12) << "e = " << e << ", step " << step;
        }
    }
}

TEST(SecularPropagator, Grid12AddsModelDerivatives) {
    const StateVector       s0 = elementsToState({7000 * km, 0.001, deg2rad(97.0), 0.1, 0.2, 0.3});
    const SecularPropagator prop(s0);
    const auto              g = prop.propagateGrid12(0.0, 30.0, 500);
    for (std::size_t k = 0; k < g.size(); k += 50) {
        const StateVector12 d = prop.propagate12(30.0 * static_cast<double>(k));
        expectVecNear(g[k].r, d.r, 1e-6);
        expectVecNear(g[k].a, pointMassAcceleration(g[k].r) + zonalAcceleration(g[k].r, kJ2J4), 1e-13);
        expectVecNear(g[k].j, pointMassJerk(g[k].r, g[k].v) + zonalJerk(g[k].r, g[k].v, kJ2J4), 1e-15);
    }
    EXPECT_TRUE(prop.propagateGrid(0.0, 1.0, 0).empty());
}

TEST(NumericalGrid, MatchesStoppingAtEachOutput) {
    // With a tight tolerance the integrator error is far below the interpolation error being tested.
    const StateVector s0 = elementsToState({7000 * km, 0.001, deg2rad(51.6), 0.3, 0.6, 0.9});
    AdaptiveOptions   tight;
    tight.relTol = 1e-14;
    tight.absTol = 1e-11;
    const NumericalPropagator prop(Time(), s0, {std::make_shared<GravityField>(kJ2J4)}, tight);
    for (const auto& [start, step, count] : {
             std::tuple{0.0,   10.0,  1000},
             std::tuple{500.0, 37.0,  300 },
             std::tuple{0.0,   -20.0, 400 }
    }) {
        std::vector<double> dts;
        for (int k = 0; k < count; ++k) {
            dts.push_back(start + step * k);
        }
        const auto grid = prop.propagateGrid(start, step, static_cast<std::size_t>(count));
        EXPECT_LT(maxPositionDifference(grid, prop.propagate(dts)), 1e-5) << "step " << step;
        EXPECT_EQ(grid.front(), prop.propagate(start));
    }
    // The grid never takes more steps than the error control asks for, whatever the output spacing.
    (void)prop.propagateGrid(0.0, 1.0, 5401);
    const int gridSteps = prop.lastStats().accepted;
    (void)prop.propagate(5400.0);
    EXPECT_LE(gridSteps, prop.lastStats().accepted + 1);
}

TEST(NumericalGrid, DefaultToleranceInterpolationAddsLittleError) {
    // At the default tolerance the grid's error is the integrator's own error at its natural step,
    // not the interpolation's: compare against a tight-tolerance reference.
    const StateVector s0    = elementsToState({7000 * km, 0.001, deg2rad(51.6), 0.0, 0.0, 0.0});
    const auto        field = std::make_shared<GravityField>(kJ2J4);
    AdaptiveOptions   tight;
    tight.relTol  = 1e-15;
    tight.absTol  = 1e-12;
    tight.maxStep = 20.0;
    const NumericalPropagator ref(Time(), s0, {field}, tight);
    const NumericalPropagator prop(Time(), s0, {field});
    std::vector<double>       dts;
    for (int k = 0; k <= 360; ++k) {
        dts.push_back(15.0 * k);
    }
    // Endpoints of the natural steps, without interpolation, for the same arc.
    const double endpointError = norm(prop.propagate(5400.0).r - ref.propagate(5400.0).r);
    const double gridError     = maxPositionDifference(prop.propagateGrid(0.0, 15.0, 361), ref.propagate(dts));
    EXPECT_LT(gridError, std::max(3.0 * endpointError, 1e-4));
}

TEST(NumericalGrid, TwelveStateAndEdgeCases) {
    const StateVector         s0   = elementsToState({7000 * km, 0.01, deg2rad(28.5), 0.0, 1.0, 2.0});
    const NumericalPropagator prop = NumericalPropagator::j2j4(Time(), s0);
    const auto                g    = prop.propagateGrid12(0.0, 60.0, 100);
    const auto                g6   = prop.propagateGrid(0.0, 60.0, 100);
    for (std::size_t k = 0; k < g.size(); ++k) {
        EXPECT_EQ(g[k].rv(), g6[k]);
        EXPECT_EQ(g[k], prop.withDerivatives(60.0 * static_cast<double>(k), g6[k]));
    }
    EXPECT_TRUE(prop.propagateGrid(0.0, 60.0, 0).empty());
    EXPECT_EQ(prop.propagateGrid(250.0, 60.0, 1).front(), prop.propagate(250.0));
    EXPECT_THROW((void)prop.propagateGrid(0.0, 0.0, 10), Error);
}

TEST(NumericalGrid, J2J4FrameOverloadAlignsThePole) {
    // The same orbit given in GCRS and in TOD (frozen at the epoch) must propagate identically once
    // the zonal field is aligned with the true pole; ignoring the alignment in GCRS costs hundreds of
    // metres over 4 h.
    const Time        t0 = Time::fromIso8601("2026-11-02T00:00:00Z");
    const StateVector sg = elementsToState({6778 * km, 0.0005, deg2rad(51.6), 0.3, 0.5, 0.7});
    const Mat3        B  = gcrsToCelestialMatrix(Frame::TOD, t0);
    const StateVector st{B * sg.r, B * sg.v};
    const auto        g    = NumericalPropagator::j2j4(t0, sg, Frame::GCRS).propagateGrid(0.0, 600.0, 25);
    const auto        tod  = NumericalPropagator::j2j4(t0, st, Frame::TOD).propagateGrid(0.0, 600.0, 25);
    const auto        bare = NumericalPropagator::j2j4(t0, sg).propagateGrid(0.0, 600.0, 25);
    double            same = 0.0, wrong = 0.0;
    for (std::size_t k = 0; k < g.size(); ++k) {
        same  = std::max(same, norm(B.transpose() * tod[k].r - g[k].r));
        wrong = std::max(wrong, norm(bare[k].r - g[k].r));
    }
    EXPECT_LT(same, 1e-3);
    EXPECT_GT(wrong, 100.0);
    // The 12-component overload takes the same path.
    const StateVector12 s12{sg.r, sg.v, {}, {}};
    EXPECT_EQ(NumericalPropagator::j2j4(t0, s12, Frame::GCRS).propagate12(3600.0).rv(),
              NumericalPropagator::j2j4(t0, sg, Frame::GCRS).propagate(3600.0));
}
