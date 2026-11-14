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

#include <cmath>
#include <limits>
#include <memory>
#include <span>
#include <tuple>
#include <utility>
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

double maxVelocityDifference(const std::vector<StateVector>& a, const std::vector<StateVector>& b) {
    double m = 0.0;
    for (std::size_t k = 0; k < a.size(); ++k) {
        m = std::max(m, norm(a[k].v - b[k].v));
    }
    return m;
}

/// State of the mean elements `el` (mean anomaly M0 at the epoch) dt seconds later, with the angles
/// advanced at Brouwer's first-order J2 rates: written out independently of the library, in
/// extended precision.
StateVector meanElementState(const KeplerianElements& el, double M0, double dt, const ZonalParameters& zp) {
    using L   = long double;
    const L a = el.a, e = el.e, mu = zp.mu;
    const L n   = std::sqrt(mu / (a * a * a));
    const L eta = std::sqrt(1.0L - e * e);
    const L k   = 0.75L * n * zp.J2 * (zp.re / (a * eta * eta)) * (zp.re / (a * eta * eta));
    const L ci = std::cos(static_cast<L>(el.i)), si = std::sin(static_cast<L>(el.i));
    const L raan = el.raan + -2.0L * k * ci * dt;
    const L argp = el.argp + k * (5.0L * ci * ci - 1.0L) * dt;
    const L M    = M0 + (n + k * eta * (3.0L * ci * ci - 1.0L)) * dt;
    L       E    = M + e * std::sin(M);
    for (int it = 0; it < 50; ++it) {
        E -= (E - e * std::sin(E) - M) / (1.0L - e * std::cos(E));
    }
    // Perifocal position and velocity, then the rotation by argp, i and raan.
    const L     x = a * (std::cos(E) - e), y = a * eta * std::sin(E);
    const L     rdot = std::sqrt(mu * a) / (a * (1.0L - e * std::cos(E)));
    const L     vx = -rdot * std::sin(E), vy = rdot * eta * std::cos(E);
    const L     cO = std::cos(raan), sO = std::sin(raan), cw = std::cos(argp), sw = std::sin(argp);
    const L     P[3] = {cO * cw - sO * sw * ci, sO * cw + cO * sw * ci, sw * si};
    const L     Q[3] = {-cO * sw - sO * cw * ci, -sO * sw + cO * cw * ci, cw * si};
    StateVector s;
    for (std::size_t j = 0; j < 3; ++j) {
        s.r[j] = static_cast<double>(x * P[j] + y * Q[j]);
        s.v[j] = static_cast<double>(vx * P[j] + vy * Q[j]);
    }
    return s;
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

TEST(SecularPropagator, StatesAreTheTwoBodyStatesOfTheDriftingMeanElements) {
    // The truth, written out here: the node, perigee and mean anomaly advance at the first-order J2
    // rates of Brouwer (1959),
    //   raanDot = -(3/2) n J2 (Re/p)^2 cos i,  argpDot = (3/4) n J2 (Re/p)^2 (5 cos^2 i - 1),
    //   MDot    = n + (3/4) n J2 (Re/p)^2 sqrt(1 - e^2) (3 cos^2 i - 1),
    // Kepler's equation is solved by Newton's method and the state is built from the perifocal
    // vectors, all in extended precision. The velocity is the two-body velocity of the mean
    // elements (secular.h). Checked through the class (from elements and from a state) and both
    // free functions.
    const ZonalParameters zp;
    const double          eps = std::numeric_limits<double>::epsilon();
    for (const KeplerianElements& el : {
             KeplerianElements{7100 * km,  0.02,  deg2rad(63.0),  1.0, 2.0, 3.0},
             KeplerianElements{26600 * km, 0.7,   deg2rad(63.4),  4.0, 4.7, 0.1},
             KeplerianElements{7000 * km,  0.001, deg2rad(98.0),  0.2, 0.3, 5.0},
             KeplerianElements{8000 * km,  0.1,   deg2rad(140.0), 5.5, 1.2, 2.0},
    }) {
        const StateVector       s0 = elementsToState(el);
        const SecularPropagator fromElements(el, zp, SecularOrder::J2);
        const SecularPropagator fromState(s0, zp, SecularOrder::J2);
        const double            M0 = el.meanAnomaly();
        for (double dt : {0.0, 100.0, -5000.0, 86400.0, 30 * 86400.0}) {
            const StateVector truth = meanElementState(el, M0, dt, zp);
            // Rounding of the mean anomaly n dt (up to 3000 rad after 30 days) and of the angles,
            // times the size of the state, and 1 / (1 - e) for the Kepler solution near perigee.
            // Measured: 7 eps in these units (from a state, which adds the element conversion).
            const double      turns       = 1.0 + std::abs(meanMotion(el.a) * dt);
            const double      tolR        = 16.0 * eps * turns * el.a / (1.0 - el.e);
            const double      tolV        = 16.0 * eps * turns * norm(truth.v) / (1.0 - el.e);
            const StateVector viaElements = elementsToState(propagateSecular(el, dt, zp, SecularOrder::J2));
            for (const auto& [what, s] : {
                     std::pair{"class from elements", fromElements.propagate(dt)},
                     std::pair{"class from a state", fromState.propagate(dt)},
                     std::pair{"free function on a state", propagateSecular(s0, dt, zp, SecularOrder::J2)},
                     std::pair{"free function on elements", viaElements},
            }) {
                EXPECT_LT(norm(s.r - truth.r), tolR) << what << ", a = " << el.a << ", dt = " << dt;
                EXPECT_LT(norm(s.v - truth.v), tolV) << what << ", a = " << el.a << ", dt = " << dt;
            }
        }
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

TEST(SecularPropagator, GridSolvesKeplerNearPerigee) {
    // Regression: on very eccentric orbits the grid's starting guess for E (previous E plus the
    // mean-anomaly step) is far from the root near perigee, where the unguarded Chebyshev
    // iteration used to stall and return a state thousands of km off.
    for (double e : {0.75, 0.8, 0.86, 0.9, 0.95, 0.99}) {
        const SecularPropagator prop({26600 * km, e, 1.1, 0.4, 1.3, 0.0});
        for (double dM = 0.15; dM < 0.96; dM += 0.05) {
            const double step  = dM / prop.rates().meanAnomalyDot;
            const auto   grid  = prop.propagateGrid(0.0, step, 400);
            double       worst = 0.0;
            for (std::size_t k = 0; k < grid.size(); ++k) {
                worst = std::max(worst, norm(grid[k].r - prop.propagate(step * static_cast<double>(k)).r));
            }
            EXPECT_LT(worst, 1e-4) << "e = " << e << ", dM = " << dM; // metres; rounding near perigee
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

TEST(SecularPropagatorCorrected, SamePositionsAndVelocityIsTheirDerivative) {
    // Near-circular, eccentric, near-equatorial and retrograde orbits, both secular orders.
    const KeplerianElements orbits[] = {
        {7000 * km,  0.0,   deg2rad(51.6),  0.3, 0.0, 0.7},
        {7178 * km,  0.01,  deg2rad(98.6),  1.3, 2.5, 0.2},
        {16878 * km, 0.56,  deg2rad(63.4),  0.1, 4.7, 0.0},
        {7000 * km,  0.002, deg2rad(0.5),   0.0, 0.0, 1.0},
        {7000 * km,  0.001, deg2rad(140.0), 0.5, 1.0, 2.0},
    };
    for (const KeplerianElements& el : orbits) {
        for (SecularOrder order : {SecularOrder::J2, SecularOrder::J2SquaredJ4}) {
            const SecularPropagator          base(el, {}, order);
            const SecularPropagatorCorrected corrected(el, {}, order);
            EXPECT_FALSE(base.velocityIsPositionRate());
            EXPECT_TRUE(corrected.velocityIsPositionRate());
            double worstRate = 0.0, worstBase = 0.0, minBase = 1e30;
            for (double dt : {0.0, 1234.5, -5000.0, 86400.0}) {
                const StateVector s = corrected.propagate(dt);
                EXPECT_EQ(s.r, base.propagate(dt).r) << order << " dt " << dt;
                // Central difference: truncation ~|d3r/dt3| h^2 / 6 ~ 1e-4 m/s at h = 0.25 s.
                const double h  = 0.25;
                const Vec3   dr = (corrected.propagate(dt + h).r - corrected.propagate(dt - h).r) / (2.0 * h);
                worstRate       = std::max(worstRate, norm(s.v - dr));
                worstBase       = std::max(worstBase, norm(base.propagate(dt).v - dr));
                minBase         = std::min(minBase, norm(base.propagate(dt).v - dr));
            }
            EXPECT_LT(worstRate, 1e-3) << order << " a = " << el.a;
            // The original velocity misses the drift by ~1e-3 relative (metres per second).
            EXPECT_GT(minBase, 0.1) << order << " a = " << el.a;
            EXPECT_LT(worstBase, 30.0) << order << " a = " << el.a;
        }
    }
}

TEST(SecularPropagatorCorrected, EpochGridTwelveComponentsAndBaseClass) {
    const StateVector s0 = elementsToState({7000 * km, 0.001, deg2rad(97.0), 0.1, 0.2, 0.3});

    // From a state: the input position is reproduced, the velocity is not (SecularPropagator
    // reproduces both).
    const SecularPropagatorCorrected corrected(s0);
    const StateVector                atEpoch = corrected.propagate(0.0);
    expectVecNear(atEpoch.r, s0.r, 1e-6);
    EXPECT_GT(norm(atEpoch.v - s0.v), 1.0);
    EXPECT_LT(norm(atEpoch.v - s0.v), 20.0);
    expectVecNear(SecularPropagator(s0).propagate(0.0).v, s0.v, 1e-9);

    // The grid gives the corrected velocity too.
    const auto grid = corrected.propagateGrid(-300.0, 17.0, 2000);
    for (std::size_t k = 0; k < grid.size(); k += 97) {
        const StateVector direct = corrected.propagate(-300.0 + 17.0 * static_cast<double>(k));
        expectVecNear(grid[k].r, direct.r, 1e-6);
        expectVecNear(grid[k].v, direct.v, 1e-9);
    }

    // 12 components: the field's acceleration and jerk at the corrected velocity; an unmodeled
    // start reproduces the input acceleration and jerk at the epoch.
    const auto g12 = corrected.propagateGrid12(0.0, 60.0, 50);
    for (std::size_t k = 0; k < g12.size(); k += 7) {
        expectVecNear(g12[k].v, corrected.propagate(60.0 * static_cast<double>(k)).v, 1e-9);
        expectVecNear(g12[k].j, pointMassJerk(g12[k].r, g12[k].v) + zonalJerk(g12[k].r, g12[k].v, kJ2J4), 1e-15);
    }
    StateVector12 s12  = SecularPropagator(s0).propagate12(0.0);
    s12.a             += Vec3{1e-4, -2e-4, 3e-4};
    s12.j             += Vec3{1e-6, 0.0, -1e-6};
    const SecularPropagatorCorrected thrusting(s12, {}, SecularOrder::J2SquaredJ4, HigherOrderMode::UnmodeledForce);
    const StateVector12              t0 = thrusting.propagate12(0.0);
    expectVecNear(t0.a, s12.a, 1e-13);
    expectVecNear(t0.j, s12.j, 1e-15);
    expectVecNear(thrusting.accelerationResidual(), Vec3{1e-4, -2e-4, 3e-4}, 1e-15);

    // Through a base-class reference or a copy it is still corrected.
    const SecularPropagator& ref  = corrected;
    const SecularPropagator  copy = corrected;
    EXPECT_EQ(ref.propagate(500.0), corrected.propagate(500.0));
    EXPECT_EQ(copy.propagate(500.0), corrected.propagate(500.0));
    EXPECT_TRUE(copy.velocityIsPositionRate());
    EXPECT_THROW(SecularPropagatorCorrected(KeplerianElements{7000 * km, 1.2, 0.5, 0, 0, 0}), Error);
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
        const auto grid  = prop.propagateGrid(start, step, static_cast<std::size_t>(count));
        const auto stops = prop.propagate(dts);
        // Measured: 8e-7 m and 9e-10 m/s (the Hermite interpolant of the grid against integration
        // stopped at every output).
        EXPECT_LT(maxPositionDifference(grid, stops), 4e-6) << "step " << step;
        EXPECT_LT(maxVelocityDifference(grid, stops), 4e-9) << "step " << step;
        EXPECT_EQ(grid.front(), prop.propagate(start));
    }
    // The grid never takes more steps than the error control asks for, whatever the output spacing.
    AdaptiveStats gridStats, endStats;
    (void)prop.propagateGrid(0.0, 1.0, 5401, &gridStats);
    (void)prop.propagate(5400.0, &endStats);
    EXPECT_GT(gridStats.accepted, 0);
    EXPECT_LE(gridStats.accepted, endStats.accepted + 1);
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
    {
        // An empty grid reports zero work rather than leaving the caller's counts in place.
        const NumericalPropagator p =
            NumericalPropagator::j2j4(Time(), elementsToState({7000 * km, 0.001, 1.0, 0.0, 0.0, 0.0}));
        AdaptiveStats st{7, 7, 7};
        p.propagateGrid(0.0, 10.0, std::span<StateVector>{}, &st);
        EXPECT_EQ(st.evaluations, 0);
        st = {7, 7, 7};
        p.propagateGridWithCovariance(0.0, 10.0, Mat6{}, std::span<StateCovariance>{}, 1e-9, &st);
        EXPECT_EQ(st.accepted, 0);
    }
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
