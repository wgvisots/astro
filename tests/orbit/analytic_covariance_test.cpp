/**********************************************************************
 * analytic_covariance_test.cpp
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

// State transition matrices and covariances of the analytic propagators: two-body
// (propagateKeplerStm, closed form) and secular J2/J4 (SecularPropagator::propagateStm,
// differenced).

#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <vector>

#include "astro/core/angle.h"
#include "astro/core/mat6.h"
#include "astro/orbit/elements.h"
#include "astro/orbit/gravity.h"
#include "astro/orbit/numerical.h"
#include "astro/orbit/secular.h"
#include "astro/orbit/two_body.h"

using namespace astro;

namespace {

constexpr double km = 1000.0;

/// Largest difference relative to the largest element of the same 3x3 block of b, over all blocks
/// (the blocks of a state transition matrix have different units).
double maxBlockRelDiff(const Mat6& a, const Mat6& b) {
    double worst = 0.0;
    for (std::size_t br = 0; br < 2; ++br) {
        for (std::size_t bc = 0; bc < 2; ++bc) {
            double num = 0.0, den = 0.0;
            for (std::size_t i = 0; i < 3; ++i) {
                for (std::size_t j = 0; j < 3; ++j) {
                    num = std::max(num, std::abs(a[3 * br + i][3 * bc + j] - b[3 * br + i][3 * bc + j]));
                    den = std::max(den, std::abs(b[3 * br + i][3 * bc + j]));
                }
            }
            worst = std::max(worst, num / den);
        }
    }
    return worst;
}

double& component(StateVector& s, std::size_t j) {
    return j < 3 ? (&s.r.x)[j] : (&s.v.x)[j - 3];
}

/// d(propagate(s0, dt)) / d(s0) by fourth-order central differences, steps delta |r0| and delta |v0|.
template <class Propagate>
Mat6 differences(Propagate&& propagate, const StateVector& s0, double delta) {
    constexpr double offset[4] = {1.0, -1.0, 2.0, -2.0};
    constexpr double weight[4] = {8.0, -8.0, -1.0, 1.0};
    Mat6             m{};
    for (std::size_t j = 0; j < 6; ++j) {
        const double h = delta * (j < 3 ? norm(s0.r) : norm(s0.v));
        for (std::size_t p = 0; p < 4; ++p) {
            StateVector s    = s0;
            component(s, j) += offset[p] * h;
            StateVector out  = propagate(s);
            for (std::size_t i = 0; i < 6; ++i) {
                m[i][j] += weight[p] * component(out, i) / (12.0 * h);
            }
        }
    }
    return m;
}

Mat6 symplecticJ() {
    Mat6 J{};
    for (std::size_t i = 0; i < 3; ++i) {
        J[i][i + 3] = 1.0;
        J[i + 3][i] = -1.0;
    }
    return J;
}

struct Case {
    const char* name;
    StateVector s0;
    double      dt;
};

std::vector<Case> twoBodyCases() {
    const StateVector leo    = elementsToState({7000 * km, 0.001, deg2rad(51.6), 0.3, 0.6, 0.9});
    const double      period = constants::twoPi * std::sqrt(std::pow(7000 * km, 3) / constants::gmEarth);
    return {
        {"LEO, 0.3 rev",        leo,                                                                1700.0        },
        {"LEO, 10.7 rev",       leo,                                                                10.7 * period },
        {"LEO, -3.2 rev",       leo,                                                                -3.2 * period },
        {"e = 0.9, 1.5 rev",    elementsToState({60000 * km, 0.9, deg2rad(30.0), 1.0, 2.0, 0.1}),   1.5 * 146000.0},
        {"hyperbola",           elementsToState({-20000 * km, 1.5, deg2rad(30.0), 1.0, 2.0, -0.5}), 20000.0       },
        {"circular equatorial", elementsToState({42164 * km, 0.0, 0.0, 0.0, 0.0, 0.0}),             2.3 * 86164.0 },
    };
}

Mat6 diagonalCovariance(double positionSigma, double velocitySigma) {
    Mat6 P{};
    for (std::size_t i = 0; i < 3; ++i) {
        P[i][i]         = positionSigma * positionSigma;
        P[i + 3][i + 3] = velocitySigma * velocitySigma;
    }
    return P;
}

} // namespace

TEST(TwoBodyStm, MatchesFiniteDifferences) {
    for (const Case& c : twoBodyCases()) {
        const Mat6 fd = differences([&](const StateVector& s) { return propagateKepler(s, c.dt); }, c.s0, 1e-6);
        EXPECT_LT(maxBlockRelDiff(propagateKeplerStm(c.s0, c.dt).stm, fd), 1e-8) << c.name;
    }
    // Nearly parabolic (alpha ~ 1e-20 /m): the differences themselves are only good to ~1e-7.
    const StateVector parabola{
        {7000 * km, 0.0,                                                               0.0  },
        {0.0,       std::sqrt(2.0 * constants::gmEarth / (7000 * km)) * (1.0 + 1e-13), 300.0}
    };
    const Mat6 fd = differences([&](const StateVector& s) { return propagateKepler(s, 5000.0); }, parabola, 1e-6);
    EXPECT_LT(maxBlockRelDiff(propagateKeplerStm(parabola, 5000.0).stm, fd), 1e-6);
}

TEST(TwoBodyStm, IsSymplectic) {
    // Two-body motion is Hamiltonian: Phi^T J Phi = J. Phi_rv reaches 3e4 s here, so the products
    // cancel large terms.
    const Mat6 J = symplecticJ();
    for (const Case& c : twoBodyCases()) {
        const Mat6 phi = propagateKeplerStm(c.s0, c.dt).stm;
        const Mat6 s   = multiply(transpose(phi), multiply(J, phi));
        for (std::size_t i = 0; i < 6; ++i) {
            for (std::size_t j = 0; j < 6; ++j) {
                EXPECT_NEAR(s[i][j], J[i][j], 1e-8) << c.name;
            }
        }
    }
}

TEST(TwoBodyStm, ComposesAlongTheTrajectory) {
    // Phi(t1 + t2, 0) = Phi(t1 + t2, t1) Phi(t1, 0).
    for (const Case& c : twoBodyCases()) {
        const StateTransition first  = propagateKeplerStm(c.s0, 0.4 * c.dt);
        const StateTransition second = propagateKeplerStm(first.state, 0.6 * c.dt);
        EXPECT_LT(maxBlockRelDiff(multiply(second.stm, first.stm), propagateKeplerStm(c.s0, c.dt).stm), 1e-9) << c.name;
    }
}

TEST(TwoBodyStm, MatchesIntegratedVariationalEquations) {
    // A point-mass NumericalPropagator integrates Phi from the gravity gradient.
    AdaptiveOptions tight;
    tight.relTol = 1e-14;
    tight.absTol = 1e-11;
    for (const Case& c : twoBodyCases()) {
        if (c.dt < 0.0) {
            continue;
        }
        const NumericalPropagator prop(Time(), c.s0, {std::make_shared<GravityField>(std::vector<double>{})}, tight);
        std::vector<StateVector>  st(1);
        std::vector<Mat6>         phi(1);
        prop.propagateGridStm(c.dt, 1.0, st, phi);
        EXPECT_LT(maxBlockRelDiff(propagateKeplerStm(c.s0, c.dt).stm, phi[0]), 1e-7) << c.name;
    }
}

TEST(TwoBodyStm, StateAndCovariance) {
    const Mat6 P0 = diagonalCovariance(10.0, 0.01);
    for (const Case& c : twoBodyCases()) {
        // The state is propagateKepler's, bit for bit, and P = Phi P0 Phi^T (symmetric).
        const StateTransition t = propagateKeplerStm(c.s0, c.dt);
        EXPECT_EQ(t.state, propagateKepler(c.s0, c.dt)) << c.name;
        const StateCovariance pc = propagateKeplerWithCovariance(c.s0, P0, c.dt);
        EXPECT_EQ(pc.state, t.state);
        const Mat6 expected = congruence(t.stm, P0);
        for (std::size_t i = 0; i < 6; ++i) {
            for (std::size_t j = 0; j < 6; ++j) {
                EXPECT_EQ(pc.covariance[i][j], expected[i][j]);
                EXPECT_EQ(pc.covariance[i][j], pc.covariance[j][i]);
            }
        }
    }
    // dt = 0: the identity.
    const StateTransition zero = propagateKeplerStm(twoBodyCases()[0].s0, 0.0);
    EXPECT_EQ(zero.stm, identity6());
    EXPECT_EQ(propagateKeplerWithCovariance(twoBodyCases()[0].s0, P0, 0.0).covariance, P0);
}

TEST(SecularStm, ReducesToTwoBodyWithoutZonals) {
    // With J2 = J4 = 0 the secular map is the two-body map, and both velocities coincide.
    ZonalParameters none;
    none.J2            = 0.0;
    none.J4            = 0.0;
    const Case cases[] = {
        {"LEO",                            elementsToState({7000 * km, 0.001, deg2rad(51.6), 0.3, 0.6, 0.9}), 86400.0    },
        {"e = 0.74",                       elementsToState({26560 * km, 0.74, deg2rad(63.4), 0.3, 4.7, 0.2}), 86400.0    },
        {"near-circular, near-equatorial", elementsToState({7000 * km, 1e-6, deg2rad(1e-5), 0.0, 0.0, 0.9}),  86400.0    },
        {"circular equatorial, a week",    elementsToState({7000 * km, 0.0, 0.0, 0.0, 0.0, 0.9}),             7 * 86400.0},
    };
    for (const Case& c : cases) {
        const Mat6 kepler = propagateKeplerStm(c.s0, c.dt).stm;
        EXPECT_LT(maxBlockRelDiff(SecularPropagator(c.s0, none).propagateStm(c.dt).stm, kepler), 1e-6) << c.name;
        EXPECT_LT(maxBlockRelDiff(SecularPropagatorCorrected(c.s0, none).propagateStm(c.dt).stm, kepler), 1e-6)
            << c.name;
    }
}

TEST(SecularStm, IsTheLinearizationOfTheSecularMap) {
    // Phi times a small initial offset (10 m, 1 cm/s) against the difference of two propagations, over
    // a day: equal to the second-order term, (delta a / a) n dt ~ 3e-4 relative.
    const StateVector s0 = elementsToState({7000 * km, 0.001, deg2rad(51.6), 0.3, 0.6, 0.9});
    const StateVector ds{
        Vec3{6.0,   -5.0,  6.0   },
        Vec3{0.006, 0.005, -0.006}
    };
    const StateVector                s1{s0.r + ds.r, s0.v + ds.v};
    const SecularPropagatorCorrected prop(s0);
    for (double dt : {3000.0, 86400.0}) {
        const StateTransition t     = prop.propagateStm(dt);
        const StateVector     moved = SecularPropagatorCorrected(s1).propagate(dt);
        const double          d[6]  = {ds.r.x, ds.r.y, ds.r.z, ds.v.x, ds.v.y, ds.v.z};
        Vec3                  lin, linV;
        for (std::size_t j = 0; j < 6; ++j) {
            lin  += Vec3{t.stm[0][j], t.stm[1][j], t.stm[2][j]} * d[j];
            linV += Vec3{t.stm[3][j], t.stm[4][j], t.stm[5][j]} * d[j];
        }
        const Vec3 dr = moved.r - t.state.r, dv = moved.v - t.state.v;
        EXPECT_LT(norm(lin - dr), 1e-3 * norm(dr)) << dt;
        EXPECT_LT(norm(linV - dv), 1e-3 * norm(dv)) << dt;
    }
}

TEST(SecularStm, VelocityConventions) {
    // The two classes give identical positions, so identical position rows of Phi. At dt = 0 the
    // original reproduces the input state (Phi = I); the corrected velocity is not the input one.
    const StateVector                s0 = elementsToState({7178 * km, 0.01, deg2rad(98.6), 1.3, 2.5, 0.2});
    const SecularPropagator          original(s0);
    const SecularPropagatorCorrected corrected(s0);
    for (double dt : {0.0, 5000.0, 86400.0}) {
        const Mat6 a = original.propagateStm(dt).stm, b = corrected.propagateStm(dt).stm;
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 6; ++j) {
                EXPECT_EQ(a[i][j], b[i][j]) << dt;
            }
        }
    }
    const Mat6 o0 = original.propagateStm(0.0).stm; // the identity, to the differences' rounding
    for (std::size_t i = 0; i < 6; ++i) {
        for (std::size_t j = 0; j < 6; ++j) {
            EXPECT_NEAR(o0[i][j], i == j ? 1.0 : 0.0, 1e-7) << i << " " << j;
        }
    }
    const Mat6 c0 = corrected.propagateStm(0.0).stm;
    EXPECT_GT(std::abs(c0[3][3] - 1.0) + std::abs(c0[4][4] - 1.0) + std::abs(c0[5][5] - 1.0), 1e-4);
}

TEST(SecularStm, GridMatchesSingleTimesAndConstructors) {
    const KeplerianElements          mean{7000 * km, 0.002, deg2rad(51.6), 0.3, 0.6, 0.9};
    const SecularPropagatorCorrected prop(mean);
    const Mat6                       P0 = diagonalCovariance(10.0, 0.01);

    // From mean elements, the covariance refers to their state.
    EXPECT_EQ(prop.epochState(), elementsToState(mean));
    const StateVector s0 = prop.epochState();

    const std::vector<StateCovariance> grid = prop.propagateGridWithCovariance(0.0, 600.0, 145, P0);
    for (std::size_t k : {0, 1, 70, 144}) {
        const StateCovariance one = prop.propagateWithCovariance(600.0 * static_cast<double>(k), P0);
        EXPECT_LT(norm(one.state.r - grid[k].state.r), 1e-6) << k;
        for (std::size_t i = 0; i < 6; ++i) {
            for (std::size_t j = 0; j < 6; ++j) {
                const double scale = std::sqrt(one.covariance[i][i] * one.covariance[j][j]);
                EXPECT_NEAR(grid[k].covariance[i][j], one.covariance[i][j], 1e-8 * scale) << k;
            }
        }
    }

    // A state start and a ModelOnly 12-component start from the same position and velocity.
    const StateTransition fromState = SecularPropagatorCorrected(s0).propagateStm(4000.0);
    const StateTransition from12 = SecularPropagatorCorrected(StateVector12(s0.r, s0.v, Vec3{1e-3, 0.0, 0.0}, Vec3{}),
                                                              {},
                                                              SecularOrder::J2SquaredJ4)
                                       .propagateStm(4000.0);
    EXPECT_LT(maxBlockRelDiff(from12.stm, fromState.stm), 1e-12);
    EXPECT_LT(maxBlockRelDiff(prop.propagateStm(4000.0).stm, fromState.stm), 1e-6);
}

TEST(SecularStm, UnmodeledForceAddsTheDerivativeOfTheResidualTerms) {
    // With HigherOrderMode::UnmodeledForce the trajectory is the secular one plus
    //   dr = da0 dt^2/2 + dj0 dt^3/6,  dv = da0 dt + dj0 dt^2/2,
    // where da0 = a_in - a(r0) and dj0 = j_in - G(r0) v0 are the residuals against the model (point
    // mass + J2 + J4, whose jerk along the motion is its gradient G times the velocity). propagateStm
    // holds a_in and j_in fixed, so Phi is the ModelOnly Phi plus the derivative of (dr, dv), written
    // out here:
    //   d(dr)/dr0 = -G dt^2/2 - J dt^3/6,  d(dr)/dv0 = -G dt^3/6,
    //   d(dv)/dr0 = -G dt     - J dt^2/2,  d(dv)/dv0 = -G dt^2/2,
    // with J = d(G(r) v0)/dr at r0 (by a central difference of the exact gradient; it is a small
    // term). Compared block by block, relative to the largest entry of each block of the added term.
    // Measured 6e-10 (the secular parts of the two differences cancel).
    const std::vector<double> J{0.0, 0.0, constants::J2, 0.0, constants::J4};
    const StateVector         s0       = elementsToState({7000 * km, 0.001, deg2rad(51.6), 0.3, 0.6, 0.9});
    const auto                gradient = [&](const Vec3& r) { return pointMassGradient(r) + zonalGradient(r, J); };
    const Mat3                G        = gradient(s0.r);
    Mat3                      Jr;
    for (std::size_t k = 0; k < 3; ++k) {
        Vec3 d;
        d[k]         = 1.0; // m
        const Vec3 c = (gradient(s0.r + d) * s0.v - gradient(s0.r - d) * s0.v) * 0.5;
        for (std::size_t i = 0; i < 3; ++i) {
            Jr(i, k) = c[i];
        }
    }
    StateVector12 s12(s0.r, s0.v, Vec3{}, Vec3{});
    s12.a = pointMassAcceleration(s0.r) + zonalAcceleration(s0.r, J) + Vec3{2e-4, -1e-4, 3e-4};
    s12.j = G * s0.v + Vec3{1e-6, 2e-6, -1e-6};
    const SecularPropagator modelOnly(s12, {}, SecularOrder::J2SquaredJ4, HigherOrderMode::ModelOnly);
    const SecularPropagator unmodeled(s12, {}, SecularOrder::J2SquaredJ4, HigherOrderMode::UnmodeledForce);
    for (const double dt : {600.0, -3000.0}) {
        const Mat6 added = [&] {
            const Mat6 a = unmodeled.propagateStm(dt).stm, b = modelOnly.propagateStm(dt).stm;
            Mat6       d{};
            for (std::size_t i = 0; i < 6; ++i) {
                for (std::size_t j = 0; j < 6; ++j) {
                    d[i][j] = a[i][j] - b[i][j];
                }
            }
            return d;
        }();
        const double t2 = dt * dt / 2.0, t3 = dt * dt * dt / 6.0;
        Mat6         truth{};
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 3; ++j) {
                truth[i][j]         = -G(i, j) * t2 - Jr(i, j) * t3;
                truth[i][j + 3]     = -G(i, j) * t3;
                truth[i + 3][j]     = -G(i, j) * dt - Jr(i, j) * t2;
                truth[i + 3][j + 3] = -G(i, j) * t2;
            }
        }
        EXPECT_LT(maxBlockRelDiff(added, truth), 1e-8) << "dt = " << dt;
    }
}
