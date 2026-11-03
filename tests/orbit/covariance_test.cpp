/**********************************************************************
 * covariance_test.cpp
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
#include <vector>

#include "astro/core/angle.h"
#include "astro/core/error.h"
#include "astro/core/mat6.h"
#include "astro/ephemeris/low_precision.h"
#include "astro/orbit/elements.h"
#include "astro/orbit/gravity.h"
#include "astro/orbit/numerical.h"
#include "test_util.h"

using namespace astro;

namespace {

constexpr double km = 1000.0;

const std::vector<double> kJ2J4{0.0, 0.0, constants::J2, 0.0, constants::J4};

StateVector leo() {
    return elementsToState({7000 * km, 0.001, deg2rad(51.6), 0.3, 0.6, 0.9});
}

double maxRelDiff(const Mat3& a, const Mat3& b) {
    double num = 0.0, den = 0.0;
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            num = std::max(num, std::abs(a(i, j) - b(i, j)));
            den = std::max(den, std::abs(b(i, j)));
        }
    }
    return num / den;
}

/// Largest difference relative to the largest element of the same 3x3 block of b, over all blocks
/// (the blocks of a state transition matrix have different units).
double maxBlockRelDiff(const Mat6& a, const Mat6& b) {
    double worst = 0.0;
    for (std::size_t br = 0; br < 2; ++br) {
        for (std::size_t bc = 0; bc < 2; ++bc) {
            worst = std::max(worst, maxRelDiff(block(a, br, bc), block(b, br, bc)));
        }
    }
    return worst;
}

/// Gradient of an acceleration function by central differences.
template <class F>
Mat3 numericalGradient(F&& accel, const Vec3& r, double h) {
    Mat3 G;
    for (std::size_t j = 0; j < 3; ++j) {
        Vec3 d;
        (j == 0 ? d.x : j == 1 ? d.y : d.z) = h;
        const Vec3 c                        = (accel(r + d) - accel(r - d)) / (2.0 * h);
        G(0, j)                             = c.x;
        G(1, j)                             = c.y;
        G(2, j)                             = c.z;
    }
    return G;
}

} // namespace

TEST(Mat6, Helpers) {
    Mat6 a{};
    for (std::size_t i = 0; i < 6; ++i) {
        for (std::size_t j = 0; j < 6; ++j) {
            a[i][j] = static_cast<double>(i * 6 + j) - 10.0;
        }
    }
    EXPECT_EQ(multiply(identity6(), a), a);
    EXPECT_EQ(transpose(transpose(a)), a);
    const Mat6 c = congruence(a, identity6());
    EXPECT_EQ(c, transpose(c));
    EXPECT_EQ(c, multiply(a, transpose(a)));
    EXPECT_EQ(block(a, 1, 0)(2, 1), a[5][1]);
}

TEST(Gravity, GradientsAreExactAndSymmetric) {
    for (const Vec3 r : {
             Vec3{7000 * km, 0,          0        },
             Vec3{3000 * km, -2000 * km, 6000 * km},
             Vec3{1,         -2,         7100 * km}
    }) {
        const Mat3 Gp = pointMassGradient(r);
        EXPECT_LT(maxRelDiff(Gp, numericalGradient([](const Vec3& x) { return pointMassAcceleration(x); }, r, 1.0)),
                  1e-8);
        const std::vector<double> J{0.0, 0.0, constants::J2, constants::J3, constants::J4, -2.3e-7, 5.4e-7};
        const Mat3                Gz = zonalGradient(r, J);
        EXPECT_LT(maxRelDiff(Gz, numericalGradient([&](const Vec3& x) { return zonalAcceleration(x, J); }, r, 1.0)),
                  1e-7);
        EXPECT_LT(maxRelDiff(Gz, Gz.transpose()), 1e-13);
        // The jerk is the gradient applied to the velocity.
        const Vec3 v{1.2 * km, 6.8 * km, 2.1 * km};
        astro::test::expectVecNear(Gz * v, zonalJerk(r, v, J), 1e-12 * norm(zonalJerk(r, v, J)));
    }
}

TEST(Gravity, ForceModelPartialsMatchDefaultDifferences) {
    const StateVector  s = leo();
    Mat3               G, D, Gfd, Dfd;
    const GravityField field(kJ2J4);
    const Vec3         a = field.partials(0.0, s.r, s.v, G, D);
    EXPECT_EQ(a, field.acceleration(0.0, s.r, s.v));
    (void)field.ForceModel::partials(0.0, s.r, s.v, Gfd, Dfd);
    EXPECT_LT(maxRelDiff(G, Gfd), 1e-6);
    EXPECT_EQ(D, Mat3{});

    const auto           eph = std::make_shared<LowPrecisionEphemeris>();
    const ThirdBodyForce moon(Body::Moon, constants::gmMoon, eph, Time());
    const Vec3           am = moon.partials(100.0, s.r, s.v, G, D);
    astro::test::expectVecNear(am, moon.acceleration(100.0, s.r, s.v), 1e-20);
    (void)moon.ForceModel::partials(100.0, s.r, s.v, Gfd, Dfd);
    EXPECT_LT(maxRelDiff(G, Gfd), 1e-5);
}

TEST(Covariance, EquinoctialJacobianIsTheLinearization) {
    // J dx against the element change for small dx, including circular, near-equatorial, eccentric,
    // polar and (non-equatorial) retrograde orbits.
    const KeplerianElements orbits[] = {
        {7000 * km,  0.001, deg2rad(51.6),  0.3, 0.6, 0.9},
        {7000 * km,  0.0,   deg2rad(0.3),   0.0, 0.0, 1.0},
        {26000 * km, 0.6,   deg2rad(63.4),  1.0, 4.7, 0.5},
        {7100 * km,  0.002, deg2rad(90.0),  2.0, 0.3, 0.1},
        {7000 * km,  0.001, deg2rad(140.0), 0.5, 1.0, 2.0}
    };
    for (const auto& el : orbits) {
        const StateVector s          = elementsToState(el);
        const Mat6        J          = stateToEquinoctialJacobian(s);
        const double      dx[6]      = {0.7, -0.4, 0.2, 3e-4, -2e-4, 5e-4}; // m, m/s
        StateVector       t          = s;
        t.r                         += Vec3{dx[0], dx[1], dx[2]};
        t.v                         += Vec3{dx[3], dx[4], dx[5]};
        const EquinoctialElements a = stateToEquinoctial(s), b = stateToEquinoctial(t);
        const double d[6] = {b.a - a.a, b.h - a.h, b.k - a.k, b.p - a.p, b.q - a.q, wrapPi(b.lambda - a.lambda)};
        for (std::size_t i = 0; i < 6; ++i) {
            double lin = 0.0, scale = 0.0;
            for (std::size_t j = 0; j < 6; ++j) {
                lin   += J[i][j] * dx[j];
                scale += std::abs(J[i][j] * dx[j]);
            }
            EXPECT_NEAR(lin, d[i], 1e-5 * scale + 1e-15) << "element " << i << ", e = " << el.e << ", i = " << el.i;
        }
    }
}

TEST(Covariance, StmMatchesFiniteDifferences) {
    // Columns of Phi against central differences of two propagations with perturbed initial states.
    const StateVector s0    = leo();
    const auto        field = std::make_shared<GravityField>(kJ2J4);
    AdaptiveOptions   tight;
    tight.relTol = 1e-14;
    tight.absTol = 1e-11;
    const NumericalPropagator prop(Time(), s0, {field}, tight);
    const double              t = 3000.0;
    std::vector<StateVector>  st(1);
    std::vector<Mat6>         phi(1);
    prop.propagateGridStm(t, 1.0, st, phi);
    Mat6 fd{};
    for (std::size_t j = 0; j < 6; ++j) {
        const double h  = j < 3 ? 1.0 : 1e-3; // 1 m, 1 mm/s
        StateVector  sp = s0, sm = s0;
        double*      cp      = j < 3 ? &sp.r.x : &sp.v.x;
        double*      cm      = j < 3 ? &sm.r.x : &sm.v.x;
        cp[j % 3]           += h;
        cm[j % 3]           -= h;
        const StateVector a  = NumericalPropagator(Time(), sp, {field}, tight).propagate(t);
        const StateVector b  = NumericalPropagator(Time(), sm, {field}, tight).propagate(t);
        const double d[6] = {a.r.x - b.r.x, a.r.y - b.r.y, a.r.z - b.r.z, a.v.x - b.v.x, a.v.y - b.v.y, a.v.z - b.v.z};
        for (std::size_t i = 0; i < 6; ++i) {
            fd[i][j] = d[i] / (2.0 * h);
        }
    }
    EXPECT_LT(maxBlockRelDiff(phi[0], fd), 1e-6);
}

TEST(Covariance, StmIsSymplectic) {
    // Gravity alone is Hamiltonian: Phi^T J Phi = J with J = [[0, I], [-I, 0]]. The identity
    // cancels products of large terms (Phi_rv reaches thousands of seconds), so check the integrated
    // variational equations at exact stops with a tight tolerance; interpolation is tested below.
    AdaptiveOptions tight;
    tight.relTol = 1e-14;
    tight.absTol = 1e-11;
    const NumericalPropagator prop(Time(), leo(), {std::make_shared<GravityField>(kJ2J4)}, tight);
    Mat6                      J{};
    for (std::size_t i = 0; i < 3; ++i) {
        J[i][i + 3] = 1.0;
        J[i + 3][i] = -1.0;
    }
    for (double t : {600.0, 3600.0, 14400.0}) {
        std::vector<StateVector> st(1);
        std::vector<Mat6>        phi(1);
        prop.propagateGridStm(t, 1.0, st, phi);
        const Mat6 s = multiply(transpose(phi[0]), multiply(J, phi[0]));
        for (std::size_t i = 0; i < 6; ++i) {
            for (std::size_t j = 0; j < 6; ++j) {
                EXPECT_NEAR(s[i][j], J[i][j], 1e-8) << t;
            }
        }
    }
}

TEST(Covariance, InterpolatedStmMatchesExactStops) {
    // count = 1 integrates exactly to the requested time; dense grids interpolate between steps.
    for (double relTol : {1e-12, 1e-9}) {
        AdaptiveOptions opt;
        opt.relTol = relTol;
        opt.absTol = relTol * 1e3;
        const NumericalPropagator prop(Time(), leo(), {std::make_shared<GravityField>(kJ2J4)}, opt);
        std::vector<StateVector>  st(721);
        std::vector<Mat6>         phi(721);
        prop.propagateGridStm(0.0, 10.0, st, phi);
        double worst = 0.0;
        for (std::size_t k = 1; k < st.size(); k += 37) {
            std::vector<StateVector> s1(1);
            std::vector<Mat6>        p1(1);
            prop.propagateGridStm(10.0 * static_cast<double>(k), 1.0, s1, p1);
            worst = std::max(worst, maxBlockRelDiff(phi[k], p1[0]));
        }
        EXPECT_LT(worst, relTol < 1e-10 ? 2e-8 : 2e-6) << relTol; // measured 1.9e-9 and 4.2e-7
    }
}

TEST(Covariance, PropagatesThroughStm) {
    const NumericalPropagator prop = NumericalPropagator::j2j4(Time(), leo());
    Mat6                      P0{};
    for (std::size_t i = 0; i < 3; ++i) {
        P0[i][i]         = 100.0; // (10 m)^2
        P0[i + 3][i + 3] = 1e-4;  // (1 cm/s)^2
        P0[i][i + 3]     = 0.05;  // some position/velocity correlation
        P0[i + 3][i]     = 0.05;
    }
    const auto               pc = prop.propagateGridWithCovariance(0.0, 60.0, 241, P0);
    std::vector<StateVector> st(241);
    std::vector<Mat6>        phi(241);
    prop.propagateGridStm(0.0, 60.0, st, phi);
    const auto grid = prop.propagateGrid(0.0, 60.0, 241);
    for (std::size_t k = 0; k < pc.size(); ++k) {
        EXPECT_EQ(pc[k].state, grid[k]); // the extra equations do not change the steps
        EXPECT_EQ(pc[k].covariance, congruence(phi[k], P0));
        EXPECT_EQ(pc[k].covariance, transpose(pc[k].covariance));
    }
    EXPECT_EQ(pc[0].covariance, P0);
    // Along-track uncertainty grows: position variance after 4 h far exceeds the initial 100 m^2.
    const Mat6& P = pc.back().covariance;
    EXPECT_GT(P[0][0] + P[1][1] + P[2][2], 10.0 * 300.0);
}

TEST(Covariance, ProcessNoiseMatchesWhiteAccelerationModel) {
    // Over a short interval gravity barely couples in, and white acceleration noise of PSD q gives
    // Q = [[q t^3/3, q t^2/2], [q t^2/2, q t]] per axis.
    const NumericalPropagator prop = NumericalPropagator::j2j4(Time(), leo());
    const double              q    = 1e-9; // m^2/s^3
    const double              t    = 20.0;
    const auto                pc   = prop.propagateGridWithCovariance(0.0, t / 4, 5, Mat6{}, q);
    const Mat6&               Q    = pc.back().covariance;
    for (std::size_t i = 0; i < 3; ++i) {
        EXPECT_NEAR(Q[i][i] / (q * t * t * t / 3.0), 1.0, 1e-3);
        EXPECT_NEAR(Q[i][i + 3] / (q * t * t / 2.0), 1.0, 1e-3);
        EXPECT_NEAR(Q[i + 3][i + 3] / (q * t), 1.0, 1e-3);
    }
    EXPECT_EQ(Q, transpose(Q));
    EXPECT_THROW((void)prop.propagateGridWithCovariance(0.0, 1.0, 2, Mat6{}, -1.0), Error);
    std::vector<StateVector> st(2);
    std::vector<Mat6>        phi(3);
    EXPECT_THROW(prop.propagateGridStm(0.0, 1.0, st, phi), Error);
}
