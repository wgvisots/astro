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

#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <string>
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

TEST(Gravity, PointMassAndJ2GradientsMatchTheirClosedForms) {
    // Point mass: mu / |r|^3 (3 rhat rhat^T - I). J2: with a_i = k x_i (5 z^2 / r^7 - c_i / r^5),
    // k = (3/2) J2 mu Re^2 and c = (1, 1, 3), differentiating term by term gives
    //   da_i/dx_j = k [(5 z^2 / r^7 - c_i / r^5) delta_ij
    //                  + x_i (10 z delta_jz / r^7 - 35 z^2 x_j / r^9 + 5 c_i x_j / r^7)].
    // Both in extended precision. Allowed: 16 rounding errors of the largest entry (measured 3.6).
    using L      = long double;
    const L mu   = constants::gmEarth;
    const L k    = 1.5L * constants::J2 * mu * constants::earthRadiusEgm2008 * constants::earthRadiusEgm2008;
    const L c[3] = {1.0L, 1.0L, 3.0L};
    for (const Vec3 r : {
             Vec3{7000 * km,  0,          0         },
             Vec3{3000 * km,  -2000 * km, 6000 * km },
             Vec3{1,          -2,         7100 * km },
             Vec3{-4100 * km, 5200 * km,  -2500 * km}
    }) {
        const L x[3] = {r.x, r.y, r.z};
        const L rm = std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]), z = x[2];
        Mat3    Gp, Gz;
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 3; ++j) {
                const L d = i == j ? 1.0L : 0.0L, dz = j == 2 ? 1.0L : 0.0L;
                Gp(i, j) = static_cast<double>(mu / (rm * rm * rm) * (3.0L * x[i] * x[j] / (rm * rm) - d));
                Gz(i, j) = static_cast<double>(
                    k
                    * ((5.0L * z * z / std::pow(rm, 7) - c[i] / std::pow(rm, 5)) * d
                       + x[i]
                             * (10.0L * z * dz / std::pow(rm, 7) - 35.0L * z * z * x[j] / std::pow(rm, 9)
                                + 5.0L * c[i] * x[j] / std::pow(rm, 7))));
            }
        }
        const double eps = std::numeric_limits<double>::epsilon();
        EXPECT_LT(maxRelDiff(pointMassGradient(r), Gp), 16.0 * eps) << r.x << " " << r.y << " " << r.z;
        EXPECT_LT(maxRelDiff(zonalGradient(r, {0.0, 0.0, constants::J2}), Gz), 16.0 * eps)
            << r.x << " " << r.y << " " << r.z;
    }
}

TEST(Gravity, ZonalGradientIsTheDerivativeOfTheAccelerationAndSymmetric) {
    // J2 to J6 against a fourth-order central difference of zonalAcceleration (itself checked against
    // the gradient of the potential in orbit_test.cpp) with steps of 1 km: truncation (h / r)^4 ~
    // 1e-15 and rounding eps |a| / (h |G|) ~ 1e-12 of the largest entry. Measured 1.6e-12.
    const std::vector<double> J{0.0, 0.0, constants::J2, constants::J3, constants::J4, -2.3e-7, 5.4e-7};
    const auto                zonal = [&](const Vec3& x) { return zonalAcceleration(x, J); };
    for (const Vec3 r : {
             Vec3{7000 * km, 0,          0        },
             Vec3{3000 * km, -2000 * km, 6000 * km},
             Vec3{1,         -2,         7100 * km}
    }) {
        const double h  = 1000.0;
        const Mat3   G1 = numericalGradient(zonal, r, h), G2 = numericalGradient(zonal, r, 2.0 * h);
        const Mat3   fd = (G1 * 4.0 - G2) * (1.0 / 3.0); // Richardson: removes the h^2 term
        const Mat3   Gz = zonalGradient(r, J);
        EXPECT_LT(maxRelDiff(Gz, fd), 1e-11) << r.x << " " << r.y << " " << r.z;
        EXPECT_LT(maxRelDiff(Gz, Gz.transpose()), 1e-15);
        // The jerk is the gradient applied to the velocity.
        const Vec3 v{1.2 * km, 6.8 * km, 2.1 * km};
        astro::test::expectVecNear(Gz * v, zonalJerk(r, v, J), 1e-14 * norm(zonalJerk(r, v, J)));
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

namespace {

/// The orbits of the equinoctial Jacobian tests: circular, near-equatorial, eccentric, polar,
/// retrograde (not equatorial), and a mean longitude 1e-9 rad below 2 pi, so that the perturbed
/// states of the differences fall on both sides of the wrap.
std::vector<KeplerianElements> jacobianOrbits() {
    return {
        {7000 * km,  0.001, deg2rad(51.6),  0.3, 0.6, 0.9},
        {7000 * km,  0.0,   deg2rad(0.3),   0.0, 0.0, 1.0},
        {26000 * km, 0.6,   deg2rad(63.4),  1.0, 4.7, 0.5},
        {7100 * km,  0.002, deg2rad(90.0),  2.0, 0.3, 0.1},
        {7000 * km,  0.001, deg2rad(140.0), 0.5, 1.0, 2.0},
        elementsFromMeanAnomaly(8000 * km, 0.05, deg2rad(30.0), 3.0, 3.0, constants::twoPi - 6.0 - 1e-9),
    };
}

/// Analytic partial derivatives of a, h, k, p, q (rows 0 to 4; row 5 is left zero) with respect to
/// (r, v), written out from the definitions:
///   a from vis-viva:  da/dr = 2 a^2 r / |r|^3,  da/dv = 2 a^2 v / mu;
///   p = w_x / (1 + w_z), q = -w_y / (1 + w_z) with w the unit angular momentum (this is
///   tan(i/2) (sin raan, cos raan));
///   k = e . f, h = e . g with e the eccentricity vector (v x (r x v)) / mu - r / |r| and f, g the
///   equinoctial unit vectors f = (1 - p^2 + q^2, 2pq, -2p) / D, g = (2pq, 1 + p^2 - q^2, 2q) / D,
///   D = 1 + p^2 + q^2 (Broucke & Cefola 1972).
Mat6 analyticEquinoctialPartials(const StateVector& s, double mu) {
    const Vec3   r = s.r, v = s.v, hv = cross(r, v);
    const double rm = norm(r), hm = norm(hv);
    const double a = 1.0 / (2.0 / rm - norm2(v) / mu);
    const Vec3   w = hv / hm;
    const double p = w.x / (1.0 + w.z), q = -w.y / (1.0 + w.z), D = 1.0 + p * p + q * q;
    const Vec3   e = cross(v, hv) / mu - r / rm;
    const Vec3   f = Vec3{1.0 - p * p + q * q, 2.0 * p * q, -2.0 * p} / D;
    const Vec3   g = Vec3{2.0 * p * q, 1.0 + p * p - q * q, 2.0 * q} / D;
    // Derivatives of f and g with respect to p and q.
    const Vec3 fp = (Vec3{-2.0 * p, 2.0 * q, -2.0} - f * (2.0 * p)) / D;
    const Vec3 fq = (Vec3{2.0 * q, 2.0 * p, 0.0} - f * (2.0 * q)) / D;
    const Vec3 gp = (Vec3{2.0 * q, 2.0 * p, 0.0} - g * (2.0 * p)) / D;
    const Vec3 gq = (Vec3{2.0 * p, -2.0 * q, 2.0} - g * (2.0 * q)) / D;
    Mat6       T{};
    for (std::size_t j = 0; j < 6; ++j) {
        Vec3 u;
        u[j % 3] = 1.0;
        // Derivatives of the angular momentum, the eccentricity vector and a along component j.
        const Vec3   dh = j < 3 ? cross(u, v) : cross(r, u);
        const Vec3   de = j < 3 ? cross(v, cross(u, v)) / mu - (u / rm - r * (r[j] / (rm * rm * rm)))
                                : (cross(u, hv) + cross(v, cross(r, u))) / mu;
        const double da = j < 3 ? 2.0 * a * a * r[j] / (rm * rm * rm) : 2.0 * a * a * v[j - 3] / mu;
        const Vec3   dw = (dh - w * dot(w, dh)) / hm;
        const double dp = (dw.x * (1.0 + w.z) - w.x * dw.z) / ((1.0 + w.z) * (1.0 + w.z));
        const double dq = -(dw.y * (1.0 + w.z) - w.y * dw.z) / ((1.0 + w.z) * (1.0 + w.z));
        T[0][j]         = da;
        T[1][j]         = dot(de, g) + dot(e, gp * dp + gq * dq);
        T[2][j]         = dot(de, f) + dot(e, fp * dp + fq * dq);
        T[3][j]         = dp;
        T[4][j]         = dq;
    }
    return T;
}

/// J x for a state-space direction x = (dr, dv), and the sum of the magnitudes of its terms (the
/// scale of the rounding and truncation errors of the sum, which may cancel).
struct Applied {
    std::array<double, 6> value{}, scale{};
};
Applied apply(const Mat6& J, const Vec3& dr, const Vec3& dv) {
    Applied out;
    for (std::size_t i = 0; i < 6; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            out.value[i] += J[i][j] * dr[j] + J[i][j + 3] * dv[j];
            out.scale[i] += std::abs(J[i][j] * dr[j]) + std::abs(J[i][j + 3] * dv[j]);
        }
    }
    return out;
}

} // namespace

TEST(Covariance, EquinoctialJacobianMatchesTheAnalyticPartials) {
    // Rows a, h, k, p, q against the analytic partials above. The Jacobian is a central difference
    // documented to ~1e-8 relative. Each entry is compared in units of the change of the element
    // for a relative change of the state, |dX/dr| |r| or |dX/dv| |v|, relative to the largest such
    // change in its row. Measured: 2.6e-9.
    for (const KeplerianElements& el : jacobianOrbits()) {
        const StateVector s = elementsToState(el);
        const Mat6        J = stateToEquinoctialJacobian(s);
        const Mat6        T = analyticEquinoctialPartials(s, constants::gmEarth);
        for (std::size_t i = 0; i < 5; ++i) {
            double scale = 0.0, worst = 0.0;
            for (std::size_t j = 0; j < 6; ++j) {
                const double unit = j < 3 ? norm(s.r) : norm(s.v);
                scale             = std::max(scale, std::abs(T[i][j]) * unit);
                worst             = std::max(worst, std::abs(J[i][j] - T[i][j]) * unit);
            }
            EXPECT_LT(worst, 1e-8 * scale) << "element " << i << ", e = " << el.e << ", i = " << el.i;
        }
    }
}

TEST(Covariance, EquinoctialJacobianMovesTheMeanLongitudeAsTheOrbitMoves) {
    // The mean-longitude row against motions of the state whose effect on lambda is known exactly:
    //  - along the two-body flow (v, -mu r / |r|^3), lambda advances at the mean motion n and
    //    a, h, k, p, q stay constant;
    //  - turning the state about z by one radian turns the node, so lambda, by one radian; h, k and
    //    p, q turn with it (dh = k, dk = -h, dp = q, dq = -p);
    //  - turning it about its own orbit normal turns the perigee, so lambda, by one radian;
    //  - scaling r by (1 + x) and v by (1 - x / 2) scales a by (1 + x) and leaves every angle.
    // Tolerance: relative to the sum of the magnitudes of the terms of J x (along the flow the terms
    // of the a row cancel). Measured: 2.6e-8, in lambda where it is near 2 pi and so known only to a
    // rounding error of 2 pi; 1e-8 is the documented accuracy of the differences.
    const double mu = constants::gmEarth;
    for (const KeplerianElements& el : jacobianOrbits()) {
        const StateVector         s = elementsToState(el);
        const Mat6                J = stateToEquinoctialJacobian(s);
        const EquinoctialElements q = stateToEquinoctial(s);
        const Vec3                z{0.0, 0.0, 1.0};
        const Vec3                w  = cross(s.r, s.v) / norm(cross(s.r, s.v));
        const double              rm = norm(s.r);
        const double              n  = std::sqrt(mu / (q.a * q.a * q.a));
        const std::string         at = "e = " + std::to_string(el.e) + ", i = " + std::to_string(el.i);
        // clang-format off
        const struct {
            const char*           what;
            Applied               change;
            std::array<double, 6> expected;
        } motions[] = {
            {"two-body flow",     apply(J, s.v, s.r * (-mu / (rm * rm * rm))), {0.0, 0.0,  0.0,  0.0, 0.0,  n  }},
            {"turn about z",      apply(J, cross(z, s.r), cross(z, s.v)),      {0.0, q.k,  -q.h, q.q, -q.p, 1.0}},
            {"turn about normal", apply(J, cross(w, s.r), cross(w, s.v)),      {0.0, q.k,  -q.h, 0.0, 0.0,  1.0}},
            {"scale",             apply(J, s.r, s.v * -0.5),                   {q.a, 0.0,  0.0,  0.0, 0.0,  0.0}},
        };
        // clang-format on
        for (const auto& m : motions) {
            for (std::size_t i = 0; i < 6; ++i) {
                EXPECT_NEAR(m.change.value[i], m.expected[i], 1e-7 * m.change.scale[i])
                    << m.what << ", element " << i << ", " << at;
            }
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

TEST(Covariance, PropagatedCovarianceMatchesPerturbedTrajectories) {
    // The truth does not use the variational equations: write P0 = L L^T (Cholesky) and propagate
    // the initial state displaced by plus and minus each column of L; half the difference of each
    // pair is Phi times that column to second order, so P(t) = sum_k (Phi L_k) (Phi L_k)^T. The
    // displacements are one standard deviation (10 m, 1 cm/s), and the reference propagations use a
    // tight tolerance. Differences are compared in units of the standard deviations,
    // |dP_ij| / sqrt(P_ii P_jj). Measured: 1.8e-8 over 4 h.
    const StateVector         s0   = leo();
    const NumericalPropagator prop = NumericalPropagator::j2j4(Time(), s0);
    Mat6                      P0{};
    for (std::size_t i = 0; i < 3; ++i) {
        P0[i][i]         = 100.0; // (10 m)^2
        P0[i + 3][i + 3] = 1e-4;  // (1 cm/s)^2
        P0[i][i + 3]     = 0.05;  // some position/velocity correlation
        P0[i + 3][i]     = 0.05;
    }
    Mat6 L{};
    for (std::size_t j = 0; j < 6; ++j) {
        double d = P0[j][j];
        for (std::size_t k = 0; k < j; ++k) {
            d -= L[j][k] * L[j][k];
        }
        L[j][j] = std::sqrt(d);
        for (std::size_t i = j + 1; i < 6; ++i) {
            double x = P0[i][j];
            for (std::size_t k = 0; k < j; ++k) {
                x -= L[i][k] * L[j][k];
            }
            L[i][j] = x / L[j][j];
        }
    }
    AdaptiveOptions tight;
    tight.relTol            = 1e-14;
    tight.absTol            = 1e-11;
    const std::size_t count = 241;
    const double      step  = 60.0;
    const auto        field = std::make_shared<GravityField>(kJ2J4);
    const auto        pc    = prop.propagateGridWithCovariance(0.0, step, count, P0);
    std::vector<Mat6> truth(count, Mat6{});
    for (std::size_t k = 0; k < 6; ++k) {
        const Vec3 dr{L[0][k], L[1][k], L[2][k]}, dv{L[3][k], L[4][k], L[5][k]};
        const auto plus =
            NumericalPropagator(Time(), {s0.r + dr, s0.v + dv}, {field}, tight).propagateGrid(0.0, step, count);
        const auto minus =
            NumericalPropagator(Time(), {s0.r - dr, s0.v - dv}, {field}, tight).propagateGrid(0.0, step, count);
        for (std::size_t t = 0; t < count; ++t) {
            const Vec3   a = (plus[t].r - minus[t].r) * 0.5, b = (plus[t].v - minus[t].v) * 0.5;
            const double x[6] = {a.x, a.y, a.z, b.x, b.y, b.z};
            for (std::size_t i = 0; i < 6; ++i) {
                for (std::size_t j = 0; j < 6; ++j) {
                    truth[t][i][j] += x[i] * x[j];
                }
            }
        }
    }
    const auto grid  = prop.propagateGrid(0.0, step, count);
    double     worst = 0.0;
    for (std::size_t t = 0; t < count; ++t) {
        EXPECT_EQ(pc[t].state, grid[t]); // the extra equations do not change the steps
        for (std::size_t i = 0; i < 6; ++i) {
            for (std::size_t j = 0; j < 6; ++j) {
                worst = std::max(worst,
                                 std::abs(pc[t].covariance[i][j] - truth[t][i][j])
                                     / std::sqrt(truth[t][i][i] * truth[t][j][j]));
            }
        }
        EXPECT_EQ(pc[t].covariance, transpose(pc[t].covariance));
    }
    EXPECT_LT(worst, 1e-7);
    EXPECT_EQ(pc[0].covariance, P0);
    // Along-track uncertainty grows: position variance after 4 h far exceeds the initial 300 m^2.
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

TEST(Covariance, ProcessNoiseGrowsInBothDirections) {
    // Before the epoch Q is the noise accumulated over [t, 0]: the same variances as after it, with
    // the position-velocity correlation reversed. It must never shrink the covariance.
    const NumericalPropagator prop = NumericalPropagator::j2j4(Time(), leo());
    const double              q    = 1e-9;
    const double              t    = 20.0;
    const Mat6                Q    = prop.propagateGridWithCovariance(-t, 1.0, 1, Mat6{}, q)[0].covariance;
    for (std::size_t i = 0; i < 3; ++i) {
        EXPECT_NEAR(Q[i][i] / (q * t * t * t / 3.0), 1.0, 1e-3);
        EXPECT_NEAR(Q[i][i + 3] / (-q * t * t / 2.0), 1.0, 1e-3);
        EXPECT_NEAR(Q[i + 3][i + 3] / (q * t), 1.0, 1e-3);
    }

    // A grid that crosses the epoch matches propagating to each of its times directly.
    Mat6 P0{};
    for (std::size_t i = 0; i < 3; ++i) {
        P0[i][i]         = 100.0;
        P0[i + 3][i + 3] = 1e-4;
    }
    const double qLeo = 1e-8;
    const auto   grid = prop.propagateGridWithCovariance(-3600.0, 600.0, 13, P0, qLeo);
    for (std::size_t k = 0; k < grid.size(); ++k) {
        const double tk     = -3600.0 + 600.0 * static_cast<double>(k);
        const Mat6   direct = prop.propagateGridWithCovariance(tk, 1.0, 1, P0, qLeo)[0].covariance;
        // Differences in units of the standard deviations: |dP_ij| / sqrt(P_ii P_jj).
        double worst = 0.0;
        for (std::size_t i = 0; i < 6; ++i) {
            for (std::size_t j = 0; j < 6; ++j) {
                worst = std::max(worst,
                                 std::abs(grid[k].covariance[i][j] - direct[i][j])
                                     / std::sqrt(direct[i][i] * direct[j][j]));
            }
        }
        EXPECT_LT(worst, 1e-7) << "t = " << tk;
        const Mat6 noNoise = prop.propagateGridWithCovariance(tk, 1.0, 1, P0, 0.0)[0].covariance;
        for (std::size_t i = 0; i < 6; ++i) {
            EXPECT_GE(direct[i][i], noNoise[i][i]) << "t = " << tk << ", i = " << i;
        }
    }
}
