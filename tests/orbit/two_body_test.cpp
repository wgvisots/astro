/**********************************************************************
 * two_body_test.cpp
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

// Universal-variable two-body propagation against an extended-precision reference, on every kind of
// conic and in particular where the iteration is hardest: short steps on very eccentric orbits and
// near-parabolic fly-bys.

#include <cmath>
#include <limits>
#include <random>

#include <gtest/gtest.h>

#include "astro/core/angle.h"
#include "astro/core/constants.h"
#include "astro/core/error.h"
#include "astro/orbit/elements.h"
#include "astro/orbit/two_body.h"

using namespace astro;

namespace {

using Real = long double;

constexpr double kMu = constants::gmEarth;

struct Vec3L {
    Real x, y, z;
};
Vec3L lift(const Vec3& a) {
    return {a.x, a.y, a.z};
}
Real dotL(const Vec3L& a, const Vec3L& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
Vec3L combine(Real f, const Vec3L& a, Real g, const Vec3L& b) {
    return {f * a.x + g * b.x, f * a.y + g * b.y, f * a.z + g * b.z};
}
double distance(const Vec3& a, const Vec3L& b) {
    const Vec3L d{a.x - b.x, a.y - b.y, a.z - b.z};
    return static_cast<double>(std::sqrt(dotL(d, d)));
}

/// Stumpff functions c0..c3 in extended precision: series near zero, closed forms elsewhere.
void stumpffL(Real z, Real& c0, Real& c1, Real& c2, Real& c3) {
    if (std::abs(z) < 0.5L) {
        Real t2 = 0.5L, t3 = 1.0L / 6.0L;
        c2 = 0.0L;
        c3 = 0.0L;
        for (int j = 0; j < 30; ++j) {
            c2 += t2;
            c3 += t3;
            t2 *= -z / ((3 + 2 * j) * (4.0L + 2 * j));
            t3 *= -z / ((4 + 2 * j) * (5.0L + 2 * j));
        }
        c1 = 1.0L - z * c3;
        c0 = 1.0L - z * c2;
    } else if (z > 0.0L) {
        const Real q = std::sqrt(z);
        c0           = std::cos(q);
        c1           = std::sin(q) / q;
        c2           = (1.0L - c0) / z;
        c3           = (1.0L - c1) / z;
    } else {
        const Real q = std::sqrt(-z);
        c0           = std::cosh(q);
        c1           = std::sinh(q) / q;
        c2           = (1.0L - c0) / z;
        c3           = (1.0L - c1) / z;
    }
}

struct StateL {
    Vec3L r, v;
};

/// Two-body propagation in extended precision. Universal variables, but with the universal Kepler
/// equation r0 U1 + sigma U2 + U3 = sqrt(mu) dt solved by bisection (its left side increases with
/// the universal anomaly), so it shares neither the arithmetic nor the iteration of the library.
StateL referenceKepler(const StateVector& s, double dt, double mu) {
    const Vec3L r0 = lift(s.r), v0 = lift(s.v);
    const Real  r0m = std::sqrt(dotL(r0, r0)), rootMu = std::sqrt(static_cast<Real>(mu));
    const Real  sigma = dotL(r0, v0) / rootMu, alpha = 2.0L / r0m - dotL(v0, v0) / mu, goal = rootMu * dt;
    Real        c0 = 0, c1 = 0, c2 = 0, c3 = 0;
    const auto  kepler = [&](Real chi) {
        stumpffL(alpha * chi * chi, c0, c1, c2, c3);
        return r0m * chi * c1 + sigma * chi * chi * c2 + chi * chi * chi * c3 - goal;
    };
    Real lo = 0.0L, hi = 0.0L;
    if (dt > 0.0) {
        for (hi = 1.0L; kepler(hi) < 0.0L; hi *= 2.0L) {
        }
    } else {
        for (lo = -1.0L; kepler(lo) > 0.0L; lo *= 2.0L) {
        }
    }
    Real chi = 0.5L * (lo + hi);
    while (chi != lo && chi != hi) {
        (kepler(chi) > 0.0L ? hi : lo) = chi;
        chi                            = 0.5L * (lo + hi);
    }
    kepler(chi);
    const Real U1 = chi * c1, U2 = chi * chi * c2;
    const Real r = r0m * c0 + sigma * U1 + U2;
    const Real f = 1.0L - U2 / r0m, g = (r0m * U1 + sigma * U2) / rootMu;
    const Real fDot = -rootMu * U1 / (r * r0m), gDot = 1.0L - U2 / r;
    return {combine(f, r0, g, v0), combine(fDot, r0, gDot, v0)};
}

/// A state with perigee radius rp, eccentricity e (any conic) and true anomaly nu, in a plane
/// turned by the Euler angles (3-1-3).
StateVector conicState(double rp, double e, double nu, double node, double inc, double argp) {
    const double p = rp * (1.0 + e), h = std::sqrt(kMu * p), r = p / (1.0 + e * std::cos(nu));
    const Vec3   rPf{r * std::cos(nu), r * std::sin(nu), 0.0};
    const Vec3   vPf{-kMu / h * std::sin(nu), kMu / h * (e + std::cos(nu)), 0.0};
    const auto   turn = [&](const Vec3& x) {
        const Vec3 a{std::cos(argp) * x.x - std::sin(argp) * x.y, std::sin(argp) * x.x + std::cos(argp) * x.y, x.z};
        const Vec3 b{a.x, std::cos(inc) * a.y - std::sin(inc) * a.z, std::sin(inc) * a.y + std::cos(inc) * a.z};
        return Vec3{std::cos(node) * b.x - std::sin(node) * b.y, std::sin(node) * b.x + std::cos(node) * b.y, b.z};
    };
    return {turn(rPf), turn(vPf)};
}

} // namespace

TEST(Stumpff, ReturnsBothFunctions) {
    // c2 and c3 come back as one struct, the members named as the functions.
    const Stumpff zero = stumpff(0.0);
    EXPECT_DOUBLE_EQ(zero.c2, 0.5);
    EXPECT_DOUBLE_EQ(zero.c3, 1.0 / 6.0);
    const Stumpff ellipse = stumpff(4.0); // sqrt(psi) = 2
    EXPECT_DOUBLE_EQ(ellipse.c2, (1.0 - std::cos(2.0)) / 4.0);
    EXPECT_DOUBLE_EQ(ellipse.c3, (2.0 - std::sin(2.0)) / 8.0);
    const Stumpff hyperbola = stumpff(-4.0);
    EXPECT_DOUBLE_EQ(hyperbola.c2, (std::cosh(2.0) - 1.0) / 4.0);
    EXPECT_DOUBLE_EQ(hyperbola.c3, (std::sinh(2.0) - 2.0) / 8.0);
    const Stumpff none{};
    EXPECT_EQ(none.c2, 0.0);
    EXPECT_EQ(none.c3, 0.0);
}

TEST(Stumpff, MatchesExtendedPrecision) {
    // c2 and c3 over elliptic (psi > 0) and hyperbolic (psi < 0) arguments, from 1e-12 to 400. Allowed:
    // a few rounding errors of the value, or of 1 / psi where the closed forms subtract (psi > 1; c2
    // and c3 pass through zero at psi = (2 pi k)^2).
    for (double sign : {1.0, -1.0}) {
        for (int k = -120; k <= 26; ++k) {
            const double psi    = sign * std::pow(10.0, 0.1 * k);
            const auto [c2, c3] = stumpff(psi);
            Real r0 = 0, r1 = 0, r2 = 0, r3 = 0;
            stumpffL(psi, r0, r1, r2, r3);
            const double floorC = std::abs(psi) > 1.0 ? 1.0 / std::abs(psi) : 0.0;
            const double eps    = std::numeric_limits<double>::epsilon();
            EXPECT_NEAR(c2, static_cast<double>(r2), 4.0 * eps * std::max(std::abs(c2), floorC)) << psi;
            EXPECT_NEAR(c3, static_cast<double>(r3), 4.0 * eps * std::max(std::abs(c3), floorC)) << psi;
        }
    }
}

TEST(PropagateKepler, ShortStepOnAnEccentricOrbit) {
    // Regression: this threw "did not converge". The Stumpff functions lost ten digits for small
    // arguments, so the iteration never reached its step tolerance.
    const StateVector s0  = elementsToState({70000e3, 0.9, 1.0, 0.3, 0.4, degToRad(73.0)});
    const StateVector s   = propagateKepler(s0, 5.5);
    const StateL      ref = referenceKepler(s0, 5.5, kMu);
    EXPECT_LT(distance(s.r, ref.r), 1e-13 * norm(s0.r));
    EXPECT_LT(distance(s.v, ref.v), 1e-13 * norm(s0.v));
}

TEST(PropagateKepler, MatchesExtendedPrecisionOnEveryConic) {
    // Perigee radius 6600-8000 km, any true anomaly, steps of 0.01 s to 1e5 s forwards and backwards:
    // circular to strongly hyperbolic, with the near-parabolic cases on both sides of e = 1. Before
    // the fix up to 1 % of the near-parabolic cases threw.
    //
    // Tolerance: 1e-12 of the size of the state. The largest errors measured are 1.5e-13 in position
    // (multi-revolution ellipses, where the period is known to a rounding error per revolution).
    std::mt19937_64                        gen(20261112);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    for (double e : {0.0, 0.5, 0.8, 0.9, 0.95, 0.99, 0.999, 0.999999, 1.0, 1.000001, 1.001, 1.5, 6.0, 50.0}) {
        for (int k = 0; k < 1500; ++k) {
            const double rp = 6600e3 + 1400e3 * uniform(gen);
            // An ellipse takes any true anomaly; an open orbit stays inside its asymptotes.
            const double nuMax = e < 1.0 ? constants::pi : 0.95 * std::acos(-1.0 / e);
            const double nu    = (2.0 * uniform(gen) - 1.0) * nuMax;
            const double node = constants::twoPi * uniform(gen), inc = std::acos(2.0 * uniform(gen) - 1.0);
            const double argp = constants::twoPi * uniform(gen);
            const double dt   = (uniform(gen) < 0.5 ? 1.0 : -1.0) * std::pow(10.0, -2.0 + 7.0 * uniform(gen));

            const StateVector s0 = conicState(rp, e, nu, node, inc, argp);
            StateVector       s;
            ASSERT_NO_THROW(s = propagateKepler(s0, dt)) << "e " << e << " nu " << nu << " dt " << dt;
            const StateL ref    = referenceKepler(s0, dt, kMu);
            const double rScale = norm(s0.r) + norm(s.r), vScale = norm(s0.v) + norm(s.v);
            ASSERT_LT(distance(s.r, ref.r), 1e-12 * rScale) << "e " << e << " nu " << nu << " dt " << dt;
            ASSERT_LT(distance(s.v, ref.v), 1e-12 * vScale) << "e " << e << " nu " << nu << " dt " << dt;
        }
    }
}

TEST(PropagateKepler, FlyByIsSymmetricAboutPeriapsis) {
    // A truth that does not use universal variables: an orbit is symmetric about its line of apsides.
    // Starting at true anomaly -nu, after twice the time from periapsis to +nu (the hyperbolic Kepler
    // equation, M = e sinh H - H) the state is the mirror image of the start.
    for (double e : {1.001, 1.2, 3.0}) {
        for (double nu : {0.3, 1.5, 2.2}) {
            if (nu >= 0.95 * std::acos(-1.0 / e)) {
                continue;
            }
            const double rp = 7000e3, a = rp / (e - 1.0); // |semi-major axis|
            const double H  = 2.0 * std::atanh(std::sqrt((e - 1.0) / (e + 1.0)) * std::tan(0.5 * nu));
            const double dt = 2.0 * (e * std::sinh(H) - H) / std::sqrt(kMu / (a * a * a));

            // In the orbit plane (no rotation): periapsis on +x, so the mirror flips y and v_x.
            const StateVector in  = conicState(rp, e, -nu, 0.0, 0.0, 0.0);
            const StateVector out = propagateKepler(in, dt);
            const double      tol = 1e-12 * (e < 1.01 ? 100.0 : 1.0); // H loses digits as e -> 1
            EXPECT_NEAR(out.r.x, in.r.x, tol * norm(in.r)) << e << " " << nu;
            EXPECT_NEAR(out.r.y, -in.r.y, tol * norm(in.r)) << e << " " << nu;
            EXPECT_NEAR(out.v.x, -in.v.x, tol * norm(in.v)) << e << " " << nu;
            EXPECT_NEAR(out.v.y, in.v.y, tol * norm(in.v)) << e << " " << nu;
        }
    }
}

TEST(PropagateKeplerStm, ShortStepOnAnEccentricOrbit) {
    // The same failure through the state-transition-matrix entry point. The matrix of a Hamiltonian
    // flow is symplectic: Phi^T J Phi = J.
    const StateVector     s0 = elementsToState({70000e3, 0.9, 1.0, 0.3, 0.4, degToRad(73.0)});
    const StateTransition st = propagateKeplerStm(s0, 5.5);
    const StateVector     s  = propagateKepler(s0, 5.5);
    EXPECT_EQ(st.state.r.x, s.r.x);
    EXPECT_EQ(st.state.v.z, s.v.z);
    // (Phi^T J Phi)_{ij} = sum_k Phi_{k,i} Phi_{k+3,j} - Phi_{k+3,i} Phi_{k,j}; position and velocity
    // blocks differ in units by seconds, so scale each entry by its own magnitude.
    for (std::size_t i = 0; i < 6; ++i) {
        for (std::size_t j = 0; j < 6; ++j) {
            double value = 0.0, scale = 0.0;
            for (std::size_t k = 0; k < 3; ++k) {
                const double a = st.stm(k, i) * st.stm(k + 3, j), b = st.stm(k + 3, i) * st.stm(k, j);
                value += a - b;
                scale += std::abs(a) + std::abs(b);
            }
            const double expected = (j == i + 3) ? 1.0 : (i == j + 3) ? -1.0 : 0.0;
            EXPECT_NEAR(value, expected, 1e-13 * std::max(scale, 1.0)) << i << " " << j;
        }
    }
}

TEST(PropagateKepler, RejectsInputItCannotPropagate) {
    const StateVector s0  = elementsToState({7000e3, 0.01, 1.0, 0.3, 0.4, 0.5});
    const double      nan = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW((void)propagateKepler(s0, nan), Error);
    EXPECT_THROW((void)propagateKepler(s0, std::numeric_limits<double>::infinity()), Error);
    EXPECT_THROW((void)propagateKepler(s0, 60.0, 0.0), Error);
    EXPECT_THROW((void)propagateKepler(s0, 60.0, -kMu), Error);
    EXPECT_THROW((void)propagateKepler({Vec3{}, s0.v}, 60.0), Error);
    EXPECT_THROW((void)propagateKepler(
                     {
                         Vec3{nan, 0.0, 0.0},
                         s0.v
    },
                     60.0),
                 Error);
    EXPECT_THROW((void)propagateKepler(
                     {
                         s0.r,
                         Vec3{0.0, nan, 0.0}
    },
                     60.0),
                 Error);
}
