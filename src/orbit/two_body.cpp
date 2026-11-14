/**********************************************************************
 * two_body.cpp
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

#include "astro/orbit/two_body.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <initializer_list>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

#include "astro/core/error.h"
#include "astro/orbit/gravity.h"
#include "astro/orbit/numerical.h"

namespace astro {

using constants::pi;

namespace {

constexpr double reciprocalFactorial(int n) noexcept {
    double f = 1.0;
    for (int i = 2; i <= n; ++i) {
        f *= i;
    }
    return 1.0 / f;
}

} // namespace

void stumpff(double psi, double& c2, double& c3) noexcept {
    if (std::abs(psi) < 1.0) {
        // The series c2 = sum (-psi)^k / (2k+2)!, c3 = sum (-psi)^k / (2k+3)!, to 1/20! and 1/21!. The
        // closed forms below subtract nearly equal numbers here: 1 - cos(sqrt(psi)) is good to one
        // rounding error of 1, a relative error of 2e-16 / psi in c2 (and 6e-16 / psi in c3).
        // clang-format off
        constexpr double a[10] = {reciprocalFactorial(2),  reciprocalFactorial(4),  reciprocalFactorial(6),
                                  reciprocalFactorial(8),  reciprocalFactorial(10), reciprocalFactorial(12),
                                  reciprocalFactorial(14), reciprocalFactorial(16), reciprocalFactorial(18),
                                  reciprocalFactorial(20)};
        constexpr double b[10] = {reciprocalFactorial(3),  reciprocalFactorial(5),  reciprocalFactorial(7),
                                  reciprocalFactorial(9),  reciprocalFactorial(11), reciprocalFactorial(13),
                                  reciprocalFactorial(15), reciprocalFactorial(17), reciprocalFactorial(19),
                                  reciprocalFactorial(21)};
        // clang-format on
        double s2 = a[9], s3 = b[9];
        for (int k = 8; k >= 0; --k) {
            s2 = a[k] - psi * s2;
            s3 = b[k] - psi * s3;
        }
        c2 = s2;
        c3 = s3;
    } else if (psi > 0.0) {
        const double sp = std::sqrt(psi);
        c2              = (1.0 - std::cos(sp)) / psi;
        c3              = (sp - std::sin(sp)) / (sp * psi);
    } else {
        const double sp = std::sqrt(-psi);
        c2              = (1.0 - std::cosh(sp)) / psi;
        c3              = (std::sinh(sp) - sp) / (sp * -psi);
    }
}

namespace {

/// The universal anomaly of a two-body propagation by dt (Vallado's KEPLER).
struct UniversalSolution {
    double chi{0.0};         ///< universal anomaly for dtReduced
    double dtReduced{0.0};   ///< dt, reduced modulo the period for an ellipse
    double revolutions{0.0}; ///< whole periods removed: dt = dtReduced + revolutions * period
};

UniversalSolution solveUniversal(const Vec3& r0, const Vec3& v0, double dt, double mu) {
    const double r0m = norm(r0);
    if (!(mu > 0.0) || !std::isfinite(mu) || !std::isfinite(dt) || !(r0m > 0.0) || !std::isfinite(r0m)
        || !std::isfinite(norm2(v0))) {
        throw Error("propagateKepler: the state and dt must be finite, the position non-zero and mu positive");
    }
    const double smu   = std::sqrt(mu);
    const double rdotv = dot(r0, v0);
    const double alpha = -norm2(v0) / mu + 2.0 / r0m; // 1/a

    UniversalSolution out;
    double            chi;
    if (alpha > 1e-12) { // ellipse: reduce dt modulo the period for robustness
        const double period = 2.0 * pi / (smu * alpha * std::sqrt(alpha));
        if (std::abs(dt) > period) {
            const double reduced = std::fmod(dt, period);
            out.revolutions      = std::round((dt - reduced) / period);
            dt                   = reduced;
        }
        chi = smu * dt * alpha;
    } else if (alpha < -1e-12) { // hyperbola
        const double a  = 1.0 / alpha;
        const double sg = dt > 0.0 ? 1.0 : -1.0;
        chi             = sg * std::sqrt(-a)
            * std::log((-2.0 * mu * alpha * dt) / (rdotv + sg * std::sqrt(-mu * a) * (1.0 - r0m * alpha)));
    } else { // parabola
        const Vec3   h = cross(r0, v0);
        const double p = norm2(h) / mu;
        const double s = 0.5 * std::atan2(1.0, 3.0 * std::sqrt(mu / (p * p * p)) * dt);
        const double w = std::atan(std::cbrt(std::tan(s)));
        chi            = std::sqrt(p) * 2.0 / std::tan(2.0 * w);
    }

    // The hyperbolic starter is meant for long times: for short ones its logarithm is negative, or not
    // a number. Then start from zero, the known end of the bracket below.
    if (!std::isfinite(chi) || std::signbit(chi) != std::signbit(dt)) {
        chi = 0.0;
    }

    // Laguerre-Conway iteration on the universal Kepler equation (robust for all conic types):
    //   F(chi) = (r0.v0/sqrt(mu)) chi^2 c2 + (1 - alpha r0) chi^3 c3 + r0 chi - sqrt(mu) dt
    //
    // F increases with chi (its derivative is the radius) and F(0) = -sqrt(mu) dt, so the root lies
    // between 0 and any chi where F has the sign of dt. The iterates keep such a bracket [lo, hi], and
    // a step that would leave it is replaced by a point inside: the iteration cannot fail to end.
    //
    // It ends when the step is below rounding, or when F itself is. The second matters near a close
    // periapsis, where the terms of F cancel: their rounding error divided by the small radius is then
    // a larger part of chi than any fixed step tolerance, and no iterate can do better.
    constexpr double kStepTol  = 1e-13;
    constexpr double kNoise    = 4.0 * std::numeric_limits<double>::epsilon();
    constexpr double kInfinity = std::numeric_limits<double>::infinity();
    constexpr double nL        = 5.0;
    const bool       forward   = dt > 0.0;
    const double     sig0      = rdotv / smu;
    const double     k3        = 1.0 - alpha * r0m;
    const double     t0        = smu * dt;
    double           lo        = forward ? 0.0 : -kInfinity;
    double           hi        = forward ? kInfinity : 0.0;
    bool             converged = false;
    for (int it = 0; it < 200 && !converged; ++it) {
        const double psi = chi * chi * alpha;
        double       c2 = 0.0, c3 = 0.0;
        stumpff(psi, c2, c3);
        const double t2 = sig0 * chi * chi * c2, t3 = k3 * chi * chi * chi * c3, t1 = r0m * chi;
        const double F = t2 + t3 + t1 - t0;
        if (F == 0.0) { // also the whole number of periods: dt reduced to 0, chi = 0
            converged = true;
            break;
        }
        // Past the root if F is positive, or if it has overflowed (then chi is far beyond the root).
        const bool farSide = std::isfinite(F) ? F > 0.0 : forward;
        if (farSide) {
            hi = std::min(hi, chi);
        } else {
            lo = std::max(lo, chi);
        }

        double next = std::numeric_limits<double>::quiet_NaN();
        if (std::isfinite(F)) {
            // Rounding error of F: a few units in the last place of each term, and of c2 and c3
            // where they pass through zero (psi > 1, where they are below 2 / psi and 1 / psi).
            const double floorC = std::abs(psi) > 1.0 ? 1.0 / std::abs(psi) : 1.0;
            const double noise  = kNoise
                               * (std::abs(t2) + std::abs(t3) + std::abs(t1) + std::abs(t0)
                                  + floorC * (std::abs(sig0) + std::abs(k3 * chi)) * chi * chi);
            if (std::abs(F) <= noise) {
                converged = true;
                break;
            }
            const double dF   = sig0 * chi * (1.0 - psi * c3) + k3 * chi * chi * c2 + r0m; // = r
            const double ddF  = sig0 * (1.0 - psi * c2) + k3 * chi * (1.0 - psi * c3);
            const double disc = std::sqrt(std::abs((nL - 1.0) * (nL - 1.0) * dF * dF - nL * (nL - 1.0) * F * ddF));
            // Laguerre's step; Newton's if the discriminant has overflowed or the denominator is zero.
            const double denom = dF + (dF >= 0.0 ? disc : -disc);
            const double step  = std::isfinite(disc) && denom != 0.0 ? -nL * F / denom : -F / dF;
            if (std::abs(step) <= kStepTol * std::abs(chi)) {
                chi       += step;
                converged  = true;
                break;
            }
            next = chi + step;
        }
        if (!(next > lo && next < hi)) {
            // Only rounding or overflow brings the step here. With both ends known take the middle;
            // otherwise the root is still further out than chi.
            next      = std::isfinite(lo) && std::isfinite(hi) ? 0.5 * (lo + hi) : 2.0 * chi;
            converged = next == lo || next == hi; // the bracket has closed to neighbouring numbers
        }
        chi = next;
    }
    if (!converged) {
        throw Error("propagateKepler: universal-variable iteration did not converge");
    }
    out.chi       = chi;
    out.dtReduced = dt;
    return out;
}

/// The propagated state from the universal anomaly (Lagrange coefficients f, g, fdot, gdot).
StateVector universalState(const Vec3& r0, const Vec3& v0, const UniversalSolution& u, double mu) {
    const double r0m   = norm(r0);
    const double smu   = std::sqrt(mu);
    const double rdotv = dot(r0, v0);
    const double alpha = -norm2(v0) / mu + 2.0 / r0m;
    const double chi   = u.chi;
    const double psi   = chi * chi * alpha;
    double       c2 = 0.0, c3 = 0.0;
    stumpff(psi, c2, c3);
    const double r    = chi * chi * c2 + rdotv / smu * chi * (1.0 - psi * c3) + r0m * (1.0 - psi * c2);
    const double f    = 1.0 - chi * chi / r0m * c2;
    const double g    = u.dtReduced - chi * chi * chi / smu * c3;
    const double gdot = 1.0 - chi * chi / r * c2;
    const double fdot = smu / (r * r0m) * chi * (psi * c3 - 1.0);
    return {f * r0 + g * v0, fdot * r0 + gdot * v0};
}

/// Universal functions U0..U3 at (chi, alpha) and their partial derivatives with respect to
/// alpha at fixed chi, dU_k/dalpha = (k U_{k+2} - chi U_{k+1}) / 2.
struct UniversalFunctions {
    double U[4];
    double dAlpha[4];
};

UniversalFunctions universalFunctions(double chi, double alpha) {
    // Stumpff functions c_k(z), z = alpha chi^2, with c_k + z c_{k+2} = 1/k!.
    const double z = alpha * chi * chi;
    double       c0, c1, c2, c3, c4 = 0.0, c5 = 0.0;
    const bool   series = std::abs(z) < 1.0;
    if (series) { // c4 and c5 by their series, the lower ones by the recurrence: no cancellation
        double t4 = 1.0 / 24.0, t5 = 1.0 / 120.0;
        for (int j = 0; j < 16; ++j) {
            c4 += t4;
            c5 += t5;
            t4 *= -z / ((5.0 + 2 * j) * (6.0 + 2 * j));
            t5 *= -z / ((6.0 + 2 * j) * (7.0 + 2 * j));
        }
        c3 = 1.0 / 6.0 - z * c5;
        c2 = 0.5 - z * c4;
        c1 = 1.0 - z * c3;
        c0 = 1.0 - z * c2;
    } else if (z > 0.0) {
        const double sq = std::sqrt(z);
        c0              = std::cos(sq);
        c1              = std::sin(sq) / sq;
        c2              = (1.0 - c0) / z;
        c3              = (1.0 - c1) / z;
    } else {
        const double sq = std::sqrt(-z);
        c0              = std::cosh(sq);
        c1              = std::sinh(sq) / sq;
        c2              = (1.0 - c0) / z;
        c3              = (1.0 - c1) / z;
    }
    UniversalFunctions u;
    const double       chi2 = chi * chi;
    u.U[0]                  = c0;
    u.U[1]                  = chi * c1;
    u.U[2]                  = chi2 * c2;
    u.U[3]                  = chi2 * chi * c3;
    u.dAlpha[0]             = -0.5 * chi * u.U[1];
    u.dAlpha[1]             = 0.5 * (u.U[3] - chi * u.U[2]);
    if (series) {
        const double U4 = chi2 * chi2 * c4, U5 = chi2 * chi2 * chi * c5;
        u.dAlpha[2] = 0.5 * (2.0 * U4 - chi * u.U[3]);
        u.dAlpha[3] = 0.5 * (3.0 * U5 - chi * U4);
    } else { // the same, with U4 and U5 eliminated (their leading terms cancel for large chi)
        u.dAlpha[2] = (chi * u.U[1] - 2.0 * u.U[2]) / (2.0 * alpha);
        u.dAlpha[3] = (chi * u.U[2] - 3.0 * u.U[3]) / (2.0 * alpha);
    }
    return u;
}

using Grad = std::array<double, 6>; // gradient with respect to (r0, v0)

Grad combine(std::initializer_list<std::pair<double, const Grad*>> terms) noexcept {
    Grad g{};
    for (const auto& [c, v] : terms) {
        for (std::size_t j = 0; j < 6; ++j) {
            g[j] += c * (*v)[j];
        }
    }
    return g;
}

} // namespace

StateVector propagateKepler(const StateVector& s0, double dt, double mu) {
    if (dt == 0.0) {
        return s0;
    }
    return universalState(s0.r, s0.v, solveUniversal(s0.r, s0.v, dt, mu), mu);
}

StateTransition propagateKeplerStm(const StateVector& s0, double dt, double mu) {
    if (dt == 0.0) {
        return {s0, identity6()};
    }
    const Vec3&             r0  = s0.r;
    const Vec3&             v0  = s0.v;
    const UniversalSolution sol = solveUniversal(r0, v0, dt, mu);
    StateTransition         out{universalState(r0, v0, sol, mu), {}};

    // Phi by differentiating r = F r0 + G v0 and v = Fdot r0 + Gdot v0 (Lagrange coefficients in
    // universal variables) with respect to (r0, v0). The coefficients depend on r0 = |r0|,
    // sigma0 = r0.v0 / sqrt(mu), alpha = 1/a and the universal anomaly chi, which is defined
    // implicitly by Kepler's equation sqrt(mu) dt = r0 U1 + sigma0 U2 + U3. chi must be that of the
    // full dt: the period depends on alpha, so the whole revolutions removed by solveUniversal
    // carry the along-track drift (chi grows by 2 pi / sqrt(alpha) per revolution).
    const double r0m   = norm(r0);
    const double smu   = std::sqrt(mu);
    const double sig0  = dot(r0, v0) / smu;
    const double alpha = -norm2(v0) / mu + 2.0 / r0m;
    const double chi   = sol.chi + (sol.revolutions != 0.0 ? sol.revolutions * 2.0 * pi / std::sqrt(alpha) : 0.0);
    const UniversalFunctions uf = universalFunctions(chi, alpha);
    const double             U0 = uf.U[0], U1 = uf.U[1], U2 = uf.U[2];

    const double r03 = r0m * r0m * r0m;
    const Grad   gr0{r0.x / r0m, r0.y / r0m, r0.z / r0m, 0.0, 0.0, 0.0};
    const Grad   gsig{v0.x / smu, v0.y / smu, v0.z / smu, r0.x / smu, r0.y / smu, r0.z / smu};
    const Grad   galpha{-2.0 * r0.x / r03,
                      -2.0 * r0.y / r03,
                      -2.0 * r0.z / r03,
                      -2.0 * v0.x / mu,
                      -2.0 * v0.y / mu,
                      -2.0 * v0.z / mu};

    // Kepler's equation f = r0 U1 + sigma0 U2 + U3 - sqrt(mu) dt = 0, with df/dchi = r.
    const double r      = r0m * U0 + sig0 * U1 + U2;
    const double fAlpha = r0m * uf.dAlpha[1] + sig0 * uf.dAlpha[2] + uf.dAlpha[3];
    const Grad   gchi   = combine({
        {-U1 / r,     &gr0   },
        {-U2 / r,     &gsig  },
        {-fAlpha / r, &galpha}
    });
    const Grad   gU0    = combine({
        {-alpha * U1,  &gchi  },
        {uf.dAlpha[0], &galpha}
    });
    const Grad   gU1    = combine({
        {U0,           &gchi  },
        {uf.dAlpha[1], &galpha}
    });
    const Grad   gU2    = combine({
        {U1,           &gchi  },
        {uf.dAlpha[2], &galpha}
    });
    const Grad   gr     = combine({
        {U0,   &gr0 },
        {r0m,  &gU0 },
        {U1,   &gsig},
        {sig0, &gU1 },
        {1.0,  &gU2 }
    });

    const double F     = 1.0 - U2 / r0m;
    const double G     = (r0m * U1 + sig0 * U2) / smu;
    const double Fdot  = -smu * U1 / (r * r0m);
    const double Gdot  = 1.0 - U2 / r;
    const Grad   gF    = combine({
        {-1.0 / r0m,       &gU2},
        {U2 / (r0m * r0m), &gr0}
    });
    const Grad   gG    = combine({
        {U1 / smu,   &gr0 },
        {r0m / smu,  &gU1 },
        {U2 / smu,   &gsig},
        {sig0 / smu, &gU2 }
    });
    const Grad   gFdot = combine({
        {-smu / (r * r0m),           &gU1},
        {smu * U1 / (r * r * r0m),   &gr },
        {smu * U1 / (r * r0m * r0m), &gr0}
    });
    const Grad   gGdot = combine({
        {-1.0 / r,     &gU2},
        {U2 / (r * r), &gr }
    });

    const double rc[3] = {r0.x, r0.y, r0.z}, vc[3] = {v0.x, v0.y, v0.z};
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 6; ++j) {
            out.stm[i][j]     = rc[i] * gF[j] + vc[i] * gG[j];
            out.stm[i + 3][j] = rc[i] * gFdot[j] + vc[i] * gGdot[j];
        }
        out.stm[i][i]         += F;
        out.stm[i][i + 3]     += G;
        out.stm[i + 3][i]     += Fdot;
        out.stm[i + 3][i + 3] += Gdot;
    }
    return out;
}

StateCovariance propagateKeplerWithCovariance(const StateVector& s0, const Mat6& P0, double dt, double mu) {
    const StateTransition t = propagateKeplerStm(s0, dt, mu);
    return {t.state, congruence(t.stm, P0)};
}

Propagation12 propagateKepler(const StateVector12& s0, double dt, double mu, HigherOrderMode mode) {
    Propagation12 out;
    out.accelerationResidual = s0.a - pointMassAcceleration(s0.r, mu);
    out.jerkResidual         = s0.j - pointMassJerk(s0.r, s0.v, mu); // point-mass jerk does not depend on a
    const bool unmodeled =
        mode == HigherOrderMode::UnmodeledForce && (out.accelerationResidual != Vec3{} || out.jerkResidual != Vec3{});
    if (unmodeled) {
        const NumericalPropagator prop(Time(), s0, {std::make_shared<GravityField>(CentralBody::pointMass(mu))}, mode);
        out.state = prop.propagate12(dt);
    } else {
        const StateVector s = propagateKepler(s0.rv(), dt, mu);
        out.state           = {s.r, s.v, pointMassAcceleration(s.r, mu), pointMassJerk(s.r, s.v, mu)};
    }
    return out;
}

} // namespace astro
