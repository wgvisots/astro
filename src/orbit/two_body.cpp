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
#include <cmath>
#include <memory>
#include <vector>

#include "astro/core/error.h"
#include "astro/orbit/gravity.h"
#include "astro/orbit/numerical.h"

namespace astro {

using constants::pi;

void stumpff(double psi, double& c2, double& c3) noexcept {
    if (psi > 1e-6) {
        const double sp = std::sqrt(psi);
        c2              = (1.0 - std::cos(sp)) / psi;
        c3              = (sp - std::sin(sp)) / (sp * psi);
    } else if (psi < -1e-6) {
        const double sp = std::sqrt(-psi);
        c2              = (1.0 - std::cosh(sp)) / psi;
        c3              = (std::sinh(sp) - sp) / (sp * -psi);
    } else {
        c2 = 0.5 - psi / 24.0 + psi * psi / 720.0;
        c3 = 1.0 / 6.0 - psi / 120.0 + psi * psi / 5040.0;
    }
}

StateVector propagateKepler(const StateVector& s0, double dt, double mu) {
    if (dt == 0.0) {
        return s0;
    }
    const Vec3&  r0    = s0.r;
    const Vec3&  v0    = s0.v;
    const double r0m   = norm(r0);
    const double smu   = std::sqrt(mu);
    const double rdotv = dot(r0, v0);
    const double alpha = -norm2(v0) / mu + 2.0 / r0m; // 1/a

    double chi;
    if (alpha > 1e-12) { // ellipse: reduce dt modulo the period for robustness
        const double period = 2.0 * pi / (smu * alpha * std::sqrt(alpha));
        if (std::abs(dt) > period) {
            dt = std::fmod(dt, period);
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

    if (!std::isfinite(chi)) {
        chi = smu * dt / r0m; // safe small-time starter
    }

    // Laguerre-Conway iteration on the universal Kepler equation (robust for all conic types):
    //   F(chi) = (r0.v0/sqrt(mu)) chi^2 c2 + (1 - alpha r0) chi^3 c3 + r0 chi - sqrt(mu) dt
    const double     sig0 = rdotv / smu;
    double           c2 = 0.0, c3 = 0.0, psi = 0.0, r = r0m;
    bool             converged = false;
    constexpr double nL        = 5.0;
    for (int it = 0; it < 200; ++it) {
        psi = chi * chi * alpha;
        stumpff(psi, c2, c3);
        const double F      = sig0 * chi * chi * c2 + (1.0 - alpha * r0m) * chi * chi * chi * c3 + r0m * chi - smu * dt;
        const double dF     = sig0 * chi * (1.0 - psi * c3) + (1.0 - alpha * r0m) * chi * chi * c2 + r0m; // = r
        const double ddF    = sig0 * (1.0 - psi * c2) + (1.0 - alpha * r0m) * chi * (1.0 - psi * c3);
        const double disc   = std::sqrt(std::abs((nL - 1.0) * (nL - 1.0) * dF * dF - nL * (nL - 1.0) * F * ddF));
        const double denom  = dF + (dF >= 0.0 ? disc : -disc);
        const double dchi   = denom != 0.0 ? -nL * F / denom : -F / dF;
        chi                += dchi;
        if (std::abs(dchi) < 1e-13 * (1.0 + std::abs(chi))) {
            converged = true;
            break;
        }
    }
    if (!converged) {
        throw Error("propagateKepler: universal-variable iteration did not converge");
    }
    psi = chi * chi * alpha;
    stumpff(psi, c2, c3);
    r                 = chi * chi * c2 + rdotv / smu * chi * (1.0 - psi * c3) + r0m * (1.0 - psi * c2);
    const double f    = 1.0 - chi * chi / r0m * c2;
    const double g    = dt - chi * chi * chi / smu * c3;
    const double gdot = 1.0 - chi * chi / r * c2;
    const double fdot = smu / (r * r0m) * chi * (psi * c3 - 1.0);
    return {f * r0 + g * v0, fdot * r0 + gdot * v0};
}

// ------------------------------------------------------------------ Izzo (2015)

namespace {

struct Izzo {
    double lambda;

    static double hypergeometricF(double z, double tol) {
        double Sj = 1.0, Cj = 1.0, err = 1.0;
        int    j = 0;
        while (err > tol && j < 1000) {
            const double Cj1  = Cj * (3.0 + j) * (1.0 + j) / (2.5 + j) * z / (j + 1);
            Sj               += Cj1;
            err               = std::abs(Cj1);
            Cj                = Cj1;
            ++j;
        }
        return Sj;
    }

    double x2tofLagrange(double x, int N) const {
        const double a = 1.0 / (1.0 - x * x);
        if (a > 0.0) {
            const double alfa = 2.0 * std::acos(x);
            double       beta = 2.0 * std::asin(std::sqrt(lambda * lambda / a));
            if (lambda < 0.0) {
                beta = -beta;
            }
            return a * std::sqrt(a) * ((alfa - std::sin(alfa)) - (beta - std::sin(beta)) + 2.0 * pi * N) / 2.0;
        }
        const double alfa = 2.0 * std::acosh(x);
        double       beta = 2.0 * std::asinh(std::sqrt(-lambda * lambda / a));
        if (lambda < 0.0) {
            beta = -beta;
        }
        return -a * std::sqrt(-a) * ((beta - std::sinh(beta)) - (alfa - std::sinh(alfa))) / 2.0;
    }

    double x2tof(double x, int N) const {
        constexpr double battin = 0.01, lagrange = 0.2;
        const double     dist = std::abs(x - 1.0);
        if (dist < lagrange && dist > battin) {
            return x2tofLagrange(x, N);
        }
        const double K   = lambda * lambda;
        const double E   = x * x - 1.0;
        const double rho = std::abs(E);
        const double z   = std::sqrt(1.0 + K * E);
        if (dist < battin) {
            const double eta = z - lambda * x;
            const double S1  = 0.5 * (1.0 - lambda - x * eta);
            const double Q   = 4.0 / 3.0 * hypergeometricF(S1, 1e-11);
            return (eta * eta * eta * Q + 4.0 * lambda * eta) / 2.0 + N * pi / std::pow(rho, 1.5);
        }
        const double y = std::sqrt(rho);
        const double g = x * z - lambda * E;
        double       d;
        if (E < 0.0) {
            d = N * pi + std::acos(std::clamp(g, -1.0, 1.0));
        } else {
            const double f = y * (z - lambda * x);
            d              = std::log(f + g);
        }
        return (x - lambda * z - d / y) / E;
    }

    void dTdx(double x, double T, double& DT, double& DDT, double& DDDT) const {
        const double l2 = lambda * lambda, l3 = l2 * lambda;
        const double umx2 = 1.0 - x * x;
        const double y    = std::sqrt(1.0 - l2 * umx2);
        const double y2 = y * y, y3 = y2 * y;
        DT   = 1.0 / umx2 * (3.0 * T * x - 2.0 + 2.0 * l3 * x / y);
        DDT  = 1.0 / umx2 * (3.0 * T + 5.0 * x * DT + 2.0 * (1.0 - l2) * l3 / y3);
        DDDT = 1.0 / umx2 * (7.0 * x * DDT + 8.0 * DT - 6.0 * (1.0 - l2) * l2 * l3 * x / (y3 * y2));
    }

    double householder(double T, double x0, int N, double eps, int maxIter) const {
        for (int it = 0; it < maxIter; ++it) {
            const double tof = x2tof(x0, N);
            double       DT, DDT, DDDT;
            dTdx(x0, tof, DT, DDT, DDDT);
            const double delta = tof - T;
            const double DT2   = DT * DT;
            const double xnew =
                x0 - delta * (DT2 - delta * DDT / 2.0) / (DT * (DT2 - delta * DDT) + DDDT * delta * delta / 6.0);
            const double err = std::abs(x0 - xnew);
            x0               = xnew;
            if (err <= eps) {
                break;
            }
        }
        return x0;
    }
};

} // namespace

std::vector<LambertSolution>
lambertIzzo(const Vec3& r1, const Vec3& r2, double tof, double mu, bool prograde, int maxRevolutions) {
    if (tof <= 0.0) {
        throw Error("lambertIzzo: time of flight must be positive");
    }
    if (mu <= 0.0) {
        throw Error("lambertIzzo: gravitational parameter must be positive");
    }
    const double R1 = norm(r1), R2 = norm(r2);
    const double c   = norm(r2 - r1);
    const double s   = 0.5 * (c + R1 + R2);
    const Vec3   ir1 = r1 / R1, ir2 = r2 / R2;
    Vec3         ih = cross(ir1, ir2);
    if (norm(ih) == 0.0) {
        throw Error("lambertIzzo: r1 and r2 are collinear; transfer plane undefined");
    }
    ih                   = normalized(ih);
    const double lambda2 = 1.0 - c / s;
    double       lambda  = std::sqrt(std::max(lambda2, 0.0));
    Vec3         it1, it2;
    if (ih.z < 0.0) { // transfer angle > 180 deg as seen from +z
        lambda = -lambda;
        it1    = cross(ir1, ih);
        it2    = cross(ir2, ih);
    } else {
        it1 = cross(ih, ir1);
        it2 = cross(ih, ir2);
    }
    it1 = normalized(it1);
    it2 = normalized(it2);
    if (!prograde) {
        lambda = -lambda;
        it1    = -it1;
        it2    = -it2;
    }
    const double lambda3 = lambda * lambda2;
    const double T       = std::sqrt(2.0 * mu / (s * s * s)) * tof;
    const Izzo   iz{lambda};

    // Maximum number of revolutions for which a solution exists.
    int          Nmax = static_cast<int>(std::floor(T / pi));
    const double T00  = std::acos(lambda) + lambda * std::sqrt(1.0 - lambda2);
    const double T0   = T00 + Nmax * pi;
    const double T1   = 2.0 / 3.0 * (1.0 - lambda3);
    if (Nmax > 0 && T < T0) {
        // Halley iterations for the minimum-time point of the Nmax family.
        double xOld = 0.0, xNew = 0.0, Tmin = T0;
        for (int it = 0; it < 13; ++it) {
            double DT, DDT, DDDT;
            iz.dTdx(xOld, Tmin, DT, DDT, DDDT);
            if (DT != 0.0) {
                xNew = xOld - DT * DDT / (DDT * DDT - DT * DDDT / 2.0);
            }
            const double err = std::abs(xOld - xNew);
            if (err < 1e-13) {
                break;
            }
            Tmin = iz.x2tof(xNew, Nmax);
            xOld = xNew;
        }
        if (Tmin > T) {
            Nmax -= 1;
        }
    }
    Nmax = std::min(Nmax, maxRevolutions);

    std::vector<double> xs;
    std::vector<int>    revs;
    std::vector<bool>   right;
    // Zero-revolution solution.
    double x0;
    if (T >= T00) {
        x0 = -(T - T00) / (T - T00 + 4.0);
    } else if (T <= T1) {
        x0 = T1 * (T1 - T) / (2.0 / 5.0 * (1.0 - lambda2 * lambda3) * T) + 1.0;
    } else {
        x0 = std::pow(T / T00, 0.69314718055994529 / std::log(T1 / T00)) - 1.0;
    }
    xs.push_back(iz.householder(T, x0, 0, 1e-5, 15));
    revs.push_back(0);
    right.push_back(false);
    for (int i = 1; i <= Nmax; ++i) {
        double tmp = std::pow((i * pi + pi) / (8.0 * T), 2.0 / 3.0);
        xs.push_back(iz.householder(T, (tmp - 1.0) / (tmp + 1.0), i, 1e-8, 15));
        revs.push_back(i);
        right.push_back(false);
        tmp = std::pow((8.0 * T) / (i * pi), 2.0 / 3.0);
        xs.push_back(iz.householder(T, (tmp - 1.0) / (tmp + 1.0), i, 1e-8, 15));
        revs.push_back(i);
        right.push_back(true);
    }

    const double                 gamma = std::sqrt(mu * s / 2.0);
    const double                 rho   = (R1 - R2) / c;
    const double                 sigma = std::sqrt(1.0 - rho * rho);
    std::vector<LambertSolution> out;
    for (std::size_t k = 0; k < xs.size(); ++k) {
        const double x   = xs[k];
        const double y   = std::sqrt(1.0 - lambda2 + lambda2 * x * x);
        const double vr1 = gamma * ((lambda * y - x) - rho * (lambda * y + x)) / R1;
        const double vr2 = -gamma * ((lambda * y - x) + rho * (lambda * y + x)) / R2;
        const double vt  = gamma * sigma * (y + lambda * x);
        out.push_back({vr1 * ir1 + (vt / R1) * it1, vr2 * ir2 + (vt / R2) * it2, revs[k], right[k]});
    }
    return out;
}

LambertSolution lambert(const Vec3& r1, const Vec3& r2, double tof, double mu, bool prograde) {
    return lambertIzzo(r1, r2, tof, mu, prograde, 0).front();
}

Propagation12 propagateKepler(const StateVector12& s0, double dt, double mu, HigherOrderMode mode) {
    Propagation12 out;
    out.accelerationResidual = s0.a - pointMassAcceleration(s0.r, mu);
    out.jerkResidual         = s0.j - pointMassJerk(s0.r, s0.v, mu); // point-mass jerk does not depend on a
    const bool unmodeled =
        mode == HigherOrderMode::UnmodeledForce && (out.accelerationResidual != Vec3{} || out.jerkResidual != Vec3{});
    if (unmodeled) {
        const NumericalPropagator prop(Time(), s0, {std::make_shared<GravityField>(std::vector<double>{}, mu)}, mode);
        out.state = prop.propagate12(dt);
    } else {
        const StateVector s = propagateKepler(s0.rv(), dt, mu);
        out.state           = {s.r, s.v, pointMassAcceleration(s.r, mu), pointMassJerk(s.r, s.v, mu)};
    }
    return out;
}

} // namespace astro
