/**********************************************************************
 * gravity_detail.h
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

#pragma once

#include <cmath>
#include <cstddef>

#include "astro/core/constants.h"
#include "astro/core/vec.h"
#include "astro/orbit/central_body.h"

// Gravity internals: the zonal recursion with its partial derivatives, and the jerks and gradients
// built on it. Header-only, so that the library and the tests share one definition without
// exporting it. Internal: not installed.
namespace astro::detail {

/// Recursion factors (2n - 1) / n and (n - 1) / n for the Legendre polynomials, tabulated to avoid
/// a division per degree.
struct LegendreFactors {
    static constexpr std::size_t size = 64;
    double                       a[size]{};
    double                       b[size]{};
    constexpr LegendreFactors() {
        for (std::size_t n = 2; n < size; ++n) {
            a[n] = (2.0 * static_cast<double>(n) - 1.0) / static_cast<double>(n);
            b[n] = (static_cast<double>(n) - 1.0) / static_cast<double>(n);
        }
    }
};
inline constexpr LegendreFactors kLegendre;

/// Acceleration and its gradient G = d(acceleration)/dr in factored form
/// G = U rhat^T + W w^T + c I, so that the jerk along a velocity v is G v.
struct GravityPartials {
    Vec3   acceleration;
    Vec3   U, W, rhat, w;
    double c{0.0};

    [[nodiscard]] Vec3 jerk(const Vec3& v) const noexcept { return U * dot(rhat, v) + W * dot(w, v) + c * v; }
    [[nodiscard]] Mat3 gradient() const noexcept { return outer(U, rhat) + outer(W, w) + Mat3::identity() * c; }
};

/// Point mass (if pointMass) plus zonal acceleration and, if WithPartials, the gradient in the
/// factored form of GravityPartials. Allocation-free: P_n(s), P_n'(s) and P_n''(s) are carried in
/// a rolling window of two degrees.
///
///   a = mu / r^2 [ cr r_hat + cz z_hat ],   s = z / r,   k_n = Jn (re / r)^n,
///   cr = -pm + sum k_n ((n + 1) P_n + s P_n'),   cz = -sum k_n P_n',   pm = 1 with the point mass.
///
/// Differentiating along a velocity v, with r' = r_hat . v and s' = w . v, w = (z_hat - s r_hat) / r:
///   cr' = Ar r' + As s',   Ar = -sum n k_n ((n + 1) P_n + s P_n') / r,  As = sum k_n ((n + 2) P_n' + s P_n'')
///   cz' = Br r' + Bs s',   Br =  sum n k_n P_n' / r,                    Bs = -sum k_n P_n''
/// which gives jerk = U (r_hat . v) + W (w . v) + c v with
///   U = -(2 g / r)(cr r_hat + cz z_hat) + g (Ar r_hat + Br z_hat) - (g cr / r) r_hat,
///   W = g (As r_hat + Bs z_hat),   c = g cr / r,   g = mu / r^2.
template <bool WithPartials>
void gravity(const Vec3&      r,
             const double*    J,
             std::size_t      size,
             double           mu,
             double           re,
             bool             pointMass,
             Vec3*            acc,
             GravityPartials* partials) noexcept {
    const double invR = 1.0 / std::sqrt(dot(r, r));
    const Vec3   rhat = r * invR;
    const double s    = rhat.z;
    const double q    = re * invR;

    double P2 = 1.0, P1 = s; // degrees n-2 and n-1 of P, P', P''
    double dP2 = 0.0, dP1 = 1.0;
    double ddP2 = 0.0, ddP1 = 0.0;
    double qn = q;
    double cr = pointMass ? -1.0 : 0.0, cz = 0.0;
    double Ar = 0.0, As = 0.0, Br = 0.0, Bs = 0.0;
    for (std::size_t n = 2; n < size; ++n) {
        const double dn  = static_cast<double>(n);
        const double an  = n < LegendreFactors::size ? kLegendre.a[n] : (2.0 * dn - 1.0) / dn;
        const double bn  = n < LegendreFactors::size ? kLegendre.b[n] : (dn - 1.0) / dn;
        const double P   = an * s * P1 - bn * P2;
        const double dP  = dP2 + (2.0 * dn - 1.0) * P1;
        double       ddP = 0.0;
        if constexpr (WithPartials) {
            ddP = ddP2 + (2.0 * dn - 1.0) * dP1;
        }
        qn *= q;
        if (J[n] != 0.0) {
            const double k   = J[n] * qn;
            const double fr  = (dn + 1.0) * P + s * dP;
            cr              += k * fr;
            cz              -= k * dP;
            if constexpr (WithPartials) {
                Ar -= dn * k * fr;
                As += k * ((dn + 2.0) * dP + s * ddP);
                Br += dn * k * dP;
                Bs -= k * ddP;
            }
        }
        P2  = P1;
        P1  = P;
        dP2 = dP1;
        dP1 = dP;
        if constexpr (WithPartials) {
            ddP2 = ddP1;
            ddP1 = ddP;
        }
    }
    const double g     = mu * invR * invR;
    const Vec3   zhat  = {0.0, 0.0, 1.0};
    const Vec3   accel = g * (cr * rhat + cz * zhat);
    if (acc) {
        *acc = accel;
    }
    if constexpr (WithPartials) {
        Ar                     *= invR;
        Br                     *= invR;
        partials->acceleration  = accel;
        partials->rhat          = rhat;
        partials->w             = (zhat - s * rhat) * invR;
        partials->c             = g * cr * invR;
        partials->U             = -2.0 * invR * accel + g * (Ar * rhat + Br * zhat) - partials->c * rhat;
        partials->W             = g * (As * rhat + Bs * zhat);
    }
}

/// Point-mass jerk (time derivative of pointMassAcceleration along the motion):
/// -mu / |r|^3 (v - 3 (r.v) r / |r|^2).
[[nodiscard]] inline Vec3 pointMassJerk(const Vec3& r, const Vec3& v, double mu = constants::gmEarth) noexcept {
    const double r2 = norm2(r), rm = std::sqrt(r2);
    return (v - r * (3.0 * dot(r, v) / r2)) * (-mu / (r2 * rm));
}

/// Gradient d(pointMassAcceleration)/dr = mu / |r|^3 (3 r_hat r_hat^T - I).
[[nodiscard]] inline Mat3 pointMassGradient(const Vec3& r, double mu = constants::gmEarth) noexcept {
    const double r2 = dot(r, r), rm = std::sqrt(r2);
    return (outer(r, r) * (3.0 / r2) - Mat3::identity()) * (mu / (r2 * rm));
}

/// Gradient d(zonalAcceleration)/dr (symmetric, exact).
[[nodiscard]] inline Mat3 zonalGradient(const Vec3& r, const CentralBody& body = CentralBody::earth()) noexcept {
    GravityPartials p;
    gravity<true>(r, body.J.data(), body.J.size(), body.mu, body.re, false, nullptr, &p);
    return p.gradient();
}

/// Time derivative of zonalAcceleration along the motion (r moving with velocity v). Exact.
[[nodiscard]] inline Vec3
zonalJerk(const Vec3& r, const Vec3& v, const CentralBody& body = CentralBody::earth()) noexcept {
    GravityPartials p;
    gravity<true>(r, body.J.data(), body.J.size(), body.mu, body.re, false, nullptr, &p);
    return p.jerk(v);
}

} // namespace astro::detail
