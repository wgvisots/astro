/**********************************************************************
 * gravity.cpp
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

#include "astro/orbit/gravity.h"

#include <cmath>

#include "gravity_detail.h"

namespace astro {

Vec3 pointMassAcceleration(const Vec3& r, double mu) noexcept {
    const double rm = norm(r);
    return r * (-mu / (rm * rm * rm));
}

Vec3 pointMassJerk(const Vec3& r, const Vec3& v, double mu) noexcept {
    const double r2 = norm2(r), rm = std::sqrt(r2);
    return (v - r * (3.0 * dot(r, v) / r2)) * (-mu / (r2 * rm));
}

namespace {

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
constexpr LegendreFactors kLegendre;

} // namespace

namespace detail {

/// Zonal (and optionally point-mass) acceleration and, if WithJerk, its time derivative along v.
/// Allocation-free: P_n(s), P_n'(s) and P_n''(s) are carried in a rolling window of two degrees.
///
///   a = mu / r^2 [ (cr - pm) r_hat + cz z_hat ],   s = z / r,   k_n = Jn (re / r)^n,
///   cr = sum k_n ((n + 1) P_n + s P_n'),     cz = -sum k_n P_n',   pm = 1 with the point mass.
template <bool WithJerk>
void gravity(const Vec3&   r,
             const Vec3&   v,
             const double* J,
             std::size_t   size,
             double        mu,
             double        re,
             bool          pointMass,
             Vec3*         acc,
             Vec3*         jerk) noexcept {
    const double invR = 1.0 / std::sqrt(dot(r, r));
    const Vec3   rhat = r * invR;
    const double s    = rhat.z;
    const double q    = re * invR;

    double P2 = 1.0, P1 = s; // degrees n-2 and n-1 of P, P', P''
    double dP2 = 0.0, dP1 = 1.0;
    double ddP2 = 0.0, ddP1 = 0.0;
    double qn = q;
    double cr = pointMass ? -1.0 : 0.0, cz = 0.0, crDot = 0.0, czDot = 0.0;
    double rDot = 0.0, sDot = 0.0;
    if constexpr (WithJerk) {
        rDot = dot(rhat, v);
        sDot = (v.z - s * rDot) * invR;
    }
    for (std::size_t n = 2; n < size; ++n) {
        const double dn  = static_cast<double>(n);
        const double an  = n < LegendreFactors::size ? kLegendre.a[n] : (2.0 * dn - 1.0) / dn;
        const double bn  = n < LegendreFactors::size ? kLegendre.b[n] : (dn - 1.0) / dn;
        const double P   = an * s * P1 - bn * P2;
        const double dP  = dP2 + (2.0 * dn - 1.0) * P1;
        double       ddP = 0.0;
        if constexpr (WithJerk) {
            ddP = ddP2 + (2.0 * dn - 1.0) * dP1;
        }
        qn *= q;
        if (J[n] != 0.0) {
            const double k   = J[n] * qn;
            const double fr  = (dn + 1.0) * P + s * dP;
            cr              += k * fr;
            cz              -= k * dP;
            if constexpr (WithJerk) {
                const double kDot  = -dn * k * rDot * invR;
                crDot             += kDot * fr + k * sDot * ((dn + 2.0) * dP + s * ddP);
                czDot             -= kDot * dP + k * sDot * ddP;
            }
        }
        P2  = P1;
        P1  = P;
        dP2 = dP1;
        dP1 = dP;
        if constexpr (WithJerk) {
            ddP2 = ddP1;
            ddP1 = ddP;
        }
    }
    const double scale = mu * invR * invR;
    if (acc) {
        *acc = scale * (cr * rhat + Vec3{0.0, 0.0, cz});
    }
    if constexpr (WithJerk) {
        const double scaleDot = -2.0 * scale * rDot * invR;
        const Vec3   rhatDot  = (v - rhat * rDot) * invR;
        *jerk =
            scaleDot * (cr * rhat + Vec3{0.0, 0.0, cz}) + scale * (crDot * rhat + cr * rhatDot + Vec3{0.0, 0.0, czDot});
    }
}

template void
gravity<false>(const Vec3&, const Vec3&, const double*, std::size_t, double, double, bool, Vec3*, Vec3*) noexcept;
template void
gravity<true>(const Vec3&, const Vec3&, const double*, std::size_t, double, double, bool, Vec3*, Vec3*) noexcept;

} // namespace detail

Vec3 zonalAcceleration(const Vec3& r, const std::vector<double>& J, double mu, double re) noexcept {
    Vec3 a;
    detail::gravity<false>(r, Vec3{}, J.data(), J.size(), mu, re, false, &a, nullptr);
    return a;
}

Vec3 zonalJerk(const Vec3& r, const Vec3& v, const std::vector<double>& J, double mu, double re) noexcept {
    Vec3 j;
    detail::gravity<true>(r, v, J.data(), J.size(), mu, re, false, nullptr, &j);
    return j;
}

void zonalAccelerationAndJerk(const Vec3&                r,
                              const Vec3&                v,
                              const std::vector<double>& J,
                              Vec3&                      acceleration,
                              Vec3&                      jerk,
                              double                     mu,
                              double                     re) noexcept {
    detail::gravity<true>(r, v, J.data(), J.size(), mu, re, false, &acceleration, &jerk);
}

Vec3 j2Acceleration(const Vec3& r, double J2, double mu, double re) noexcept {
    const double rm   = norm(r);
    const double r2   = rm * rm;
    const double z2r2 = r.z * r.z / r2;
    const double k    = -1.5 * J2 * mu * re * re / (r2 * r2 * rm);
    return {k * r.x * (1.0 - 5.0 * z2r2), k * r.y * (1.0 - 5.0 * z2r2), k * r.z * (3.0 - 5.0 * z2r2)};
}

Vec3 thirdBodyAcceleration(const Vec3& r, const Vec3& s, double muBody) noexcept {
    const Vec3   d  = s - r;
    const double dm = norm(d), sm = norm(s);
    return muBody * (d / (dm * dm * dm) - s / (sm * sm * sm));
}

} // namespace astro
