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

namespace astro {

Vec3 pointMassAcceleration(const Vec3& r, double mu) noexcept {
    const double rm = norm(r);
    return r * (-mu / (rm * rm * rm));
}

Vec3 zonalAcceleration(const Vec3& r, const std::vector<double>& J, double mu, double re) noexcept {
    // a_n = mu Jn re^n / r^(n+2) [ ((n+1) Pn(s) + s Pn'(s)) r_hat - Pn'(s) z_hat ],  s = z / r
    const double rm = norm(r);
    const Vec3 rhat = r / rm;
    const double s = rhat.z;
    const std::size_t nmax = J.size() > 0 ? J.size() - 1 : 0;
    if (nmax < 2) return {};
    // Legendre polynomials and derivatives by recursion.
    std::vector<double> P(nmax + 1), dP(nmax + 1);
    P[0] = 1.0;
    P[1] = s;
    dP[0] = 0.0;
    dP[1] = 1.0;
    for (std::size_t n = 2; n <= nmax; ++n) {
        const double dn = static_cast<double>(n);
        P[n] = ((2.0 * dn - 1.0) * s * P[n - 1] - (dn - 1.0) * P[n - 2]) / dn;
        dP[n] = dP[n - 2] + (2.0 * dn - 1.0) * P[n - 1];
    }
    double cr = 0.0, cz = 0.0;
    const double q = re / rm;
    double qn = q;  // (re / r)^n, built up
    for (std::size_t n = 2; n <= nmax; ++n) {
        qn *= q;
        if (J[n] == 0.0) continue;
        const double k = J[n] * qn;
        cr += k * ((static_cast<double>(n) + 1.0) * P[n] + s * dP[n]);
        cz -= k * dP[n];
    }
    const double scale = mu / (rm * rm);
    return scale * (cr * rhat + Vec3{0.0, 0.0, cz});
}

Vec3 j2Acceleration(const Vec3& r, double J2, double mu, double re) noexcept {
    const double rm = norm(r);
    const double r2 = rm * rm;
    const double z2r2 = r.z * r.z / r2;
    const double k = -1.5 * J2 * mu * re * re / (r2 * r2 * rm);
    return {k * r.x * (1.0 - 5.0 * z2r2), k * r.y * (1.0 - 5.0 * z2r2), k * r.z * (3.0 - 5.0 * z2r2)};
}

Vec3 thirdBodyAcceleration(const Vec3& r, const Vec3& s, double muBody) noexcept {
    const Vec3 d = s - r;
    const double dm = norm(d), sm = norm(s);
    return muBody * (d / (dm * dm * dm) - s / (sm * sm * sm));
}

} // namespace astro
