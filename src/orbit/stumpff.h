/**********************************************************************
 * stumpff.h
 *
 * Copyright 2024-2026 W. Greg Visotski
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

// The Stumpff functions of the universal-variable Kepler propagator. Header-only, so that the
// library and the tests share one definition without exporting it. Internal: not installed.
namespace astro::detail {

/// The two Stumpff functions of one argument, as stumpff() returns them.
struct Stumpff {
    double c2{0.0}; ///< c2(psi), dimensionless
    double c3{0.0}; ///< c3(psi), dimensionless
};

constexpr double reciprocalFactorial(int n) noexcept {
    double f = 1.0;
    for (int i = 2; i <= n; ++i) {
        f *= i;
    }
    return 1.0 / f;
}

/// Stumpff functions c2(psi) = (1 - cos sqrt(psi)) / psi and c3(psi) = (sqrt(psi) - sin sqrt(psi)) / psi^(3/2),
/// continued to psi <= 0 (hyperbolic functions). `psi` is dimensionless (chi^2 alpha in the universal
/// variables). Accurate to a few rounding errors for every psi: by their power series for |psi| < 1, where the
/// closed forms cancel.
[[nodiscard]] inline Stumpff stumpff(double psi) noexcept {
    Stumpff s;
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
        s.c2 = s2;
        s.c3 = s3;
    } else if (psi > 0.0) {
        const double sp = std::sqrt(psi);
        s.c2            = (1.0 - std::cos(sp)) / psi;
        s.c3            = (sp - std::sin(sp)) / (sp * psi);
    } else {
        const double sp = std::sqrt(-psi);
        s.c2            = (1.0 - std::cosh(sp)) / psi;
        s.c3            = (std::sinh(sp) - sp) / (sp * -psi);
    }
    return s;
}

} // namespace astro::detail
