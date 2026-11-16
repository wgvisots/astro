/**********************************************************************
 * chi_square.h
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

#include "astro/core/constants.h"

// Chi-square tail probabilities for the consistency tests of EstimateTracker. Header-only, so that
// the library and the tests share one definition without exporting it. Internal: not installed.
namespace astro::detail {

/// Survival function of the chi-square distribution with 3 and 6 degrees of freedom, P(X >= x).
/// 1 for x <= 0; 0 for +infinity and for any x so large that the value is below the smallest double
/// (x > 1490). NaN for NaN, so a NaN statistic never passes a p-value test.
[[nodiscard]] inline double chiSquareSurvival3(double x) noexcept {
    if (std::isnan(x)) {
        return x;
    }
    if (!(x > 0.0)) {
        return 1.0;
    }
    // Once exp(-x/2) underflows (x > 1490) the tail is below the smallest double. Returning 0 there
    // also keeps a huge or infinite x from giving infinity * 0 = NaN.
    const double e = std::exp(-0.5 * x);
    if (e == 0.0) {
        return 0.0;
    }
    return std::erfc(std::sqrt(0.5 * x)) + std::sqrt(2.0 * x / constants::pi) * e;
}

[[nodiscard]] inline double chiSquareSurvival6(double x) noexcept {
    if (std::isnan(x)) {
        return x;
    }
    if (!(x > 0.0)) {
        return 1.0;
    }
    const double h = 0.5 * x;
    const double e = std::exp(-h);
    if (e == 0.0) { // as in chiSquareSurvival3: the polynomial overflows from x = 2.7e154
        return 0.0;
    }
    return e * (1.0 + h + 0.5 * h * h);
}

} // namespace astro::detail
