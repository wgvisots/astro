/**********************************************************************
 * mat6.h
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

#include <array>
#include <cstddef>

#include "astro/core/vec.h"

/// 6x6 matrices for position/velocity covariances and state transition matrices. Element [i][j]
/// is row i, column j; rows and columns 0-2 are position, 3-5 velocity.
///
/// Mat6 is a std::array, so the operations are named functions (operators declared here would not
/// be found by argument-dependent lookup outside namespace astro).
namespace astro {

using Mat6 = std::array<std::array<double, 6>, 6>;

[[nodiscard]] constexpr Mat6 identity6() noexcept {
    Mat6 m{};
    for (std::size_t i = 0; i < 6; ++i) {
        m[i][i] = 1.0;
    }
    return m;
}

[[nodiscard]] constexpr Mat6 transpose(const Mat6& a) noexcept {
    Mat6 t{};
    for (std::size_t i = 0; i < 6; ++i) {
        for (std::size_t j = 0; j < 6; ++j) {
            t[i][j] = a[j][i];
        }
    }
    return t;
}

/// a * b
[[nodiscard]] constexpr Mat6 multiply(const Mat6& a, const Mat6& b) noexcept {
    Mat6 c{};
    for (std::size_t i = 0; i < 6; ++i) {
        for (std::size_t k = 0; k < 6; ++k) {
            const double aik = a[i][k];
            for (std::size_t j = 0; j < 6; ++j) {
                c[i][j] += aik * b[k][j];
            }
        }
    }
    return c;
}

/// a + b
[[nodiscard]] constexpr Mat6 add(const Mat6& a, const Mat6& b) noexcept {
    Mat6 c = a;
    for (std::size_t i = 0; i < 6; ++i) {
        for (std::size_t j = 0; j < 6; ++j) {
            c[i][j] += b[i][j];
        }
    }
    return c;
}

/// a * p * a^T, symmetrized ((m + m^T) / 2) so that a covariance stays exactly symmetric.
[[nodiscard]] constexpr Mat6 congruence(const Mat6& a, const Mat6& p) noexcept {
    const Mat6 ap = multiply(a, p);
    Mat6       c{};
    for (std::size_t i = 0; i < 6; ++i) {
        for (std::size_t j = i; j < 6; ++j) {
            double s = 0.0;
            for (std::size_t k = 0; k < 6; ++k) {
                s += ap[i][k] * a[j][k];
            }
            c[i][j] = s;
            c[j][i] = s;
        }
    }
    return c;
}

/// The 3x3 block at block row br, block column bc (0 = position, 1 = velocity).
[[nodiscard]] constexpr Mat3 block(const Mat6& a, std::size_t br, std::size_t bc) noexcept {
    Mat3 m;
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            m(i, j) = a[3 * br + i][3 * bc + j];
        }
    }
    return m;
}

} // namespace astro
