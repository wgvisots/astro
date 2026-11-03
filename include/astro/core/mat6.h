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
#include <cmath>
#include <cstddef>
#include <optional>

#include "astro/core/vec.h"

/// \file
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

/// Squared Mahalanobis distance d^T S^-1 d of the first n components (n = 3 or 6) of d, for
/// symmetric positive-definite S (its leading n x n block). Uses a Cholesky factorization of S
/// scaled to unit diagonal, so mixed units (m^2 and m^2/s^2) do not hurt the conditioning.
/// nullopt if S is not positive definite.
[[nodiscard]] inline std::optional<double>
mahalanobis2(const std::array<double, 6>& d, const Mat6& S, std::size_t n = 6) noexcept {
    if (n == 0 || n > 6) {
        return std::nullopt;
    }
    double scale[6];
    for (std::size_t i = 0; i < n; ++i) {
        if (!(S[i][i] > 0.0) || !std::isfinite(S[i][i])) {
            return std::nullopt;
        }
        scale[i] = 1.0 / std::sqrt(S[i][i]);
    }
    // L L^T = D S D with D = diag(scale); then solve L y = D d and return |y|^2.
    double Lm[6][6] = {};
    double y[6]     = {};
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            double sum = S[i][j] * scale[i] * scale[j];
            for (std::size_t k = 0; k < j; ++k) {
                sum -= Lm[i][k] * Lm[j][k];
            }
            if (i == j) {
                if (!(sum > 1e-14)) {
                    return std::nullopt;
                }
                Lm[i][i] = std::sqrt(sum);
            } else {
                Lm[i][j] = sum / Lm[j][j];
            }
        }
    }
    double d2 = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        double sum = d[i] * scale[i];
        for (std::size_t k = 0; k < i; ++k) {
            sum -= Lm[i][k] * y[k];
        }
        y[i]  = sum / Lm[i][i];
        d2   += y[i] * y[i];
    }
    return d2;
}

} // namespace astro
