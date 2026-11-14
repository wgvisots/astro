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
#include <ostream>
#include <string>

#include "astro/core/error.h"
#include "astro/core/vec.h"

/// \file
/// 6x6 matrices for position/velocity covariances and state transition matrices. Element (i, j) is
/// row i, column j; rows and columns 0-2 are position, 3-5 velocity.
namespace astro {

/// 6x6 matrix, row-major: a(i, j) is row i, column j. Used as Mat3 is: a(i, j), a * b, a + b,
/// a.transpose(), Mat6::identity().
struct Mat6 {
    /// Elements, m[row][column]; zero-initialised.
    std::array<std::array<double, 6>, 6> m{};

    /// 6x6 identity matrix.
    [[nodiscard]] static constexpr Mat6 identity() noexcept {
        Mat6 a;
        for (std::size_t i = 0; i < 6; ++i) {
            a.m[i][i] = 1.0;
        }
        return a;
    }

    /// @name Element access
    /// @{

    /// Element at row r, column c (0-based, no bounds check).
    [[nodiscard]] constexpr double  operator()(std::size_t r, std::size_t c) const noexcept { return m[r][c]; }
    [[nodiscard]] constexpr double& operator()(std::size_t r, std::size_t c) noexcept { return m[r][c]; }
    /// @}

    /// The 3x3 block at block row br, block column bc (0 = position, 1 = velocity).
    [[nodiscard]] constexpr Mat3 block(std::size_t br, std::size_t bc) const noexcept {
        Mat3 b;
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 3; ++j) {
                b(i, j) = m[3 * br + i][3 * bc + j];
            }
        }
        return b;
    }

    /// Transpose.
    [[nodiscard]] constexpr Mat6 transpose() const noexcept {
        Mat6 t;
        for (std::size_t i = 0; i < 6; ++i) {
            for (std::size_t j = 0; j < 6; ++j) {
                t.m[i][j] = m[j][i];
            }
        }
        return t;
    }

    /// @name Compound assignment
    /// @{

    /// Element-wise in-place addition/subtraction, and scaling by a scalar.
    constexpr Mat6& operator+=(const Mat6& o) noexcept {
        for (std::size_t i = 0; i < 6; ++i) {
            for (std::size_t j = 0; j < 6; ++j) {
                m[i][j] += o.m[i][j];
            }
        }
        return *this;
    }
    constexpr Mat6& operator-=(const Mat6& o) noexcept {
        for (std::size_t i = 0; i < 6; ++i) {
            for (std::size_t j = 0; j < 6; ++j) {
                m[i][j] -= o.m[i][j];
            }
        }
        return *this;
    }
    constexpr Mat6& operator*=(double s) noexcept {
        for (auto& r : m) {
            for (auto& e : r) {
                e *= s;
            }
        }
        return *this;
    }
    /// @}

    friend constexpr bool operator==(const Mat6&, const Mat6&) noexcept = default;
};

/// @name Mat6 element-wise arithmetic
/// @{

/// Element-wise sum and difference; scaling by a scalar. (The result is a named local, so it is
/// returned without a further copy of the 36 elements.)
[[nodiscard]] constexpr Mat6 operator+(const Mat6& a, const Mat6& b) noexcept {
    Mat6 c  = a;
    c      += b;
    return c;
}
[[nodiscard]] constexpr Mat6 operator-(const Mat6& a, const Mat6& b) noexcept {
    Mat6 c  = a;
    c      -= b;
    return c;
}
[[nodiscard]] constexpr Mat6 operator*(const Mat6& a, double s) noexcept {
    Mat6 c  = a;
    c      *= s;
    return c;
}
[[nodiscard]] constexpr Mat6 operator*(double s, const Mat6& a) noexcept {
    Mat6 c  = a;
    c      *= s;
    return c;
}
/// @}

/// Matrix product a * b.
[[nodiscard]] constexpr Mat6 operator*(const Mat6& a, const Mat6& b) noexcept {
    Mat6 c;
    for (std::size_t i = 0; i < 6; ++i) {
        for (std::size_t k = 0; k < 6; ++k) {
            const double aik = a.m[i][k];
            for (std::size_t j = 0; j < 6; ++j) {
                c.m[i][j] += aik * b.m[k][j];
            }
        }
    }
    return c;
}

/// Writes the matrix row by row as "[[r0], [r1], ..., [r5]]" using the stream's current formatting.
inline std::ostream& operator<<(std::ostream& os, const Mat6& a) {
    os << '[';
    for (std::size_t i = 0; i < 6; ++i) {
        os << (i == 0 ? "[" : ", [");
        for (std::size_t j = 0; j < 6; ++j) {
            os << (j == 0 ? "" : ", ") << a(i, j);
        }
        os << ']';
    }
    return os << ']';
}

/// a * p * a^T, returned exactly symmetric: the upper triangle is computed and mirrored into the
/// lower one, so a covariance stays a covariance.
[[nodiscard]] constexpr Mat6 congruence(const Mat6& a, const Mat6& p) noexcept {
    const Mat6 ap = a * p;
    Mat6       c;
    for (std::size_t i = 0; i < 6; ++i) {
        for (std::size_t j = i; j < 6; ++j) {
            double s = 0.0;
            for (std::size_t k = 0; k < 6; ++k) {
                s += ap(i, k) * a(j, k);
            }
            c(i, j) = s;
            c(j, i) = s;
        }
    }
    return c;
}

/// Rotate a symmetric 6x6 position/velocity covariance by a (static) 3x3 frame rotation:
/// P' = diag(R, R) P diag(R, R)^T. Only the upper-right velocity block of P is read (P must be
/// symmetric). The result is symmetric to rounding, not exactly: the two diagonal blocks are plain
/// matrix products R P_ii R^T (the off-diagonal blocks are exact transposes of each other).
[[nodiscard]] inline Mat6 rotateCovariance(const Mat6& P, const Mat3& R) noexcept {
    // Block by block: P'_ij = R P_ij R^T (the lower-left block is the transpose of the upper-right).
    const Mat3 rr = R * P.block(0, 0) * R.transpose();
    const Mat3 rv = R * P.block(0, 1) * R.transpose();
    const Mat3 vv = R * P.block(1, 1) * R.transpose();
    Mat6       out;
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            out(i, j)         = rr(i, j);
            out(i, j + 3)     = rv(i, j);
            out(i + 3, j)     = rv(j, i);
            out(i + 3, j + 3) = vv(i, j);
        }
    }
    return out;
}

/// Squared Mahalanobis distance d^T S^-1 d of the first n components (n = 3 or 6) of d, for
/// symmetric positive-definite S (its leading n x n block). Uses a Cholesky factorization of S
/// scaled to unit diagonal, so mixed units (m^2 and m^2/s^2) do not hurt the conditioning.
/// nullopt if S is not positive definite or the distance is not finite (d holds NaN or inf).
/// Only the lower triangle of S is read. A matrix that is positive definite but within rounding
/// of singular (a correlation of 1 - 1e-15, say) can be reported as not positive definite.
/// Throws astro::Error if `n` is not 1 to 6.
[[nodiscard]] inline std::optional<double>
mahalanobis2(const std::array<double, 6>& d, const Mat6& S, std::size_t n = 6) {
    if (n == 0 || n > 6) {
        throw Error("mahalanobis2: n must be 1 to 6 (got " + std::to_string(n) + ")");
    }
    double scale[6];
    for (std::size_t i = 0; i < n; ++i) {
        if (!(S(i, i) > 0.0) || !std::isfinite(S(i, i))) {
            return std::nullopt;
        }
        scale[i] = 1.0 / std::sqrt(S(i, i));
    }
    // L L^T = D S D with D = diag(scale); then solve L y = D d and return |y|^2.
    double Lm[6][6] = {};
    double y[6]     = {};
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            double sum = S(i, j) * scale[i] * scale[j];
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
    if (!std::isfinite(d2)) {
        return std::nullopt;
    }
    return d2;
}

} // namespace astro
