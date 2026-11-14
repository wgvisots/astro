/**********************************************************************
 * vec.h
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
#include <ostream>

/// \file
/// Vec3, Mat3 and the position/velocity state vectors StateVector and StateVector12.
namespace astro {

/// Three-component Cartesian vector (double precision).
struct Vec3 {
    double x{0.0}; ///< x component
    double y{0.0}; ///< y component
    double z{0.0}; ///< z component

    constexpr Vec3() noexcept = default;
    /// Vector with components (x_, y_, z_).
    constexpr Vec3(double x_, double y_, double z_) noexcept
        : x(x_)
        , y(y_)
        , z(z_) {}

    /// @name Component access
    /// @{

    /// Component i: 0 = x, 1 = y, any other index = z (no bounds check).
    [[nodiscard]] constexpr double  operator[](std::size_t i) const noexcept { return i == 0 ? x : (i == 1 ? y : z); }
    [[nodiscard]] constexpr double& operator[](std::size_t i) noexcept { return i == 0 ? x : (i == 1 ? y : z); }
    /// @}

    /// @name Compound assignment
    /// @{

    /// Component-wise in-place addition/subtraction, and scaling by a scalar.
    constexpr Vec3& operator+=(const Vec3& o) noexcept {
        x += o.x;
        y += o.y;
        z += o.z;
        return *this;
    }
    constexpr Vec3& operator-=(const Vec3& o) noexcept {
        x -= o.x;
        y -= o.y;
        z -= o.z;
        return *this;
    }
    constexpr Vec3& operator*=(double s) noexcept {
        x *= s;
        y *= s;
        z *= s;
        return *this;
    }
    constexpr Vec3& operator/=(double s) noexcept {
        x /= s;
        y /= s;
        z /= s;
        return *this;
    }
    /// @}

    friend constexpr bool operator==(const Vec3&, const Vec3&) noexcept = default;
};

/// @name Vec3 arithmetic
/// @{

/// Component-wise sum, difference and negation; scaling by a scalar.
[[nodiscard]] constexpr Vec3 operator+(Vec3 a, const Vec3& b) noexcept {
    return a += b;
}
[[nodiscard]] constexpr Vec3 operator-(Vec3 a, const Vec3& b) noexcept {
    return a -= b;
}
[[nodiscard]] constexpr Vec3 operator-(const Vec3& a) noexcept {
    return {-a.x, -a.y, -a.z};
}
[[nodiscard]] constexpr Vec3 operator*(Vec3 a, double s) noexcept {
    return a *= s;
}
[[nodiscard]] constexpr Vec3 operator*(double s, Vec3 a) noexcept {
    return a *= s;
}
[[nodiscard]] constexpr Vec3 operator/(Vec3 a, double s) noexcept {
    return a /= s;
}
/// @}

/// Dot (scalar) product a.b.
[[nodiscard]] constexpr double dot(const Vec3& a, const Vec3& b) noexcept {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
/// Cross product a x b (right-handed).
[[nodiscard]] constexpr Vec3 cross(const Vec3& a, const Vec3& b) noexcept {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
/// Squared Euclidean length a.a.
[[nodiscard]] constexpr double norm2(const Vec3& a) noexcept {
    return dot(a, a);
}
/// Euclidean length. Uses sqrt(a.a) rather than std::hypot: the three-argument hypot guards
/// against overflow far outside any physical range here and is several times slower (5.5x measured).
[[nodiscard]] inline double norm(const Vec3& a) noexcept {
    return std::sqrt(dot(a, a));
}
/// Unit vector along a. Returns the zero vector if a is zero. A not-a-number (NaN) component makes
/// every component of the result NaN. Like norm(), this squares the components, so it is meant for
/// lengths between about 1e-154 and 1e154: beyond that the sum of squares underflows or overflows and
/// the result is the zero vector or is not of unit length.
[[nodiscard]] inline Vec3 normalized(const Vec3& a) noexcept {
    const double n = norm(a);
    // Test for zero, not for n > 0: a NaN length fails every comparison and must not be taken for zero.
    return n == 0.0 ? Vec3{} : a / n;
}
/// Angle between two vectors in [0, pi] (numerically robust atan2 form).
[[nodiscard]] inline double angleBetween(const Vec3& a, const Vec3& b) noexcept {
    return std::atan2(norm(cross(a, b)), dot(a, b));
}

/// Writes the vector as "[x, y, z]" using the stream's current formatting.
inline std::ostream& operator<<(std::ostream& os, const Vec3& v) {
    return os << '[' << v.x << ", " << v.y << ", " << v.z << ']';
}

/// Six-component state vector: position and velocity.
struct StateVector {
    Vec3                  r; ///< position
    Vec3                  v; ///< velocity
    friend constexpr bool operator==(const StateVector&, const StateVector&) noexcept = default;
};

/// Twelve-component state vector: position, velocity, acceleration and jerk (e.g. m, m/s, m/s^2,
/// m/s^3). All four must be given when constructing one, so `{r, v}` always means a StateVector.
struct StateVector12 {
    Vec3 r; ///< position
    Vec3 v; ///< velocity
    Vec3 a; ///< acceleration
    Vec3 j; ///< jerk (time derivative of acceleration)

    constexpr StateVector12() noexcept = default;
    /// State from position, velocity, acceleration and jerk.
    constexpr StateVector12(const Vec3& r_, const Vec3& v_, const Vec3& a_, const Vec3& j_) noexcept
        : r(r_)
        , v(v_)
        , a(a_)
        , j(j_) {}

    /// Position and velocity only.
    [[nodiscard]] constexpr StateVector rv() const noexcept { return {r, v}; }
    friend constexpr bool               operator==(const StateVector12&, const StateVector12&) noexcept = default;
};

/// 3x3 matrix, row-major.
struct Mat3 {
    // clang-format off
    /// Elements, m[row][column]; zero-initialised.
    std::array<std::array<double, 3>, 3> m{{{0.0, 0.0, 0.0},
                                            {0.0, 0.0, 0.0},
                                            {0.0, 0.0, 0.0}}};

    constexpr Mat3() noexcept = default;
    /// Matrix from its nine elements in row-major order (aRC = row R, column C).
    constexpr Mat3(double a00, double a01, double a02,
                   double a10, double a11, double a12,
                   double a20, double a21, double a22) noexcept
        : m{{{a00, a01, a02},
             {a10, a11, a12},
             {a20, a21, a22}}} {}

    /// 3x3 identity matrix.
    [[nodiscard]] static constexpr Mat3 identity() noexcept {
        return {1.0, 0.0, 0.0,
                0.0, 1.0, 0.0,
                0.0, 0.0, 1.0};
    }
    /// Matrix whose rows are r0, r1, r2.
    [[nodiscard]] static constexpr Mat3 fromRows(const Vec3& r0, const Vec3& r1, const Vec3& r2) noexcept {
        return {r0.x, r0.y, r0.z,
                r1.x, r1.y, r1.z,
                r2.x, r2.y, r2.z};
    }
    /// Matrix whose columns are c0, c1, c2.
    [[nodiscard]] static constexpr Mat3 fromCols(const Vec3& c0, const Vec3& c1, const Vec3& c2) noexcept {
        return {c0.x, c1.x, c2.x,
                c0.y, c1.y, c2.y,
                c0.z, c1.z, c2.z};
    }
    /// Skew-symmetric cross-product matrix: skew(a) * b == cross(a, b).
    [[nodiscard]] static constexpr Mat3 skew(const Vec3& a) noexcept {
        return { 0.0, -a.z,  a.y,
                 a.z,  0.0, -a.x,
                -a.y,  a.x,  0.0};
    }
    // clang-format on

    /// @name Element access
    /// @{

    /// Element at row r, column c (0-based, no bounds check).
    [[nodiscard]] constexpr double  operator()(std::size_t r, std::size_t c) const noexcept { return m[r][c]; }
    [[nodiscard]] constexpr double& operator()(std::size_t r, std::size_t c) noexcept { return m[r][c]; }
    /// @}

    /// Row r (0-based) as a vector.
    [[nodiscard]] constexpr Vec3 row(std::size_t r) const noexcept { return {m[r][0], m[r][1], m[r][2]}; }
    /// Column c (0-based) as a vector.
    [[nodiscard]] constexpr Vec3 col(std::size_t c) const noexcept { return {m[0][c], m[1][c], m[2][c]}; }

    /// Transpose (the inverse, for a rotation matrix).
    [[nodiscard]] constexpr Mat3 transpose() const noexcept {
        return {m[0][0], m[1][0], m[2][0], m[0][1], m[1][1], m[2][1], m[0][2], m[1][2], m[2][2]};
    }
    /// Sum of the diagonal elements.
    [[nodiscard]] constexpr double trace() const noexcept { return m[0][0] + m[1][1] + m[2][2]; }
    /// Determinant (+1 for a proper rotation matrix).
    [[nodiscard]] constexpr double determinant() const noexcept {
        return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
             + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    }

    /// @name Compound assignment
    /// @{

    /// Element-wise in-place addition/subtraction, and scaling by a scalar.
    constexpr Mat3& operator+=(const Mat3& o) noexcept {
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 3; ++j) {
                m[i][j] += o.m[i][j];
            }
        }
        return *this;
    }
    constexpr Mat3& operator-=(const Mat3& o) noexcept {
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 3; ++j) {
                m[i][j] -= o.m[i][j];
            }
        }
        return *this;
    }
    constexpr Mat3& operator*=(double s) noexcept {
        for (auto& r : m) {
            for (auto& e : r) {
                e *= s;
            }
        }
        return *this;
    }
    /// @}

    friend constexpr bool operator==(const Mat3&, const Mat3&) noexcept = default;
};

/// @name Mat3 element-wise arithmetic
/// @{

/// Element-wise sum and difference; scaling by a scalar.
[[nodiscard]] constexpr Mat3 operator+(Mat3 a, const Mat3& b) noexcept {
    return a += b;
}
[[nodiscard]] constexpr Mat3 operator-(Mat3 a, const Mat3& b) noexcept {
    return a -= b;
}
[[nodiscard]] constexpr Mat3 operator*(Mat3 a, double s) noexcept {
    return a *= s;
}
[[nodiscard]] constexpr Mat3 operator*(double s, Mat3 a) noexcept {
    return a *= s;
}
/// @}

/// Matrix product a * b.
[[nodiscard]] constexpr Mat3 operator*(const Mat3& a, const Mat3& b) noexcept {
    Mat3 c;
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            c.m[i][j] = a.m[i][0] * b.m[0][j] + a.m[i][1] * b.m[1][j] + a.m[i][2] * b.m[2][j];
        }
    }
    return c;
}

/// Matrix-vector product a * v (v as a column vector).
[[nodiscard]] constexpr Vec3 operator*(const Mat3& a, const Vec3& v) noexcept {
    return {a.m[0][0] * v.x + a.m[0][1] * v.y + a.m[0][2] * v.z,
            a.m[1][0] * v.x + a.m[1][1] * v.y + a.m[1][2] * v.z,
            a.m[2][0] * v.x + a.m[2][1] * v.y + a.m[2][2] * v.z};
}

/// Outer product a * b^T.
[[nodiscard]] constexpr Mat3 outer(const Vec3& a, const Vec3& b) noexcept {
    // clang-format off
    return {a.x * b.x, a.x * b.y, a.x * b.z,
            a.y * b.x, a.y * b.y, a.y * b.z,
            a.z * b.x, a.z * b.y, a.z * b.z};
    // clang-format on
}

/// Writes the matrix row by row as "[[r0], [r1], [r2]]".
inline std::ostream& operator<<(std::ostream& os, const Mat3& a) {
    return os << '[' << a.row(0) << ", " << a.row(1) << ", " << a.row(2) << ']';
}

} // namespace astro
