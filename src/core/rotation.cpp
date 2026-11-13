/**********************************************************************
 * rotation.cpp
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

#include "astro/core/rotation.h"

#include <algorithm>
#include <cmath>

#include "../core/enum_io.h"
#include "astro/core/error.h"

namespace astro {

// clang-format off
Mat3 R1(double theta) noexcept {
    const double c = std::cos(theta), s = std::sin(theta);
    return {1.0, 0.0, 0.0,
            0.0,   c,   s,
            0.0,  -s,   c};
}

Mat3 R2(double theta) noexcept {
    const double c = std::cos(theta), s = std::sin(theta);
    return {  c, 0.0,  -s,
            0.0, 1.0, 0.0,
              s, 0.0,   c};
}

Mat3 R3(double theta) noexcept {
    const double c = std::cos(theta), s = std::sin(theta);
    return {  c,   s, 0.0,
             -s,   c, 0.0,
            0.0, 0.0, 1.0};
}
// clang-format on

Mat3 rotationAboutAxis(int axis, double theta) {
    switch (axis) {
    case 1:  return R1(theta);
    case 2:  return R2(theta);
    case 3:  return R3(theta);
    default: throw Error("rotationAboutAxis: axis must be 1, 2 or 3");
    }
}

Mat3 axisAngleToMatrix(const Vec3& axis, double theta) noexcept {
    return Quat::fromAxisAngle(axis, theta).toMatrix();
}

// ------------------------------------------------------------ Quat

Quat Quat::fromAxisAngle(const Vec3& axis, double theta) noexcept {
    const Vec3   n = astro::normalized(axis);
    const double s = std::sin(0.5 * theta);
    return {std::cos(0.5 * theta), s * n.x, s * n.y, s * n.z};
}

Mat3 Quat::toMatrix() const noexcept {
    // Passive DCM: A = (w^2 - |v|^2) I + 2 v v^T - 2 w [v x]
    const double ww = w * w, xx = x * x, yy = y * y, zz = z * z;
    const double xy = x * y, xz = x * z, yz = y * z;
    const double wx = w * x, wy = w * y, wz = w * z;
    // clang-format off
    return {ww + xx - yy - zz, 2.0 * (xy + wz),   2.0 * (xz - wy),
            2.0 * (xy - wz),   ww - xx + yy - zz, 2.0 * (yz + wx),
            2.0 * (xz + wy),   2.0 * (yz - wx),   ww - xx - yy + zz};
    // clang-format on
}

Quat Quat::fromMatrix(const Mat3& a) noexcept {
    // Shepperd's method: choose the largest of 4w^2, 4x^2, 4y^2, 4z^2.
    const double                tr = a.trace();
    const std::array<double, 4> d{tr, a(0, 0), a(1, 1), a(2, 2)};
    const auto                  idx = static_cast<std::size_t>(std::max_element(d.begin(), d.end()) - d.begin());
    Quat                        q;
    switch (idx) {
    case 0: {
        const double s = 2.0 * std::sqrt(1.0 + tr); // 4w
        q              = {0.25 * s, (a(1, 2) - a(2, 1)) / s, (a(2, 0) - a(0, 2)) / s, (a(0, 1) - a(1, 0)) / s};
        break;
    }
    case 1: {
        const double s = 2.0 * std::sqrt(1.0 + 2.0 * a(0, 0) - tr); // 4x
        q              = {(a(1, 2) - a(2, 1)) / s, 0.25 * s, (a(0, 1) + a(1, 0)) / s, (a(0, 2) + a(2, 0)) / s};
        break;
    }
    case 2: {
        const double s = 2.0 * std::sqrt(1.0 + 2.0 * a(1, 1) - tr); // 4y
        q              = {(a(2, 0) - a(0, 2)) / s, (a(0, 1) + a(1, 0)) / s, 0.25 * s, (a(1, 2) + a(2, 1)) / s};
        break;
    }
    default: {
        const double s = 2.0 * std::sqrt(1.0 + 2.0 * a(2, 2) - tr); // 4z
        q              = {(a(0, 1) - a(1, 0)) / s, (a(0, 2) + a(2, 0)) / s, (a(1, 2) + a(2, 1)) / s, 0.25 * s};
        break;
    }
    }
    if (q.w < 0.0) {
        q = {-q.w, -q.x, -q.y, -q.z};
    }
    return q.normalized();
}

double Quat::norm() const noexcept {
    return std::sqrt(w * w + x * x + y * y + z * z);
}

Quat Quat::normalized() const noexcept {
    const double n = norm();
    // Test for zero, not for n > 0: a NaN norm fails every comparison and must not be taken for zero.
    return n == 0.0 ? Quat{} : Quat{w / n, x / n, y / n, z / n};
}

Vec3 Quat::transform(const Vec3& v) const noexcept {
    // A v = q* (0,v) q
    const Quat p = conjugate() * Quat{0.0, v.x, v.y, v.z} * (*this);
    return p.vec();
}

Vec3 Quat::rotate(const Vec3& v) const noexcept {
    const Quat p = (*this) * Quat{0.0, v.x, v.y, v.z} * conjugate();
    return p.vec();
}

Quat slerp(const Quat& a, const Quat& b0, double t) noexcept {
    Quat   b = b0;
    double d = a.w * b.w + a.x * b.x + a.y * b.y + a.z * b.z;
    if (d < 0.0) {
        b = {-b.w, -b.x, -b.y, -b.z};
        d = -d;
    }
    double wa, wb;
    if (d > 0.9995) {
        wa = 1.0 - t;
        wb = t;
    } else {
        const double theta = std::acos(std::clamp(d, -1.0, 1.0));
        const double s     = std::sin(theta);
        wa                 = std::sin((1.0 - t) * theta) / s;
        wb                 = std::sin(t * theta) / s;
    }
    return Quat{wa * a.w + wb * b.w, wa * a.x + wb * b.x, wa * a.y + wb * b.y, wa * a.z + wb * b.z}.normalized();
}

// ------------------------------------------------------------ Euler

std::array<int, 3> eulerAxes(EulerSequence seq) noexcept {
    switch (seq) {
    case EulerSequence::XYZ: return {1, 2, 3};
    case EulerSequence::XZY: return {1, 3, 2};
    case EulerSequence::YXZ: return {2, 1, 3};
    case EulerSequence::YZX: return {2, 3, 1};
    case EulerSequence::ZXY: return {3, 1, 2};
    case EulerSequence::ZYX: return {3, 2, 1};
    case EulerSequence::XYX: return {1, 2, 1};
    case EulerSequence::XZX: return {1, 3, 1};
    case EulerSequence::YXY: return {2, 1, 2};
    case EulerSequence::YZY: return {2, 3, 2};
    case EulerSequence::ZXZ: return {3, 1, 3};
    case EulerSequence::ZYZ: return {3, 2, 3};
    }
    return {3, 2, 1};
}

Mat3 eulerToMatrix(EulerSequence seq, const EulerAngles& e) noexcept {
    const auto ax = eulerAxes(seq);
    return rotationAboutAxis(ax[2], e.c) * rotationAboutAxis(ax[1], e.b) * rotationAboutAxis(ax[0], e.a);
}

EulerAngles matrixToEuler(EulerSequence seq, const Mat3& dcm) noexcept {
    // Work with the active matrix M = A^T = Rot_i(a) Rot_j(b) Rot_k(c).
    const Mat3        M      = dcm.transpose();
    const auto        ax     = eulerAxes(seq);
    const std::size_t i      = static_cast<std::size_t>(ax[0] - 1);
    const std::size_t j      = static_cast<std::size_t>(ax[1] - 1);
    const bool        proper = ax[0] == ax[2];
    // Third axis distinct from i and j.
    const std::size_t k = proper ? 3 - i - j : static_cast<std::size_t>(ax[2] - 1);
    // +1 if (i, j, k) is a cyclic permutation of (0, 1, 2), else -1.
    const double     s       = ((j == (i + 1) % 3) ? 1.0 : -1.0);
    constexpr double lockTol = 1e-12;

    EulerAngles e;
    if (!proper) {
        const double sb = std::clamp(s * M(i, k), -1.0, 1.0);
        e.b             = std::asin(sb);
        if (1.0 - std::abs(sb) > lockTol) {
            e.a = std::atan2(-s * M(j, k), M(k, k));
            e.c = std::atan2(-s * M(i, j), M(i, i));
        } else {
            e.a = std::atan2(s * M(k, j), M(j, j));
            e.c = 0.0;
        }
    } else {
        const double cb = std::clamp(M(i, i), -1.0, 1.0);
        e.b             = std::acos(cb);
        if (1.0 - std::abs(cb) > lockTol) {
            e.a = std::atan2(M(j, i), -s * M(k, i));
            e.c = std::atan2(M(i, j), s * M(i, k));
        } else {
            e.a = std::atan2(s * M(k, j), M(j, j));
            e.c = 0.0;
        }
    }
    return e;
}

Mat3 orthonormalize(const Mat3& m) noexcept {
    return Quat::fromMatrix(m).toMatrix();
}

std::ostream& operator<<(std::ostream& os, EulerSequence value) {
    return detail::writeEnum(os, "EulerSequence", value);
}

} // namespace astro
