/**********************************************************************
 * rotation.h
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
#include <iosfwd>

#include <boost/describe/enum.hpp>

#include "astro/core/vec.h"

/// \file
/// Rotation conventions
/// --------------------
/// All matrices produced here are *passive* (frame / coordinate) rotations,
/// following the Explanatory Supplement and SOFA: R3(theta) rotates the
/// coordinate frame by +theta about its z axis, so a fixed vector's components
/// in the new frame are v' = R3(theta) * v. A passive matrix is the transpose of
/// the corresponding active (vector-rotating) matrix.
namespace astro {

/// Passive rotation about the x axis.
[[nodiscard]] Mat3 R1(double theta) noexcept;
/// Passive rotation about the y axis.
[[nodiscard]] Mat3 R2(double theta) noexcept;
/// Passive rotation about the z axis.
[[nodiscard]] Mat3 R3(double theta) noexcept;
/// Passive rotation about coordinate axis `axis` (1, 2 or 3).
[[nodiscard]] Mat3 rotationAboutAxis(int axis, double theta);

/// Passive rotation by `theta` about an arbitrary axis. Only the direction of `axis` is used: it may
/// have any length except zero. A zero axis defines no rotation and gives the identity matrix.
/// Not-a-number (NaN) or infinity in `axis` or `theta` gives a matrix of NaN.
[[nodiscard]] Mat3 axisAngleToMatrix(const Vec3& axis, double theta) noexcept;

/// Unit quaternion, scalar first: q = w + x i + y j + z k (Hamilton algebra).
///
/// The quaternion cos(t/2) + sin(t/2) n describes a frame rotation by angle t
/// about axis n, consistent with the passive matrices above:
/// Quat::fromAxisAngle({0,0,1}, t).toMatrix() == R3(t).
///
/// Composition: for passive matrices, toMatrix(p * q) == toMatrix(q) * toMatrix(p),
/// i.e. `p * q` means "rotate the frame by p, then by q (about the new axes)".
struct Quat {
    double w{1.0}; ///< scalar part
    double x{0.0}; ///< i component of the vector part
    double y{0.0}; ///< j component of the vector part
    double z{0.0}; ///< k component of the vector part

    constexpr Quat() noexcept = default;
    /// Quaternion w + x i + y j + z k, stored as given (not normalized).
    constexpr Quat(double w_, double x_, double y_, double z_) noexcept
        : w(w_)
        , x(x_)
        , y(y_)
        , z(z_) {}

    /// Frame rotation by `theta` (rad) about `axis`. Only the direction of `axis` is used (it is
    /// normalized internally): it may have any length except zero. A zero axis defines no rotation
    /// and gives the identity quaternion. Not-a-number (NaN) or infinity in `axis` or `theta` gives NaN.
    [[nodiscard]] static Quat fromAxisAngle(const Vec3& axis, double theta) noexcept;
    /// Quaternion from a passive rotation (direction cosine) matrix. Uses Shepperd's
    /// method; the result has w >= 0.
    [[nodiscard]] static Quat fromMatrix(const Mat3& dcm) noexcept;

    /// Passive rotation (direction cosine) matrix.
    [[nodiscard]] Mat3 toMatrix() const noexcept;

    /// Vector part (x, y, z).
    [[nodiscard]] constexpr Vec3 vec() const noexcept { return {x, y, z}; }
    /// Conjugate w - x i - y j - z k (the inverse rotation, for a unit quaternion).
    [[nodiscard]] constexpr Quat conjugate() const noexcept { return {w, -x, -y, -z}; }
    /// Euclidean norm sqrt(w^2 + x^2 + y^2 + z^2).
    [[nodiscard]] double norm() const noexcept;
    /// Unit quaternion along this one; the identity if the norm is zero. A not-a-number (NaN)
    /// component makes every component of the result NaN.
    [[nodiscard]] Quat normalized() const noexcept;

    /// Coordinates of a fixed vector in the rotated frame: toMatrix() * v.
    [[nodiscard]] Vec3 transform(const Vec3& v) const noexcept;
    /// Actively rotate vector v (in a fixed frame): toMatrix().transpose() * v.
    [[nodiscard]] Vec3 rotate(const Vec3& v) const noexcept;

    friend constexpr bool operator==(const Quat&, const Quat&) noexcept = default;
};

/// Hamilton product.
[[nodiscard]] constexpr Quat operator*(const Quat& a, const Quat& b) noexcept {
    return {a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z,
            a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
            a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x,
            a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w};
}

/// Spherical linear interpolation between unit quaternions (shortest path), t in [0, 1]: the rotation
/// that turns from a to b at a constant rate, a fraction t of the way along. Accurate to rounding
/// however small or large the angle between a and b.
[[nodiscard]] Quat slerp(const Quat& a, const Quat& b, double t) noexcept;

/// Euler angle sequence: the axes (1=x, 2=y, 3=z) of the three successive
/// frame rotations. E.g. ZYX = 3-2-1 (yaw, pitch, roll); ZXZ = 3-1-3 (classic
/// orbital Omega, i, omega).
enum class EulerSequence {
    XYZ,
    XZY,
    YXZ,
    YZX,
    ZXY,
    ZYX, // Tait-Bryan
    XYX,
    XZX,
    YXY,
    YZY,
    ZXZ,
    ZYZ, // proper Euler
};

BOOST_DESCRIBE_ENUM(EulerSequence, XYZ, XZY, YXZ, YZX, ZXY, ZYX, XYX, XZX, YXY, YZY, ZXZ, ZYZ)

/// Writes "EulerSequence::<enumerator>", or "EulerSequence::INVALID_VALUE" for any other value.
std::ostream& operator<<(std::ostream& os, EulerSequence value);

/// Axes of the sequence as 1-based axis numbers.
[[nodiscard]] std::array<int, 3> eulerAxes(EulerSequence seq) noexcept;

/// Euler angles (radians) for a sequence: first rotation by a, then b, then c.
struct EulerAngles {
    double a{0.0}; ///< first rotation angle, rad
    double b{0.0}; ///< second rotation angle, rad
    double c{0.0}; ///< third rotation angle, rad
};

/// Passive DCM for successive frame rotations: R_k(c) * R_j(b) * R_i(a)
/// where (i, j, k) are the sequence axes.
[[nodiscard]] Mat3 eulerToMatrix(EulerSequence seq, const EulerAngles& angles) noexcept;

/// Inverse of eulerToMatrix. Middle angle is in [-pi/2, pi/2] for Tait-Bryan
/// sequences and in [0, pi] for proper Euler sequences. The angles returned give `dcm` back through
/// eulerToMatrix to rounding (about 1e-15) wherever the middle angle is.
///
/// At gimbal lock only the sum or the difference of the first and third angles is defined: the third
/// angle is set to zero and the first carries both. Gimbal lock is a middle angle within rounding
/// (2e-16 rad) of +-pi/2 for a Tait-Bryan sequence, or of 0 or pi for a proper Euler sequence. Just
/// outside it the first and third angles are returned separately, but they are only as well
/// determined as the matrix allows: an error of 1e-16 in `dcm` moves them by 1e-16 divided by the
/// distance from lock.
[[nodiscard]] EulerAngles matrixToEuler(EulerSequence seq, const Mat3& dcm) noexcept;

/// Re-orthonormalize a nearly-orthogonal rotation matrix (via quaternion).
[[nodiscard]] Mat3 orthonormalize(const Mat3& m) noexcept;

} // namespace astro
