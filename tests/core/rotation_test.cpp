/**********************************************************************
 * rotation_test.cpp
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
#include <array>
#include <cmath>
#include <limits>
#include <random>
#include <string>
#include <vector>

#include <boost/describe/enum_to_string.hpp>

#include "astro/core/angle.h"
#include "astro/core/error.h"
#include "test_util.h"

using namespace astro;
using astro::test::expectMatNear;
using astro::test::expectOrthonormal;
using astro::test::expectVecNear;

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kInf = std::numeric_limits<double>::infinity();

void expectAllNaN(const Quat& q) {
    EXPECT_TRUE(std::isnan(q.w)) << "w = " << q.w;
    EXPECT_TRUE(std::isnan(q.x)) << "x = " << q.x;
    EXPECT_TRUE(std::isnan(q.y)) << "y = " << q.y;
    EXPECT_TRUE(std::isnan(q.z)) << "z = " << q.z;
}

void expectAllNaN(const Mat3& m) {
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            EXPECT_TRUE(std::isnan(m(i, j))) << "(" << i << "," << j << ") = " << m(i, j);
        }
    }
}

void expectQuatNear(const Quat& a, const Quat& b, double tol, const char* what = "") {
    EXPECT_NEAR(a.w, b.w, tol) << what << " w";
    EXPECT_NEAR(a.x, b.x, tol) << what << " x";
    EXPECT_NEAR(a.y, b.y, tol) << what << " y";
    EXPECT_NEAR(a.z, b.z, tol) << what << " z";
}

/// Frame rotation by `theta` about the unit vector `n`, from Rodrigues' formula in its passive form:
/// A = cos(theta) I + (1 - cos(theta)) n n^T - sin(theta) [n x]. For n = z this is R3(theta).
Mat3 rodrigues(const Vec3& n, double theta) {
    const double c = std::cos(theta), s = std::sin(theta);
    return c * Mat3::identity() + (1.0 - c) * outer(n, n) - s * Mat3::skew(n);
}

// The Euler tests below build and compare matrices in extended precision (long double), with the
// rotation matrices written out again from their definition, so that the truth owes nothing to the
// library and has no rounding error of its own at the level being tested.
using Mat3L = std::array<std::array<long double, 3>, 3>;

/// Passive rotation about axis 1, 2 or 3 (as R1, R2, R3 in rotation.h) from the cosine and sine of its angle.
Mat3L frameRotation(int axis, long double c, long double s) {
    // clang-format off
    switch (axis) {
    case 1:  return {{{1, 0,  0}, { 0, c, s}, {0, -s, c}}};
    case 2:  return {{{c, 0, -s}, { 0, 1, 0}, {s,  0, c}}};
    default: return {{{c, s,  0}, {-s, c, 0}, {0,  0, 1}}};
    }
    // clang-format on
}

Mat3L product(const Mat3L& a, const Mat3L& b) {
    Mat3L p{};
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            p[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j];
        }
    }
    return p;
}

/// Direction cosine matrix of an Euler sequence, R_k(c) R_j(b) R_i(a). The middle rotation is given by
/// its cosine and sine, so that gimbal lock can be set exactly.
Mat3L eulerMatrix(const std::array<int, 3>& axes, long double a, long double cosB, long double sinB, long double c) {
    return product(product(frameRotation(axes[2], std::cos(c), std::sin(c)), frameRotation(axes[1], cosB, sinB)),
                   frameRotation(axes[0], std::cos(a), std::sin(a)));
}

/// Axes of an Euler sequence read from its name (1 = x, 2 = y, 3 = z): "ZYX" is {3, 2, 1}.
std::array<int, 3> axesFromName(EulerSequence seq) {
    const std::string name = boost::describe::enum_to_string(seq, "???");
    return {name[0] - 'X' + 1, name[1] - 'X' + 1, name[2] - 'X' + 1};
}

/// The matrix rounded to double precision: the nearest a caller could hold.
Mat3 rounded(const Mat3L& m) {
    Mat3 r;
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            r(i, j) = static_cast<double>(m[i][j]);
        }
    }
    return r;
}

/// Largest difference between an element of `a` and the same element of `b`.
double worstDifference(const Mat3L& a, const Mat3& b) {
    long double worst = 0.0L;
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            worst = std::max(worst, std::abs(a[i][j] - static_cast<long double>(b(i, j))));
        }
    }
    return static_cast<double>(worst);
}

} // namespace

TEST(Rotation, PassiveConvention) {
    // Rotating the frame +90 deg about z: the old x axis appears at -y in the new frame.
    const double t = deg2rad(90.0);
    expectVecNear(R3(t) * Vec3{1, 0, 0}, {0, -1, 0}, 1e-15);
    expectVecNear(R1(t) * Vec3{0, 1, 0}, {0, 0, -1}, 1e-15);
    // R2: old z axis appears at -x; old x axis appears at +z.
    expectVecNear(R2(t) * Vec3{0, 0, 1}, {-1, 0, 0}, 1e-15);
    expectVecNear(R2(t) * Vec3{1, 0, 0}, {0, 0, 1}, 1e-15);
}

// All nine elements of R1, R2 and R3, at angles in every quadrant, against the passive rotation matrices
// of SOFA's iauRx, iauRy and iauRz (the convention of the Explanatory Supplement), written out again above
// (frameRotation) in extended precision. The only rounding is in the double cosine and sine of the
// library: half a unit in the last place, so 1.1e-16 at most (measured 5.5e-17).
TEST(Rotation, ElementaryRotationsMatchTheirDefinitionInEveryElement) {
    for (double t : {0.0, 0.3, 1.2, 2.0, 3.0, -0.7, -2.5, 10.0}) {
        const long double tL = t;
        const long double c = std::cos(tL), s = std::sin(tL);
        SCOPED_TRACE(testing::Message() << "theta " << t);
        EXPECT_LE(worstDifference(frameRotation(1, c, s), R1(t)), 1.2e-16);
        EXPECT_LE(worstDifference(frameRotation(2, c, s), R2(t)), 1.2e-16);
        EXPECT_LE(worstDifference(frameRotation(3, c, s), R3(t)), 1.2e-16);
        EXPECT_LE(worstDifference(frameRotation(1, c, s), rotationAboutAxis(1, t)), 1.2e-16);
        EXPECT_LE(worstDifference(frameRotation(2, c, s), rotationAboutAxis(2, t)), 1.2e-16);
        EXPECT_LE(worstDifference(frameRotation(3, c, s), rotationAboutAxis(3, t)), 1.2e-16);
    }
}

TEST(Rotation, CyclicConsistency) {
    // The three elementary rotations are related by cyclic permutation of axes.
    // Permutation P maps (x,y,z) -> (y,z,x); then R2(t) = P R1(t) P^T, R3(t) = P R2(t) P^T.
    const Mat3 P{0, 1, 0, 0, 0, 1, 1, 0, 0};
    for (double t : {0.1, 1.2, -2.5}) {
        expectMatNear(P.transpose() * R1(t) * P, R2(t), 1e-15, "R2");
        expectMatNear(P.transpose() * R2(t) * P, R3(t), 1e-15, "R3");
        expectOrthonormal(R1(t));
        expectOrthonormal(R2(t));
        expectOrthonormal(R3(t));
    }
}

TEST(Quat, MatchesElementaryRotations) {
    for (double t : {0.0, 0.4, -1.3, 3.0}) {
        expectMatNear(Quat::fromAxisAngle({1, 0, 0}, t).toMatrix(), R1(t), 1e-15);
        expectMatNear(Quat::fromAxisAngle({0, 1, 0}, t).toMatrix(), R2(t), 1e-15);
        expectMatNear(Quat::fromAxisAngle({0, 0, 1}, t).toMatrix(), R3(t), 1e-15);
    }
}

TEST(Quat, CompositionOrder) {
    const Quat p = Quat::fromAxisAngle({0, 0, 1}, 0.7);
    const Quat q = Quat::fromAxisAngle({1, 0, 0}, -0.2);
    expectMatNear((p * q).toMatrix(), q.toMatrix() * p.toMatrix(), 1e-15);
}

TEST(Quat, MatrixRoundTripAndTransform) {
    std::mt19937                     rng(42);
    std::normal_distribution<double> n(0.0, 1.0);
    for (int i = 0; i < 500; ++i) {
        const Quat q = Quat(n(rng), n(rng), n(rng), n(rng)).normalized();
        const Mat3 A = q.toMatrix();
        // Tolerances: a few rounding errors (worst seen 1.3e-15 for the matrix, 3.3e-16 for the quaternion and
        // 8.9e-16 for the vectors, whose components are of order 1). This test once allowed 1e-14 and 1e-13.
        expectOrthonormal(A, 3e-15);
        Quat       r  = Quat::fromMatrix(A);
        const Quat qq = q.w < 0 ? Quat{-q.w, -q.x, -q.y, -q.z} : q;
        EXPECT_NEAR(r.w, qq.w, 1e-15);
        EXPECT_NEAR(r.x, qq.x, 1e-15);
        EXPECT_NEAR(r.y, qq.y, 1e-15);
        EXPECT_NEAR(r.z, qq.z, 1e-15);
        const Vec3 v{n(rng), n(rng), n(rng)};
        expectVecNear(q.transform(v), A * v, 3e-15);
        expectVecNear(q.rotate(v), A.transpose() * v, 3e-15);
    }
}

// The axis may have any length: only its direction counts. The lengths here run from the smallest
// double (2^-1074) to near the largest. The axis is (3, -4, 12) scaled by a power of two, so its unit
// vector is (3, -4, 12) / 13 exactly. Truth: the quaternion (cos(theta/2), sin(theta/2) n), and
// Rodrigues' formula for the matrix. Before the fix the components were squared as given, which
// underflows below about 1e-162 and overflows above about 1e154, and the result was
// (cos(theta/2), 0, 0, 0): not a rotation.
TEST(Quat, AxisOfAnyLengthIsNormalized) {
    const Vec3 n{3.0 / 13.0, -4.0 / 13.0, 12.0 / 13.0};
    for (int power : {-1074, -1000, -600, -300, -1, 0, 300, 600, 1019}) {
        const Vec3 axis{std::ldexp(3.0, power), std::ldexp(-4.0, power), std::ldexp(12.0, power)};
        for (double theta : {1.0, -2.5, 3.0}) {
            SCOPED_TRACE(testing::Message() << "axis scaled by 2^" << power << ", theta " << theta);
            const double c = std::cos(0.5 * theta), s = std::sin(0.5 * theta);
            // Tolerances. Quaternion: two rounding errors of 1.1e-16, one in the unit vector and one in
            // the product. Matrix: each element is a sum of four squares or products of those
            // components, and the reference has its own rounding; 5.6e-16 was the worst seen over all
            // angles, 1e-15 allows for a different compiler.
            expectQuatNear(Quat::fromAxisAngle(axis, theta), Quat{c, s * n.x, s * n.y, s * n.z}, 2.3e-16);
            expectMatNear(axisAngleToMatrix(axis, theta), rodrigues(n, theta), 1e-15);
        }
    }
    // The shortest axis there is: the smallest double, along -z.
    expectMatNear(axisAngleToMatrix({0.0, 0.0, -5e-324}, 1.0), R3(-1.0), 1e-15, "smallest axis");
}

// A zero axis defines no rotation. Both functions are noexcept, so the documented result is the
// identity. (Before the fix: cos(theta/2)^2 times the identity, which is the zero matrix for theta = pi.)
TEST(Quat, ZeroAxisGivesTheIdentity) {
    for (double theta : {0.0, 1.0, constants::pi, -7.0}) {
        EXPECT_EQ(Quat::fromAxisAngle({0.0, 0.0, 0.0}, theta), Quat{}) << theta;
        EXPECT_EQ(Quat::fromAxisAngle({-0.0, 0.0, -0.0}, theta), Quat{}) << theta;
        EXPECT_EQ(axisAngleToMatrix({0.0, 0.0, 0.0}, theta), Mat3::identity()) << theta;
    }
}

// Not-a-number (NaN) or infinity in the axis or in the angle gives NaN, never a finite matrix.
TEST(Quat, AxisAngleKeepsNaN) {
    for (const Vec3& axis : {
             Vec3{kNaN, 0.0,   0.0 },
             Vec3{0.0,  kNaN,  0.0 },
             Vec3{0.0,  1.0,   kNaN},
             Vec3{kNaN, kNaN,  kNaN},
             Vec3{0.0,  0.0,   kInf},
             Vec3{1.0,  -kInf, 2.0 }
    }) {
        const Quat q = Quat::fromAxisAngle(axis, 1.0);
        EXPECT_TRUE(std::isnan(q.x) && std::isnan(q.y) && std::isnan(q.z)) << axis << " -> " << q.vec();
        expectAllNaN(axisAngleToMatrix(axis, 1.0));
    }
    for (double theta : {kNaN, kInf, -kInf}) {
        expectAllNaN(Quat::fromAxisAngle({0.0, 0.0, 1.0}, theta));
        expectAllNaN(Quat::fromAxisAngle({0.0, 0.0, 0.0}, theta));
        expectAllNaN(axisAngleToMatrix({0.0, 0.0, 0.0}, theta));
    }
}

// A quaternion with a not-a-number (NaN) component is not a rotation, and it is not the zero quaternion:
// normalising it must give NaN, not the identity (which the caller would take for a valid attitude).
TEST(Quat, NormalizedKeepsNaN) {
    for (const Quat& q : {
             Quat{kNaN, 0.0,  0.0,  1.0 },
             Quat{1.0,  kNaN, 0.0,  0.0 },
             Quat{1.0,  0.0,  0.0,  kNaN},
             Quat{kNaN, kNaN, kNaN, kNaN}
    }) {
        expectAllNaN(q.normalized());
    }
    // The zero quaternion keeps its documented result, the identity.
    EXPECT_EQ(Quat(0.0, 0.0, 0.0, 0.0).normalized(), Quat{});
    EXPECT_EQ(Quat(0.0, 3.0, 0.0, -4.0).normalized(), Quat(0.0, 0.6, 0.0, -0.8));
}

// The functions that finish by normalising a quaternion must pass NaN on in the same way.
TEST(Quat, NaNInGivesNaNOut) {
    Mat3 offDiagonal    = Mat3::identity();
    offDiagonal(0, 1)   = kNaN;
    Mat3 onDiagonal     = Mat3::identity();
    onDiagonal(1, 1)    = kNaN;
    const Quat identity = Quat{};
    for (const Mat3& m : {offDiagonal, onDiagonal}) {
        expectAllNaN(Quat::fromMatrix(m));
        expectAllNaN(orthonormalize(m));
    }
    expectAllNaN(slerp(identity, Quat{kNaN, 0.0, 0.0, 1.0}, 0.5));
    expectAllNaN(slerp(Quat{1.0, 0.0, kNaN, 0.0}, identity, 0.5));
    expectAllNaN(slerp(identity, Quat{0.0, 0.0, 0.0, 1.0}, kNaN));
}

// The edges of the interpolation, between two unrelated attitudes (different axes, angles in all
// quadrants). At t = 0 and t = 1 the result is the end point itself, b taken on a's side; the end points'
// components are exact doubles, so the only error is the rounding of the weights and the final
// normalisation: none was seen; 2.3e-16 allows one rounding error. Interpolating between q and -q, which are the same rotation, stays
// at that rotation for every t.
TEST(Quat, SlerpEndPointsAndTheSameRotationWrittenTwice) {
    const Quat a = Quat(0.5, -0.5, 0.5, 0.5);                       // unit: four halves
    const Quat b = Quat(3.0 / 13.0, -4.0 / 13.0, 12.0 / 13.0, 0.0); // unit: 3-4-12-13
    const Quat c = Quat(-0.6, 0.0, 0.0, -0.8);                      // unit: 3-4-5, on the far side of a
    expectQuatNear(slerp(a, b, 0.0), a, 2.3e-16, "a to b, t = 0:");
    expectQuatNear(slerp(a, b, 1.0), b, 2.3e-16, "a to b, t = 1:");
    // a.c = -0.3 - 0.4 < 0, so the shorter arc ends at -c.
    expectQuatNear(slerp(a, c, 0.0), a, 2.3e-16, "a to c, t = 0:");
    expectQuatNear(slerp(a, c, 1.0), Quat(0.6, 0.0, 0.0, 0.8), 2.3e-16, "a to c, t = 1:");
    for (double t : {0.0, 0.3, 0.5, 1.0}) {
        expectQuatNear(slerp(b, Quat(-b.w, -b.x, -b.y, -b.z), t), b, 2.3e-16, "b to -b:");
    }
}

// Truth: a and b are frame rotations about one axis n by the angles alpha and alpha + sep, so the
// interpolant is the rotation about n by alpha + t sep, known in closed form as
// (cos(angle/2), sin(angle/2) n). Multiplying everything on the left by a fixed quaternion g moves the
// arc to a general place without changing the answer. Separations run from 0 to 3 rad; before the fix
// those up to 0.063 rad were interpolated linearly, which is wrong by up to 1e-6 (1e-7 in a component
// at sep = 0.0632, t = 0.25).
TEST(Quat, SlerpIsAccurateForEveryAngle) {
    const Vec3 n{3.0 / 13.0, -4.0 / 13.0, 12.0 / 13.0}; // unit: 3-4-12-13
    const auto about = [&](double angle) {
        const double c = std::cos(0.5 * angle), s = std::sin(0.5 * angle);
        return Quat{c, s * n.x, s * n.y, s * n.z};
    };
    const Quat   g     = Quat(1.0, 2.0, -3.0, 4.0).normalized();
    const double alpha = 0.7;
    for (double sep :
         {0.0, 1e-300, 1e-12, 1e-9, 2e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 0.02, 0.0632, 0.07, 0.3, 1.0, 3.0}) {
        const Quat a = about(alpha), b = about(alpha + sep);
        for (double t : {0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0}) {
            SCOPED_TRACE(testing::Message() << "sep " << sep << ", t " << t);
            const Quat expected = about(alpha + t * sep);
            // Tolerance: both sides are double precision and each carries a few rounding errors of
            // 1.1e-16 (the reference in its angle, sine and cosine; slerp in its end points, weights and
            // sum). The largest difference over 2 million random (alpha, sep, t) was 5.6e-16.
            expectQuatNear(slerp(a, b, t), expected, 1e-15);
            expectQuatNear(slerp(g * a, g * b, t), g * expected, 1e-15, "moved by g:");
        }
    }
}

// b and -b are the same rotation, so the interpolation takes the shorter arc. With the end points more
// than pi apart as rotations (their quaternions more than pi/2 apart) that arc runs the other way round:
// by sep - 2 pi instead of sep. Same closed-form truth as above.
TEST(Quat, SlerpTakesTheShorterArc) {
    const auto aboutZ = [](double angle) { return Quat{std::cos(0.5 * angle), 0.0, 0.0, std::sin(0.5 * angle)}; };
    for (double sep : {3.2, 4.0, 6.0}) {
        for (double t : {0.0, 0.25, 0.5, 1.0}) {
            SCOPED_TRACE(testing::Message() << "sep " << sep << ", t " << t);
            expectQuatNear(slerp(aboutZ(0.0), aboutZ(sep), t), aboutZ(t * (sep - constants::twoPi)), 1e-15);
        }
    }
}

// rotation.h: "Throws astro::Error for any other `axis`" than 1, 2 or 3.
TEST(Rotation, RotationAboutAxisRejectsOtherAxes) {
    for (int axis : {0, 4, -1, 12}) {
        try {
            (void)rotationAboutAxis(axis, 0.5);
            ADD_FAILURE() << "axis " << axis << ": no error";
        } catch (const Error& e) {
            EXPECT_NE(std::string(e.what()).find("axis"), std::string::npos) << e.what();
        }
    }
}

// Whole-number components, so the truth is exact: |q|^2 = 1 + 4 + 4 + 16 = 25, and q times its conjugate is
// |q|^2 with no vector part (an identity of the Hamilton product that the library does not use).
TEST(Quat, ConjugateNormAndVectorPart) {
    const Quat q{1.0, 2.0, -2.0, 4.0};
    EXPECT_EQ(q.norm(), 5.0);
    EXPECT_EQ(q.vec(), Vec3(2.0, -2.0, 4.0));
    EXPECT_EQ(q.conjugate(), Quat(1.0, -2.0, 2.0, -4.0));
    EXPECT_EQ(q * q.conjugate(), Quat(25.0, 0.0, 0.0, 0.0));
    EXPECT_EQ(q.conjugate() * q, Quat(25.0, 0.0, 0.0, 0.0));
    // For a unit quaternion the conjugate is the inverse rotation: its matrix is the transpose, and it undoes
    // transform(). Tolerance: a few rounding errors of 1.1e-16 in each element (worst seen 2.2e-16).
    const Quat u = q.normalized();
    expectMatNear(u.conjugate().toMatrix(), u.toMatrix().transpose(), 5e-16);
    const Vec3 v{0.3, -1.2, 2.0};
    expectVecNear(u.conjugate().transform(u.transform(v)), v, 1e-15);
    expectVecNear(u.rotate(v), u.conjugate().transform(v), 1e-15);
}

// A half turn about the unit axis n is the matrix 2 n n^T - I (Rodrigues' formula with cos = -1 and sin = 0):
// the trace is -1 and the scalar part of the quaternion is zero, the case Shepperd's method exists for. The
// quaternion is (0, n) or (0, -n), the same rotation. Tolerance: the matrix carries a rounding error or two of
// 1.1e-16 per element, and the square root halves its relative size (worst seen 2.2e-16).
TEST(Quat, FromMatrixOfAHalfTurn) {
    for (const Vec3& n : {
             Vec3{1.0,        0.0,         0.0        },
             Vec3{0.0,        1.0,         0.0        },
             Vec3{0.0,        0.0,         1.0        },
             Vec3{3.0 / 13.0, -4.0 / 13.0, 12.0 / 13.0},
             Vec3{-0.6,       0.0,         0.8        },
             Vec3{2.0 / 3.0,  2.0 / 3.0,   -1.0 / 3.0 }
    }) {
        const Mat3 a = 2.0 * outer(n, n) - Mat3::identity();
        const Quat q = Quat::fromMatrix(a);
        SCOPED_TRACE(testing::Message() << "axis " << n);
        EXPECT_NEAR(q.w, 0.0, 5e-16);
        const double sign = dot(q.vec(), n) < 0.0 ? -1.0 : 1.0;
        expectVecNear(sign * q.vec(), n, 5e-16);
        expectMatNear(q.toMatrix(), a, 5e-16);
    }
}

// Three properties that fix a rotation, none of which the library uses: the axis is left where it is
// (A n = n), the trace is 1 + 2 cos(theta), and a vector v at right angles to n turns by theta in the passive
// sense, n . (v x A v) = -sin(theta). Axes of any length, angles in all quadrants. Tolerance: a few rounding
// errors in each element and in the sums (worst seen 1.0e-15, in the determinant).
TEST(Quat, AxisAngleMatrixFixesTheAxisAndTurnsByTheAngle) {
    std::mt19937                           rng(3);
    std::normal_distribution<double>       g(0.0, 1.0);
    std::uniform_real_distribution<double> angle(-constants::pi, constants::pi);
    for (int k = 0; k < 200; ++k) {
        const Vec3   axis{g(rng), g(rng), g(rng)};
        const double theta = angle(rng);
        const Vec3   n     = normalized(axis);
        const Vec3   v     = normalized(cross(n, Vec3{g(rng), g(rng), g(rng)}));
        const Mat3   a     = axisAngleToMatrix(1e3 * axis, theta);
        SCOPED_TRACE(testing::Message() << "axis " << axis << ", theta " << theta);
        expectVecNear(a * n, n, 2e-15);
        EXPECT_NEAR(a.trace(), 1.0 + 2.0 * std::cos(theta), 2e-15);
        EXPECT_NEAR(dot(n, cross(v, a * v)), -std::sin(theta), 2e-15);
        expectOrthonormal(a, 2e-15);
    }
}

// orthonormalize() takes a matrix that has drifted from a rotation and returns a rotation. Truth: the
// rotation the matrix drifted from, R, and the properties of a rotation (A A^T = I, det A = +1, to a few
// rounding errors: worst seen 1.3e-15). A drift of size e in each element moves the result by up to a few e
// (worst seen 2.7 e over these cases); an exact rotation comes back to rounding (worst seen 4.4e-16).
TEST(Quat, OrthonormalizeRestoresARotation) {
    std::mt19937                           rng(5);
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    for (int k = 0; k < 100; ++k) {
        const Mat3 r = Quat(u(rng), u(rng), u(rng), u(rng)).normalized().toMatrix();
        for (double e : {0.0, 1e-12, 1e-9, 1e-6}) {
            const Mat3 drifted = r + e * Mat3(u(rng), u(rng), u(rng), u(rng), u(rng), u(rng), u(rng), u(rng), u(rng));
            const Mat3 a       = orthonormalize(drifted);
            expectOrthonormal(a, 3e-15);
            expectMatNear(a, r, 4.0 * e + 1e-15, "drift");
        }
    }
}

class EulerTest : public ::testing::TestWithParam<EulerSequence> {};

// The axes come from the name of the sequence ("ZYX" is z, then y, then x) and the matrix from the three
// frame rotations written out in extended precision (eulerMatrix above), so nothing here is taken from the
// library's own eulerAxes or rotation matrices. Any angles are allowed, the middle one included.
// Tolerance: each library element is a sum of products of three double cosines and sines, each within
// 1.1e-16 of its true value, plus the rounding of the products and sums; the worst seen over 15000
// random angles in each sequence was 2.7e-16.
TEST_P(EulerTest, MatrixIsTheProductOfTheThreeFrameRotations) {
    const EulerSequence seq = GetParam();
    const auto          ax  = axesFromName(seq);
    EXPECT_EQ(eulerAxes(seq), ax);
    std::mt19937                           rng(7);
    std::uniform_real_distribution<double> u(-constants::pi, constants::pi);
    for (int i = 0; i < 300; ++i) {
        const double      a = u(rng), b = u(rng), c = u(rng);
        const long double bL = b;
        SCOPED_TRACE(testing::Message() << "a = " << a << ", b = " << b << ", c = " << c);
        EXPECT_LE(worstDifference(eulerMatrix(ax, a, std::cos(bL), std::sin(bL), c), eulerToMatrix(seq, {a, b, c})),
                  6e-16);
    }
}

// What a caller does: angles to a matrix with eulerToMatrix and back with matrixToEuler. Truth: the angles
// put in, with the middle angle anywhere in its range at least 0.01 rad from gimbal lock. The worst
// difference seen over 15000 random angles in each sequence was 4.4e-16 (two rounding errors of an angle
// near pi); 1e-15 allows for a different compiler. (This test once allowed 1e-9.)
TEST_P(EulerTest, RoundTripThroughTheMatrixRecoversTheAngles) {
    const EulerSequence                    seq    = GetParam();
    const auto                             ax     = axesFromName(seq);
    const bool                             proper = ax[0] == ax[2];
    std::mt19937                           rng(7);
    std::uniform_real_distribution<double> u(-3.1, 3.1);
    std::uniform_real_distribution<double> mid(proper ? 0.01 : -1.56, proper ? 3.13 : 1.56);
    for (int i = 0; i < 300; ++i) {
        const EulerAngles e{u(rng), mid(rng), u(rng)};
        const EulerAngles r = matrixToEuler(seq, eulerToMatrix(seq, e));
        SCOPED_TRACE(testing::Message() << "a = " << e.a << ", b = " << e.b << ", c = " << e.c);
        EXPECT_NEAR(r.a, e.a, 1e-15);
        EXPECT_NEAR(r.b, e.b, 1e-15);
        EXPECT_NEAR(r.c, e.c, 1e-15);
    }
    // At gimbal lock, set through the library's constants, the angles returned still give the matrix back
    // to rounding (worst seen 2.4e-16; this once allowed 1e-12).
    for (double b : proper ? std::vector<double>{0.0, constants::pi}
                           : std::vector<double>{constants::halfPi, -constants::halfPi}) {
        const Mat3 A = eulerToMatrix(seq, {0.4, b, -0.3});
        expectMatNear(eulerToMatrix(seq, matrixToEuler(seq, A)), A, 1e-15, "gimbal");
    }
}

// The angles returned for a matrix must give that matrix back, to rounding, wherever the middle angle is:
// far from gimbal lock, at it, and in the band next to it. The middle angle is set a distance d from
// lock, from 1e-18 rad to 1.5 rad, with the other two angles random. Truth: the matrix itself, rebuilt
// from the returned angles in extended precision. Before the fix the rebuilt matrix was wrong by up to
// 2.7e-6 for d below 1.4e-6, and the middle angle by up to 1e-8.
TEST_P(EulerTest, RecomposesToRoundingUpToGimbalLock) {
    const EulerSequence                    seq    = GetParam();
    const auto                             ax     = eulerAxes(seq);
    const bool                             proper = ax[0] == ax[2];
    const long double                      piL    = std::acos(-1.0L);
    std::mt19937                           rng(11);
    std::uniform_real_distribution<double> u(-3.1, 3.1);
    // Gimbal lock is at a middle angle of 0 and pi for the proper Euler sequences, +-pi/2 for the others.
    for (long double lock :
         proper ? std::array<long double, 2>{0.0L, piL} : std::array<long double, 2>{piL / 2, -piL / 2}) {
        const long double inward = lock > 1.0L ? -1.0L : 1.0L; // steps into the range of the middle angle
        for (double d : {1e-18, 1e-17,  1e-16, 2e-16, 2.5e-16, 1e-15, 1e-14, 1e-12, 1e-10, 1e-8, 1e-7,
                         1e-6,  1.4e-6, 3e-6,  1e-5,  1e-4,    1e-3,  1e-2,  0.1,   0.7,   1.5}) {
            for (int n = 0; n < 40; ++n) {
                const double      a = u(rng), c = u(rng);
                const long double b = lock + inward * d;
                SCOPED_TRACE(testing::Message() << "lock at " << static_cast<double>(lock) << ", d = " << d
                                                << ", a = " << a << ", c = " << c);
                const Mat3        A = rounded(eulerMatrix(ax, a, std::cos(b), std::sin(b), c));
                const EulerAngles r = matrixToEuler(seq, A);
                // Tolerance 1e-15. A is within 1.1e-16 of a rotation matrix (it was rounded), the
                // returned angles carry a rounding error each, and inside the lock band (d below
                // 2.2e-16) setting the third angle to zero moves an element by up to 2 d = 4.4e-16.
                // The worst seen over 4.8 million random cases, all sequences, was 5.0e-16.
                const long double bBack = r.b;
                EXPECT_LE(worstDifference(eulerMatrix(ax, r.a, std::cos(bBack), std::sin(bBack), r.c), A), 1e-15);
                // The middle angle is well determined everywhere. 1e-15 is about two rounding errors of
                // pi (4.4e-16 each); the worst seen was one.
                EXPECT_NEAR(r.b, static_cast<double>(b), 1e-15);
                // Outside the lock band the first and third angles are recovered too, to the same
                // 1e-15 (worst seen 4.4e-16). They come from elements of size d, which the rounding of A
                // keeps to 16 digits here; in a matrix with absolute errors of 1e-16 they are only
                // good to 1e-16 / d, which is gimbal lock itself.
                if (d >= 1e-15) {
                    EXPECT_NEAR(r.a, a, 1e-15);
                    EXPECT_NEAR(r.c, c, 1e-15);
                }
            }
        }
    }
}

// At gimbal lock only the sum or the difference of the first and third angles is defined. The documented
// convention is that the third angle is returned as zero (the first then carries both). "At lock" means
// exactly, and also within rounding of it: the constants halfPi and pi are not exactly pi/2 and pi, so
// cos(halfPi) and sin(pi) are about 1e-16 instead of zero.
TEST_P(EulerTest, ThirdAngleIsZeroAtGimbalLock) {
    const EulerSequence                    seq    = GetParam();
    const auto                             ax     = eulerAxes(seq);
    const bool                             proper = ax[0] == ax[2];
    std::mt19937                           rng(13);
    std::uniform_real_distribution<double> u(-3.1, 3.1);
    for (int n = 0; n < 40; ++n) {
        const double a = u(rng), c = u(rng);
        for (long double sign : {1.0L, -1.0L}) {
            // Exact lock: the cosine (or sine) of the middle angle is exactly zero.
            const Mat3L       exact = proper ? eulerMatrix(ax, a, sign, 0.0L, c) : eulerMatrix(ax, a, 0.0L, sign, c);
            const EulerAngles r     = matrixToEuler(seq, rounded(exact));
            EXPECT_EQ(r.c, 0.0);
            EXPECT_LE(worstDifference(exact, eulerToMatrix(seq, r)), 1e-15);
        }
        // The same through the library's own constants.
        for (double b : proper ? std::array<double, 2>{0.0, constants::pi}
                               : std::array<double, 2>{constants::halfPi, -constants::halfPi}) {
            EXPECT_EQ(matrixToEuler(seq, eulerToMatrix(seq, {a, b, c})).c, 0.0) << "middle angle " << b;
        }
    }
}

INSTANTIATE_TEST_SUITE_P(AllSequences,
                         EulerTest,
                         ::testing::Values(EulerSequence::XYZ,
                                           EulerSequence::XZY,
                                           EulerSequence::YXZ,
                                           EulerSequence::YZX,
                                           EulerSequence::ZXY,
                                           EulerSequence::ZYX,
                                           EulerSequence::XYX,
                                           EulerSequence::XZX,
                                           EulerSequence::YXY,
                                           EulerSequence::YZY,
                                           EulerSequence::ZXZ,
                                           EulerSequence::ZYZ));
