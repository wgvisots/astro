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

#include <cmath>
#include <limits>
#include <random>

#include "astro/core/angle.h"
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

TEST(Rotation, R2MatchesExplanatorySupplement) {
    // R2(theta) = [cos 0 -sin; 0 1 0; sin 0 cos]  (the original lib/Astro.cpp had this transposed)
    const double t = 0.3;
    const Mat3   m = R2(t);
    EXPECT_DOUBLE_EQ(m(0, 2), -std::sin(t));
    EXPECT_DOUBLE_EQ(m(2, 0), std::sin(t));
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
        expectOrthonormal(A, 1e-14);
        Quat       r  = Quat::fromMatrix(A);
        const Quat qq = q.w < 0 ? Quat{-q.w, -q.x, -q.y, -q.z} : q;
        EXPECT_NEAR(r.w, qq.w, 1e-13);
        EXPECT_NEAR(r.x, qq.x, 1e-13);
        EXPECT_NEAR(r.y, qq.y, 1e-13);
        EXPECT_NEAR(r.z, qq.z, 1e-13);
        const Vec3 v{n(rng), n(rng), n(rng)};
        expectVecNear(q.transform(v), A * v, 1e-13);
        expectVecNear(q.rotate(v), A.transpose() * v, 1e-13);
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

TEST(Quat, Slerp) {
    const Quat a = Quat::fromAxisAngle({0, 0, 1}, 0.0);
    const Quat b = Quat::fromAxisAngle({0, 0, 1}, 1.0);
    expectMatNear(slerp(a, b, 0.25).toMatrix(), R3(0.25), 1e-14);
    expectMatNear(slerp(a, b, 1.0).toMatrix(), R3(1.0), 1e-14);
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

class EulerTest : public ::testing::TestWithParam<EulerSequence> {};

TEST_P(EulerTest, DefinitionAndRoundTrip) {
    const EulerSequence                    seq    = GetParam();
    const auto                             ax     = eulerAxes(seq);
    const bool                             proper = ax[0] == ax[2];
    std::mt19937                           rng(7);
    std::uniform_real_distribution<double> u(-3.1, 3.1);
    std::uniform_real_distribution<double> mid(proper ? 0.01 : -1.55, proper ? 3.13 : 1.55);
    for (int i = 0; i < 300; ++i) {
        const EulerAngles e{u(rng), mid(rng), u(rng)};
        const Mat3        A = eulerToMatrix(seq, e);
        expectMatNear(A,
                      rotationAboutAxis(ax[2], e.c) * rotationAboutAxis(ax[1], e.b) * rotationAboutAxis(ax[0], e.a),
                      1e-15);
        const EulerAngles r = matrixToEuler(seq, A);
        EXPECT_NEAR(r.a, e.a, 1e-9);
        EXPECT_NEAR(r.b, e.b, 1e-9);
        EXPECT_NEAR(r.c, e.c, 1e-9);
    }
    // Gimbal lock still reproduces the matrix.
    for (double b : proper ? std::vector<double>{0.0, constants::pi}
                           : std::vector<double>{constants::halfPi, -constants::halfPi}) {
        const EulerAngles e{0.4, b, -0.3};
        const Mat3        A = eulerToMatrix(seq, e);
        expectMatNear(eulerToMatrix(seq, matrixToEuler(seq, A)), A, 1e-12, "gimbal");
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
