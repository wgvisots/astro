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

#include <random>

#include "astro/core/angle.h"
#include "test_util.h"

using namespace astro;
using astro::test::expectMatNear;
using astro::test::expectOrthonormal;
using astro::test::expectVecNear;

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

TEST(Quat, Slerp) {
    const Quat a = Quat::fromAxisAngle({0, 0, 1}, 0.0);
    const Quat b = Quat::fromAxisAngle({0, 0, 1}, 1.0);
    expectMatNear(slerp(a, b, 0.25).toMatrix(), R3(0.25), 1e-14);
    expectMatNear(slerp(a, b, 1.0).toMatrix(), R3(1.0), 1e-14);
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
