/**********************************************************************
 * vec_test.cpp
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

#include "astro/core/vec.h"

#include <cmath>
#include <limits>
#include <sstream>

#include "astro/core/constants.h"
#include "test_util.h"

using namespace astro;
using astro::test::expectVecNear;

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

} // namespace

TEST(Vec3, NormalizedGivesTheUnitVector) {
    // 3-4-12-13 is a Pythagorean quadruple, so the length is exactly 13 and the truth is known in
    // closed form. One rounding error in each quotient: 1.1e-16.
    expectVecNear(normalized(Vec3{3.0, -4.0, 12.0}), {3.0 / 13.0, -4.0 / 13.0, 12.0 / 13.0}, 1.2e-16);
    // The zero vector has no direction: the documented result is the zero vector.
    EXPECT_EQ(normalized(Vec3{}), Vec3{});
    EXPECT_EQ(normalized(Vec3{-0.0, 0.0, -0.0}), Vec3{});
}

// A vector with a not-a-number (NaN) component has no direction either, but it is not the zero vector:
// the NaN must come back out, not be replaced by zeros (which the caller would take for a valid answer).
TEST(Vec3, NormalizedKeepsNaN) {
    for (const Vec3& v : {
             Vec3{kNaN, 1.0,  2.0 },
             Vec3{1.0,  kNaN, 2.0 },
             Vec3{1.0,  2.0,  kNaN},
             Vec3{kNaN, 0.0,  0.0 },
             Vec3{kNaN, kNaN, kNaN}
    }) {
        const Vec3 n = normalized(v);
        EXPECT_TRUE(std::isnan(n.x)) << v << " -> " << n;
        EXPECT_TRUE(std::isnan(n.y)) << v << " -> " << n;
        EXPECT_TRUE(std::isnan(n.z)) << v << " -> " << n;
    }
}

// The tests below use small whole numbers and halves, whose sums, differences and products are exact in
// double precision, so the truth is worked out by hand and compared for equality.

TEST(Vec3, ComponentsByIndex) {
    Vec3 v{1.0, 2.0, 3.0};
    EXPECT_EQ(v[0], 1.0);
    EXPECT_EQ(v[1], 2.0);
    EXPECT_EQ(v[2], 3.0);
    EXPECT_EQ(v[7], 3.0); // any other index is z, as vec.h says
    v[0] = -1.0;
    v[1] = -2.0;
    v[2] = -3.0;
    EXPECT_EQ(v, Vec3(-1.0, -2.0, -3.0));
}

TEST(Vec3, ArithmeticIsComponentWise) {
    const Vec3 a{1.0, -2.0, 3.5}, b{0.5, 4.0, -1.0};
    EXPECT_EQ(a + b, Vec3(1.5, 2.0, 2.5));
    EXPECT_EQ(a - b, Vec3(0.5, -6.0, 4.5));
    EXPECT_EQ(-a, Vec3(-1.0, 2.0, -3.5));
    EXPECT_EQ(a * 2.0, Vec3(2.0, -4.0, 7.0));
    EXPECT_EQ(2.0 * a, Vec3(2.0, -4.0, 7.0));
    EXPECT_EQ(a / 2.0, Vec3(0.5, -1.0, 1.75));
    Vec3 c  = a;
    c      += b;
    EXPECT_EQ(c, a + b);
    c -= b;
    EXPECT_EQ(c, a);
    c *= -4.0;
    EXPECT_EQ(c, Vec3(-4.0, 8.0, -14.0));
    c /= -4.0;
    EXPECT_EQ(c, a);
}

TEST(Vec3, DotCrossAndNorm) {
    // (1, 2, 3) . (4, -5, 6) = 4 - 10 + 18 = 12.
    EXPECT_EQ(dot(Vec3{1.0, 2.0, 3.0}, Vec3{4.0, -5.0, 6.0}), 12.0);
    // (1, 2, 3) x (4, 5, 6) = (2*6 - 3*5, 3*4 - 1*6, 1*5 - 2*4) = (-3, 6, -3).
    EXPECT_EQ(cross(Vec3{1.0, 2.0, 3.0}, Vec3{4.0, 5.0, 6.0}), Vec3(-3.0, 6.0, -3.0));
    // Right-handed: x cross y is z, and so on round the cycle; the other order changes the sign.
    const Vec3 x{1.0, 0.0, 0.0}, y{0.0, 1.0, 0.0}, z{0.0, 0.0, 1.0};
    EXPECT_EQ(cross(x, y), z);
    EXPECT_EQ(cross(y, z), x);
    EXPECT_EQ(cross(z, x), y);
    EXPECT_EQ(cross(y, x), -z);
    // 3-4-12-13: the length is exactly 13.
    EXPECT_EQ(norm2(Vec3{3.0, -4.0, 12.0}), 169.0);
    EXPECT_EQ(norm(Vec3{3.0, -4.0, 12.0}), 13.0);
    EXPECT_EQ(norm(Vec3{}), 0.0);
}

// angleBetween uses atan2 of |a x b| and a . b, which keeps every digit for tiny angles and for angles near
// pi, where acos and asin of a dot product lose half of them. Truth: the angle between (1, 0, 0) and
// (cos t, sin t, 0) is t, here from 1e-12 rad (where acos(cos t) gives 0) to pi - 1e-9 rad. Tolerance:
// a rounding error or two of t (worst seen 0, that is, equal).
TEST(Vec3, AngleBetweenIsAccurateAtEveryAngle) {
    const Vec3 x{1.0, 0.0, 0.0};
    for (double t : {1e-12, 1e-8, 1e-3, 0.5, 1.0, 2.0, 3.0, constants::pi - 1e-6, constants::pi - 1e-9}) {
        const double expected = t;
        EXPECT_NEAR(angleBetween(x, Vec3{std::cos(t), std::sin(t), 0.0}), expected, 4.5e-16 * t) << t;
        // The lengths of the vectors do not matter.
        EXPECT_NEAR(angleBetween(1e3 * x, Vec3{1e-3 * std::cos(t), 1e-3 * std::sin(t), 0.0}), expected, 4.5e-16 * t)
            << t;
    }
    EXPECT_EQ(angleBetween(x, Vec3{0.0, 2.0, 0.0}), constants::halfPi);
    EXPECT_EQ(angleBetween(x, -x), constants::pi);
    EXPECT_EQ(angleBetween(x, x), 0.0);
}

TEST(Vec3, StreamsAsBracketedList) {
    std::ostringstream os;
    os << Vec3{1.0, -2.5, 3.0};
    EXPECT_EQ(os.str(), "[1, -2.5, 3]");
}

TEST(StateVector12, PositionAndVelocityPart) {
    const StateVector12 s{
        {1.0,  2.0,  3.0 },
        {4.0,  5.0,  6.0 },
        {7.0,  8.0,  9.0 },
        {10.0, 11.0, 12.0}
    };
    EXPECT_EQ(s.rv(),
              (StateVector{
                  {1.0, 2.0, 3.0},
                  {4.0, 5.0, 6.0}
    }));
}

// clang-format off
// A matrix with no two elements alike, so that any element read from the wrong place shows:
//   [ 1  2  3 ]
//   [ 4  5  6 ]
//   [ 7  8 10 ]
const Mat3 kM{1.0, 2.0,  3.0,
              4.0, 5.0,  6.0,
              7.0, 8.0, 10.0};
// clang-format on

TEST(Mat3, ElementsRowsAndColumns) {
    EXPECT_EQ(kM(0, 1), 2.0); // row 0, column 1
    EXPECT_EQ(kM(1, 0), 4.0);
    EXPECT_EQ(kM(2, 2), 10.0);
    EXPECT_EQ(kM.row(1), Vec3(4.0, 5.0, 6.0));
    EXPECT_EQ(kM.col(1), Vec3(2.0, 5.0, 8.0));
    EXPECT_EQ(Mat3::fromRows({1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 10.0}), kM);
    EXPECT_EQ(Mat3::fromCols({1.0, 4.0, 7.0}, {2.0, 5.0, 8.0}, {3.0, 6.0, 10.0}), kM);
    EXPECT_EQ(kM.transpose(), Mat3(1.0, 4.0, 7.0, 2.0, 5.0, 8.0, 3.0, 6.0, 10.0));
    EXPECT_EQ(Mat3{}, Mat3(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0));
    EXPECT_EQ(Mat3::identity(), Mat3(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0));
    Mat3 m;
    m(2, 0) = 5.0;
    EXPECT_EQ(m.col(0), Vec3(0.0, 0.0, 5.0));
}

TEST(Mat3, TraceAndDeterminant) {
    EXPECT_EQ(kM.trace(), 16.0); // 1 + 5 + 10
    // Expanded along the first row: 1 (5*10 - 6*8) - 2 (4*10 - 6*7) + 3 (4*8 - 5*7) = 2 + 4 - 9 = -3.
    EXPECT_EQ(kM.determinant(), -3.0);
    EXPECT_EQ(kM.transpose().determinant(), -3.0);
    EXPECT_EQ(Mat3::identity().determinant(), 1.0);
}

TEST(Mat3, ArithmeticAndProducts) {
    const Mat3 i = Mat3::identity();
    EXPECT_EQ(kM + i, Mat3(2.0, 2.0, 3.0, 4.0, 6.0, 6.0, 7.0, 8.0, 11.0));
    EXPECT_EQ(kM - i, Mat3(0.0, 2.0, 3.0, 4.0, 4.0, 6.0, 7.0, 8.0, 9.0));
    EXPECT_EQ(kM * 0.5, Mat3(0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0));
    EXPECT_EQ(0.5 * kM, kM * 0.5);
    Mat3 m  = kM;
    m      += i;
    m      -= 2.0 * i;
    m      *= -1.0;
    EXPECT_EQ(m, Mat3(0.0, -2.0, -3.0, -4.0, -4.0, -6.0, -7.0, -8.0, -9.0));
    // Row times column, worked out by hand: row 0 of kM with column 0 of kM is 1*1 + 2*4 + 3*7 = 30, ...
    EXPECT_EQ(kM * kM, Mat3(30.0, 36.0, 45.0, 66.0, 81.0, 102.0, 109.0, 134.0, 169.0));
    EXPECT_EQ(kM * i, kM);
    EXPECT_EQ(i * kM, kM);
    // (1, 2, 3) . rows: 1 + 4 + 9 = 14, 4 + 10 + 18 = 32, 7 + 16 + 30 = 53.
    EXPECT_EQ(kM * Vec3(1.0, 2.0, 3.0), Vec3(14.0, 32.0, 53.0));
}

TEST(Mat3, OuterProductAndSkewMatrix) {
    // a b^T: element (i, j) is a_i b_j.
    EXPECT_EQ(outer(Vec3{1.0, 2.0, 3.0}, Vec3{4.0, 5.0, 6.0}), Mat3(4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 12.0, 15.0, 18.0));
    // skew(a) b is a x b: (1, 2, 3) x (4, 5, 6) = (-3, 6, -3), as above.
    const Mat3 s = Mat3::skew({1.0, 2.0, 3.0});
    EXPECT_EQ(s * Vec3(4.0, 5.0, 6.0), Vec3(-3.0, 6.0, -3.0));
    EXPECT_EQ(s.transpose(), -1.0 * s); // antisymmetric
    EXPECT_EQ(s, Mat3(0.0, -3.0, 2.0, 3.0, 0.0, -1.0, -2.0, 1.0, 0.0));
}

TEST(Mat3, StreamsRowByRow) {
    std::ostringstream os;
    os << kM;
    EXPECT_EQ(os.str(), "[[1, 2, 3], [4, 5, 6], [7, 8, 10]]");
}
