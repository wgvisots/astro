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
