/**********************************************************************
 * angle_huge_test.cpp
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

// Sexagesimal text with hundreds of digits: a number, but possibly not an angle a double can hold.
// The parsers return a finite angle or throw. (Found by fuzzing; see fuzz/angle.cpp.)

#include <cmath>
#include <string>

#include <gtest/gtest.h>

#include "astro/core/angle.h"
#include "astro/core/constants.h"
#include "astro/core/error.h"

using namespace astro;

TEST(Angle, HugeSexagesimalValuesAreFiniteOrRefused) {
    // 6.6...e300 hours. Hours were turned into degrees (times 15) before radians, which overflowed
    // to infinity although the angle in radians, 1.7e300, is an ordinary double.
    const std::string hours = std::string(301, '6');
    EXPECT_NEAR(parseHms(hours + "h") / (6.0 / 9.0 * 1e301 * constants::pi / 12.0), 1.0, 1e-14);
    EXPECT_TRUE(std::isfinite(parseDms(hours + "d")));
    // 1.79e308 hours and 9e307 minutes: each is a double, their sum is not.
    const std::string tooMuch = "179" + std::string(306, '0') + ":" + "9" + std::string(307, '0');
    EXPECT_THROW((void)parseHms(tooMuch), Error);
    EXPECT_THROW((void)parseDms(tooMuch), Error);
}
