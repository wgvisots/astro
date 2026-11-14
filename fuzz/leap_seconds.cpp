/**********************************************************************
 * leap_seconds.cpp
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

// Fuzz target: leap-second files (tai-utc.dat, Leap_Second.dat, leap-seconds.list).

#include <cmath>

#include "astro/time/leap_seconds.h"
#include "common.h"

using namespace astro;

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size) {
    try {
        const LeapSecondTable table = LeapSecondTable::parse(fuzz::asText(data, size));
        // A table that parses answers for any date from 1960 on with a finite number.
        for (int year : {1960, 1965, 1972, 1999, 2017, 2100, 9999}) {
            fuzz::require(std::isfinite(table.taiMinusUtc(year, 7, 1)));
        }
    } catch (const Error&) {
    }
    return 0;
}
