/**********************************************************************
 * eop.cpp
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

// Fuzz target: Earth orientation files (finals2000A and C04).

#include <cmath>

#include "astro/frames/eop.h"
#include "common.h"

using namespace astro;

namespace {

void use(const EopTable& table) {
    // A table that parses interpolates to finite values everywhere inside its range.
    const double first = table.firstMjd(), last = table.lastMjd();
    for (double f : {0.0, 0.25, 0.5, 1.0}) {
        const EopValues v = table.at(first + f * (last - first));
        fuzz::require(std::isfinite(v.xp + v.yp + v.dut1 + v.lod + v.dX + v.dY));
    }
}

} // namespace

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size) {
    const std::string_view text = fuzz::asText(data, size);
    for (EopInterpolation interp : {EopInterpolation::Linear, EopInterpolation::Lagrange4}) {
        try {
            use(EopTable::parseFinals2000A(text, interp));
        } catch (const Error&) {
        }
        try {
            use(EopTable::parseC04(text, interp));
        } catch (const Error&) {
        }
    }
    return 0;
}
