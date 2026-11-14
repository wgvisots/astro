/**********************************************************************
 * iso8601.cpp
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

// Fuzz target: ISO-8601 times.

#include <cmath>
#include <string>

#include "astro/time/time.h"
#include "common.h"

using namespace astro;

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size) {
    try {
        const Time t = Time::fromIso8601(fuzz::asText(data, size));
        // Whatever parses can be written, and its text parses back to the same instant (to the
        // nanosecond the text carries).
        const std::string text = t.iso8601(9);
        const Time        back = Time::fromIso8601(text);
        fuzz::require(std::abs(back - t) < 2e-9);
    } catch (const Error&) {
    }
    return 0;
}
