/**********************************************************************
 * angle.cpp
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

// Fuzz target: sexagesimal angles ("12:34:56.7", "-12d 34' 56.7\"").

#include <cmath>

#include "astro/core/angle.h"
#include "common.h"

using namespace astro;

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data, std::size_t size) {
    const std::string_view text = fuzz::asText(data, size);
    try {
        fuzz::require(std::isfinite(parseHms(text)));
    } catch (const Error&) {
    }
    try {
        fuzz::require(std::isfinite(parseDms(text)));
    } catch (const Error&) {
    }
    return 0;
}
