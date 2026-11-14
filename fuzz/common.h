/**********************************************************************
 * common.h
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

#pragma once

// Shared pieces of the fuzz targets.
//
// A target feeds the bytes it is given to one parser. The contract it checks is the one in
// SECURITY.md: malformed input must end in astro::Error (or in a valid result), never in a crash,
// a hang, a memory error or undefined behaviour. So astro::Error is caught and ignored; anything
// else that escapes, and anything the sanitizers report, stops the fuzzer with the input saved.

#include <cstddef>
#include <cstdint>
#include <string_view>

#include "astro/core/error.h"

namespace astro::fuzz {

/// The input bytes as text.
inline std::string_view asText(const std::uint8_t* data, std::size_t size) {
    return {reinterpret_cast<const char*>(data), size};
}

/// Stops the fuzzer (the input is saved as a crash) if a property that must hold does not.
inline void require(bool condition) {
    if (!condition) {
        __builtin_trap();
    }
}

} // namespace astro::fuzz
