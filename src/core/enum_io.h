/**********************************************************************
 * enum_io.h
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

#include <ostream>
#include <string>
#include <string_view>

#include <boost/describe/enum_to_string.hpp>

namespace astro::detail {

/// Writes "type::Enumerator" as one string (so std::setw applies to the whole token), or
/// "type::INVALID_VALUE" if value is not a described enumerator of E.
template <class E>
std::ostream& writeEnum(std::ostream& os, std::string_view type, E value) {
    const char* name = boost::describe::enum_to_string(value, nullptr);
    std::string s(type);
    s.append("::").append(name ? name : "INVALID_VALUE");
    return os << s;
}

/// Bare enumerator name, or fallback if value is not a described enumerator of E.
template <class E>
std::string_view enumName(E value, std::string_view fallback) noexcept {
    const char* name = boost::describe::enum_to_string(value, nullptr);
    return name ? std::string_view(name) : fallback;
}

} // namespace astro::detail
