/**********************************************************************
 * error.h
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

#include <stdexcept>
#include <string>

/// \file
/// astro::Error, the exception thrown for bad input and missing data.
namespace astro {

/// Base exception for all errors raised by the astro library: invalid input,
/// missing data (e.g. no EOP for the requested epoch) or non-convergence.
class Error : public std::runtime_error {
public:
    /// Error with message `what` (returned by what()).
    explicit Error(const std::string& what)
        : std::runtime_error(what) {}
};

} // namespace astro
