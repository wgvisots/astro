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
/// astro::Error, the exception thrown for bad input and missing data, and the error policy of the library.
///
/// Every function follows the same rules, so a caller can tell from a declaration how it reports a failure:
///
/// 1. **A failure throws astro::Error**: invalid input, missing data (no Earth orientation data for an
///    epoch, say) or a solution that does not converge. A subclass exists only where it carries data the
///    caller can act on: today only Sgp4Error, with its Sgp4Status.
/// 2. **A message names the function first**, in the form `"Function: what is wrong (got <value>)"`, for
///    example `"gcrsToCelestialMatrix: invalid Frame value (got 99)"`. The text after the name is for
///    people and may change; do not parse it.
/// 3. **A `noexcept` function evaluates a formula and does not check its input.** Not-a-number (NaN) or
///    infinity in gives NaN out, and outside its documented domain it returns whatever the formula gives.
///    A function that checks its input throws, and is not `noexcept`. (FrameTransformer::gcrsTo and
///    FrameTransformer::matrix return NaN for a value that is not a Frame enumerator: they are formulas on
///    a hot path.)
/// 4. **An outcome that is normal but has no value is a `std::optional`**: predictImpact when there is no
///    impact, propagateToEvent when the event does not happen, mahalanobis2 for a covariance that is not
///    positive definite.
/// 5. **A value that was not computed is a documented NaN** in a field of the result, such as
///    LookAngles::rangeRate when no velocity was given.
/// 6. **No out-parameters.** A function with more than one result returns a struct (Stumpff,
///    YearMonthDay, Sgp4Result).
/// 7. **A non-throwing variant exists only where a caller must handle a failure at every sample of a
///    loop** (today only SGP4). It is named `try...`, is `noexcept`, and returns a struct holding the
///    value and a status `enum class` (Sgp4::tryPropagate and Sgp4Result).
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
