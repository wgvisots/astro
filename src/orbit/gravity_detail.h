/**********************************************************************
 * gravity_detail.h
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

#include <cstddef>

#include "astro/core/vec.h"

namespace astro::detail {

/// Point mass (if pointMass) plus zonal J[2..size-1] acceleration and, if WithJerk, jerk, in one
/// pass sharing the square root and Legendre recursion. Defined in gravity.cpp.
template <bool WithJerk>
void gravity(const Vec3&   r,
             const Vec3&   v,
             const double* J,
             std::size_t   size,
             double        mu,
             double        re,
             bool          pointMass,
             Vec3*         acc,
             Vec3*         jerk) noexcept;

} // namespace astro::detail
