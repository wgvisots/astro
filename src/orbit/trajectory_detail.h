/**********************************************************************
 * trajectory_detail.h
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

#include <string_view>

#include "astro/frames/transform.h"
#include "astro/orbit/trajectory.h"

// Frame checks shared by the functions that take a Trajectory. Each throws astro::Error with a
// message that starts with `function` and names the frame it got and the frame it needs.
namespace astro::detail {

/// Throws unless `t` is in `needed`. `what` names the argument ("the satellite trajectory").
void requireFrame(const Trajectory& t, Frame needed, std::string_view function, std::string_view what);

/// Throws unless `a` and `b` are in the same frame. `whatA` and `whatB` name the arguments.
void requireSameFrame(const Trajectory& a,
                      const Trajectory& b,
                      std::string_view  function,
                      std::string_view  whatA,
                      std::string_view  whatB);

/// Throws if `t` is in an Earth-fixed frame (see isEarthFixed). `what` names the argument.
void requireInertial(const Trajectory& t, std::string_view function, std::string_view what);

} // namespace astro::detail
