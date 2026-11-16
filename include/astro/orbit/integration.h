/**********************************************************************
 * integration.h
 *
 * Copyright 2024-2026 W. Greg Visotski
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

/// \file
/// Settings and work statistics of the adaptive integrator behind NumericalPropagator.
namespace astro {

/// Settings for the adaptive Runge-Kutta-Fehlberg 7(8) integrator that NumericalPropagator (and
/// the propagators built on it) uses. Steps and times are in seconds.
///
/// **What the error control cannot see.** Fehlberg's error estimate compares derivative evaluations
/// at the two ends of a step. A part of the acceleration that depends on time alone (a thrust
/// profile a(t), say) cancels in it exactly, whatever its shape, so the step is chosen for the rest
/// of the forces only. For such a force that varies smoothly, set maxStep to a fraction of the time
/// over which it changes; where it jumps (an engine switching on or off), end the propagation at the
/// jump and start a new one there. Forces that depend on the state are not affected.
struct AdaptiveOptions {
    double relTol{1e-12};      ///< relative error tolerance per component
    double absTol{1e-9};       ///< absolute error tolerance per component (m and m/s)
    double initialStep{60.0};  ///< size of the first trial step (clipped to the span and to maxStep)
    double minStep{1e-6};      ///< smallest step; throws astro::Error if the error is still too large
    double maxStep{1e30};      ///< largest step, the first included
    int    maxSteps{10000000}; ///< limit on attempted steps; throws astro::Error when exceeded
    /// Only the first errorComponents components enter the error control (e.g. the state but not
    /// an attached state transition matrix). Default: all.
    std::size_t errorComponents{static_cast<std::size_t>(-1)};
};

/// Work done by an adaptive integration.
struct AdaptiveStats {
    int accepted{0};    ///< accepted steps
    int rejected{0};    ///< rejected (retried) steps
    int evaluations{0}; ///< evaluations of the derivative function (13 per attempted step)
};

} // namespace astro
