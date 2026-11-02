/**********************************************************************
 * state12.h
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

#include <iosfwd>

#include <boost/describe/enum.hpp>

#include "astro/core/vec.h"

/// Options and results shared by the propagators that take a 12-component state
/// (position, velocity, acceleration, jerk).
namespace astro {

/// What a propagator does with the acceleration and jerk of an input StateVector12.
///
/// Under the propagator's force model, position and velocity alone determine the trajectory, and
/// the model also determines acceleration and jerk. The two modes differ in how they treat input
/// acceleration and jerk that disagree with the model.
enum class HigherOrderMode {
    /// Propagate from position and velocity only. The input acceleration and jerk are compared
    /// with the model (see Propagation12::accelerationResidual) but do not affect the result.
    ModelOnly,
    /// Treat the difference between the input and modeled acceleration as an extra, unmodeled
    /// inertial force (thrust, drag, ...): da(t) = da0 + dj0 * t, where da0 and dj0 are the
    /// acceleration and jerk residuals at the epoch. The propagated state then reproduces the
    /// input acceleration and jerk exactly at dt = 0.
    UnmodeledForce,
};

BOOST_DESCRIBE_ENUM(HigherOrderMode, ModelOnly, UnmodeledForce)

/// Writes "HigherOrderMode::<enumerator>", or "HigherOrderMode::INVALID_VALUE" for any other value.
std::ostream& operator<<(std::ostream& os, HigherOrderMode value);

/// Result of propagating a StateVector12.
struct Propagation12 {
    StateVector12 state;                ///< state at epoch + dt; a and j include any unmodeled force
    Vec3          accelerationResidual; ///< input acceleration minus the force model's, at the epoch
    Vec3          jerkResidual;         ///< input jerk minus the force model's, at the epoch
};

} // namespace astro
