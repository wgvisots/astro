/**********************************************************************
 * state_covariance.h
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

#include "astro/core/mat6.h"
#include "astro/core/vec.h"

/// \file
/// A state with its covariance, and a state with its state transition matrix.
namespace astro {

/// A state and its 6x6 position/velocity covariance at one time.
struct StateCovariance {
    StateVector state;      ///< position (m) and velocity (m/s)
    Mat6        covariance; ///< covariance of (r, v) in the same frame (m^2, m^2/s, m^2/s^2)
};

/// A state and its state transition matrix Phi(t, epoch) = d(state at t) / d(state at the epoch),
/// rows and columns ordered (r, v): P(t) = Phi P0 Phi^T.
struct StateTransition {
    StateVector state; ///< position (m) and velocity (m/s)
    Mat6        stm;   ///< Phi(t, epoch)
};

} // namespace astro
