/**********************************************************************
 * orbital_frames.h
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

#include <array>

#include "astro/core/vec.h"

/// Satellite-centered orbital frames. Each function returns the rotation from
/// the inertial frame of `s` into the orbital frame: v_frame = M * v_inertial
/// (the rows of M are the frame's unit vectors).
namespace astro {

/// RIC / RSW: Radial (r), In-track (completes the triad), Cross-track (angular momentum).
[[nodiscard]] Mat3 ricMatrix(const StateVector& s) noexcept;
/// LVLH (CCSDS/NASA convention): z = -r (nadir), y = -h (negative orbit normal), x = y cross z.
[[nodiscard]] Mat3 lvlhMatrix(const StateVector& s) noexcept;
/// VNC: Velocity, Normal (angular momentum), Co-normal (V cross N).
[[nodiscard]] Mat3 vncMatrix(const StateVector& s) noexcept;
/// NTW: in-plane Normal to velocity (N = T cross W), Tangential (velocity), W (angular momentum).
[[nodiscard]] Mat3 ntwMatrix(const StateVector& s) noexcept;

/// Position and velocity of `deputy` relative to `chief`, expressed in the chief's
/// rotating RIC frame (velocity includes the frame rotation, omega = h / r^2).
[[nodiscard]] StateVector relativeStateRic(const StateVector& chief, const StateVector& deputy) noexcept;
/// Inverse of relativeStateRic: inertial deputy state from a relative RIC state.
[[nodiscard]] StateVector deputyFromRelativeRic(const StateVector& chief, const StateVector& relRic) noexcept;

using Mat6 = std::array<std::array<double, 6>, 6>;

/// Rotate a 6x6 position/velocity covariance by a (static) 3x3 frame rotation:
/// P' = diag(R, R) P diag(R, R)^T.
[[nodiscard]] Mat6 rotateCovariance(const Mat6& P, const Mat3& R) noexcept;

} // namespace astro
