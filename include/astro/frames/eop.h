/**********************************************************************
 * eop.h
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

namespace astro {

/// Earth orientation parameters at one epoch. Angles in radians, times in seconds.
struct EopValues {
    double xp{0.0};   ///< polar motion x (rad)
    double yp{0.0};   ///< polar motion y (rad)
    double dut1{0.0}; ///< UT1 - UTC (s)
    double lod{0.0};  ///< excess length of day (s)
    double dX{0.0};   ///< celestial pole offset dX w.r.t. IAU 2006/2000A (rad)
    double dY{0.0};   ///< celestial pole offset dY w.r.t. IAU 2006/2000A (rad)
};

/// Source of Earth orientation parameters, evaluated at a UTC epoch.
class EopProvider {
public:
    virtual ~EopProvider() = default;
    /// EOP at the given MJD (UTC). May throw astro::Error if out of range.
    [[nodiscard]] virtual EopValues at(double mjdUtc) const = 0;
};

/// Returns the same values for every epoch (e.g. from a bulletin, or all zero).
class ConstantEop final : public EopProvider {
public:
    constexpr ConstantEop() noexcept = default;
    constexpr explicit ConstantEop(const EopValues& v) noexcept : v_(v) {}
    [[nodiscard]] EopValues at(double /*mjdUtc*/) const override { return v_; }

private:
    EopValues v_{};
};

} // namespace astro
