/**********************************************************************
 * helmert.h
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

/// \file
/// Seven- and fourteen-parameter Helmert (similarity) transformations between
/// terrestrial reference frames, in the linearized IERS form
///     X2 = X1 + T + s X1 + R X1,
///     R  = [ 0  -r3  r2 ;  r3  0  -r1 ;  -r2  r1  0 ]   (position-vector convention)
/// The coordinate-frame convention (used by e.g. some EPSG/NGA publications)
/// has the opposite rotation signs.
namespace astro {

enum class HelmertConvention {
    PositionVector,
    CoordinateFrame
};

BOOST_DESCRIBE_ENUM(HelmertConvention, PositionVector, CoordinateFrame)

/// Writes "HelmertConvention::<enumerator>", or "HelmertConvention::INVALID_VALUE" for any other value.
std::ostream& operator<<(std::ostream& os, HelmertConvention value);

struct HelmertParams {
    Vec3              t;      ///< translation (m)
    Vec3              r;      ///< rotation angles (rad)
    double            s{0.0}; ///< scale (unitless, e.g. 1.5 ppb -> 1.5e-9)
    HelmertConvention convention{HelmertConvention::PositionVector};
};

/// Time-dependent (14-parameter) Helmert transformation: parameters at epoch
/// `epoch` (decimal years) plus their rates per year.
struct HelmertTimeDependent {
    HelmertParams p;
    Vec3          tRate;      ///< m / yr
    Vec3          rRate;      ///< rad / yr
    double        sRate{0.0}; ///< 1 / yr
    double        epoch{2010.0};

    /// The 7 parameters at a given epoch (decimal years).
    [[nodiscard]] HelmertParams at(double decimalYear) const noexcept;
};

[[nodiscard]] Vec3 helmertTransform(const Vec3& x, const HelmertParams& p) noexcept;
/// First-order inverse transformation (all parameters negated).
[[nodiscard]] HelmertParams inverse(const HelmertParams& p) noexcept;
/// Transform position and velocity with a time-dependent Helmert at the given epoch.
[[nodiscard]] StateVector
helmertTransform(const StateVector& x, const HelmertTimeDependent& h, double decimalYear) noexcept;

} // namespace astro
