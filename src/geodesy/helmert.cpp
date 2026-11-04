/**********************************************************************
 * helmert.cpp
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

#include "astro/geodesy/helmert.h"

#include "../core/enum_io.h"

namespace astro {

namespace {
Mat3 rotationMatrix(const Vec3& r, HelmertConvention c) {
    const Mat3 R = Mat3::skew(r); // [0 -r3 r2; r3 0 -r1; -r2 r1 0]
    return c == HelmertConvention::PositionVector ? R : R.transpose();
}
} // namespace

HelmertParams HelmertTimeDependent::at(double decimalYear) const noexcept {
    const double dt = decimalYear - epoch;
    return {p.t + tRate * dt, p.r + rRate * dt, p.s + sRate * dt, p.convention};
}

Vec3 helmertTransform(const Vec3& x, const HelmertParams& p) noexcept {
    return x + p.t + p.s * x + rotationMatrix(p.r, p.convention) * x;
}

HelmertParams inverse(const HelmertParams& p) noexcept {
    return {-p.t, -p.r, -p.s, p.convention};
}

StateVector helmertTransform(const StateVector& x, const HelmertTimeDependent& h, double decimalYear) noexcept {
    const HelmertParams p = h.at(decimalYear);
    StateVector         out;
    out.r = helmertTransform(x.r, p);
    // V2 = V1 + Tdot + sdot X1 + Rdot X1 (+ s V1 + R V1, retained for consistency)
    out.v = x.v + h.tRate + h.sRate * x.r + rotationMatrix(h.rRate, p.convention) * x.r + p.s * x.v
          + rotationMatrix(p.r, p.convention) * x.v;
    return out;
}

std::ostream& operator<<(std::ostream& os, HelmertConvention value) {
    return detail::writeEnum(os, "HelmertConvention", value);
}

} // namespace astro
