/**********************************************************************
 * gravity.cpp
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

#include "astro/orbit/gravity.h"

#include <cmath>

#include "gravity_detail.h"

namespace astro {

Vec3 pointMassAcceleration(const Vec3& r, double mu) noexcept {
    const double rm = norm(r);
    return r * (-mu / (rm * rm * rm));
}

Vec3 zonalAcceleration(const Vec3& r, const CentralBody& body) noexcept {
    Vec3 a;
    detail::gravity<false>(r, body.J.data(), body.J.size(), body.mu, body.re, false, &a, nullptr);
    return a;
}

Vec3 thirdBodyAcceleration(const Vec3& r, const Vec3& s, double muBody) noexcept {
    const Vec3   d  = s - r;
    const double dm = norm(d), sm = norm(s);
    return muBody * (d / (dm * dm * dm) - s / (sm * sm * sm));
}

} // namespace astro
