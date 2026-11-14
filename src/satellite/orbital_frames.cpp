/**********************************************************************
 * orbital_frames.cpp
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

#include "astro/satellite/orbital_frames.h"

namespace astro {

Mat3 ricMatrix(const StateVector& s) noexcept {
    const Vec3 R = normalized(s.r);
    const Vec3 C = normalized(cross(s.r, s.v));
    const Vec3 I = cross(C, R);
    return Mat3::fromRows(R, I, C);
}

Mat3 lvlhMatrix(const StateVector& s) noexcept {
    const Vec3 z = -normalized(s.r);
    const Vec3 y = -normalized(cross(s.r, s.v));
    const Vec3 x = cross(y, z);
    return Mat3::fromRows(x, y, z);
}

Mat3 vncMatrix(const StateVector& s) noexcept {
    const Vec3 V = normalized(s.v);
    const Vec3 N = normalized(cross(s.r, s.v));
    const Vec3 C = cross(V, N);
    return Mat3::fromRows(V, N, C);
}

Mat3 ntwMatrix(const StateVector& s) noexcept {
    const Vec3 T = normalized(s.v);
    const Vec3 W = normalized(cross(s.r, s.v));
    const Vec3 N = cross(T, W);
    return Mat3::fromRows(N, T, W);
}

StateVector relativeStateRic(const StateVector& chief, const StateVector& deputy) noexcept {
    const Mat3 M     = ricMatrix(chief);
    const Vec3 omega = cross(chief.r, chief.v) / norm2(chief.r);
    const Vec3 dr    = deputy.r - chief.r;
    const Vec3 dv    = deputy.v - chief.v;
    return {M * dr, M * (dv - cross(omega, dr))};
}

StateVector relativeRicToDeputy(const StateVector& chief, const StateVector& rel) noexcept {
    const Mat3 Mt    = ricMatrix(chief).transpose();
    const Vec3 omega = cross(chief.r, chief.v) / norm2(chief.r);
    const Vec3 dr    = Mt * rel.r;
    return {chief.r + dr, chief.v + Mt * rel.v + cross(omega, dr)};
}

} // namespace astro
