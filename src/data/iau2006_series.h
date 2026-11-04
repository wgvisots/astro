/**********************************************************************
 * iau2006_series.h  (internal)
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

#include <cstddef>
#include <cstdint>

namespace astro::detail {

/// One term of a Poisson series: (sinCoef sin(ARG) + cosCoef cos(ARG)) * t^power,
/// ARG = sum n[i] * F_i over the 14 fundamental arguments. Coefficients in microarcseconds.
struct SeriesTerm {
    std::int8_t power;
    std::int8_t n[14];
    double      sinCoef;
    double      cosCoef;
};

/// IAU 2000B luni-solar term (0.1 microarcsecond units).
struct Nut2000BTerm {
    int    n[5];        ///< l, l', F, D, Om
    double ps, pst, pc; ///< longitude: sin, t*sin, cos
    double ec, ect, es; ///< obliquity: cos, t*cos, sin
};

extern const double kCipXPoly[6];
extern const double kCipYPoly[6];
extern const double kSXY2Poly[6];

extern const SeriesTerm  kCipX[];
extern const std::size_t kCipXCount;
extern const SeriesTerm  kCipY[];
extern const std::size_t kCipYCount;
extern const SeriesTerm  kSXY2[];
extern const std::size_t kSXY2Count;
extern const SeriesTerm  kNutLongitude[];
extern const std::size_t kNutLongitudeCount;
extern const SeriesTerm  kNutObliquity[];
extern const std::size_t kNutObliquityCount;

extern const Nut2000BTerm kNut2000B[77];

} // namespace astro::detail
