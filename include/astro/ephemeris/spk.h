/**********************************************************************
 * spk.h
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

#include <memory>
#include <string>
#include <vector>

#include "astro/ephemeris/ephemeris.h"

/// \file
/// Reader for JPL SPK ephemeris kernels (.bsp; DAF segment types 2 and 3).
namespace astro {

/// Description of one SPK segment.
struct SpkSegmentInfo {
    int         target{0};
    int         center{0};
    int         frame{0};     ///< NAIF frame id (1 = J2000/ICRF)
    int         type{0};      ///< SPK data type
    double      startEt{0.0}; ///< TDB seconds past J2000
    double      endEt{0.0};
    std::string name;
};

/// Reader for JPL/NAIF SPK (.bsp) ephemeris files, e.g. DE440/DE441.
///
/// Supports Chebyshev segment types 2 (position) and 3 (position + velocity),
/// either byte order. Bodies are chained through their centers, so e.g. the
/// Moon relative to the Earth is resolved via the Earth-Moon barycenter. When
/// several segments cover an epoch, the one latest in the file wins (NAIF rule).
/// Records are read on demand; the file stays open. Thread-safe.
class SpkEphemeris final : public Ephemeris {
public:
    explicit SpkEphemeris(const std::string& path);
    ~SpkEphemeris() override;
    SpkEphemeris(SpkEphemeris&&) noexcept;
    SpkEphemeris& operator=(SpkEphemeris&&) noexcept;

    [[nodiscard]] StateVector state(Body target, Body center, const Time& t) const override;
    /// State by NAIF id at an ephemeris time (TDB seconds past J2000); meters, m/s.
    [[nodiscard]] StateVector stateById(int target, int center, double et) const;

    [[nodiscard]] const std::vector<SpkSegmentInfo>& segments() const noexcept;
    [[nodiscard]] bool                               bigEndian() const noexcept;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

/// TDB seconds past J2000.0 for a time in any scale.
[[nodiscard]] double ephemerisTime(const Time& t);

} // namespace astro
