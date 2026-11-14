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
    int         target{0};    ///< NAIF id of the target body
    int         center{0};    ///< NAIF id of the center body
    int         frame{0};     ///< NAIF frame id (1 = J2000/ICRF)
    int         type{0};      ///< SPK data type
    double      startEt{0.0}; ///< TDB seconds past J2000
    double      endEt{0.0};   ///< end of coverage, TDB seconds past J2000
    std::string name;         ///< segment name from the file (trailing blanks removed)
};

/// Reader for JPL/NAIF SPK (.bsp) ephemeris files, e.g. DE440/DE441.
///
/// Supports Chebyshev segment types 2 (position) and 3 (position + velocity),
/// either byte order, in the J2000 reference frame (NAIF frame id 1; in the JPL
/// planetary kernels its axes are those of the ICRF). Bodies are chained through
/// their centers, so e.g. the Moon relative to the Earth is resolved via the
/// Earth-Moon barycenter. When several segments cover an epoch, the one latest in
/// the file wins (NAIF rule). Records are read on demand; the file stays open.
/// Thread-safe.
///
/// A kernel may also hold segments of other types or in other frames. They are
/// listed by segments() and do no harm until one is needed for a state: state()
/// and stateById() then throw astro::Error. A segment in another frame is never
/// returned as if it were in the ICRF.
///
/// A body is resolved from its own segments only: nothing stands in for a body
/// the kernel lacks. DE430 and DE440 have the Mars system barycenter
/// (Body::MarsBarycenter, NAIF id 4) but not Mars (Body::Mars, id 499), so
/// Body::Mars throws with those kernels, and the message names the barycenter.
class SpkEphemeris final : public Ephemeris {
public:
    /// Open the kernel at `path` and read its segment summaries. Throws astro::Error if the file
    /// cannot be opened or read, is not a DAF/SPK file in IEEE format, or has no segments. The file
    /// may begin with the id word "DAF/SPK " or with the "NAIF/DAF" of older kernels (de405.bsp, for
    /// one); for those, which do not name their byte order, it is found from the file record.
    explicit SpkEphemeris(const std::string& path);
    ~SpkEphemeris() override;
    /// Move; the moved-from object may only be destroyed or assigned to.
    SpkEphemeris(SpkEphemeris&&) noexcept;
    /// Move; the moved-from object may only be destroyed or assigned to.
    SpkEphemeris& operator=(SpkEphemeris&&) noexcept;

    /// As stateById(), with the NAIF ids of `target` and `center` and the ephemeris time of `t`.
    /// Throws astro::Error as stateById() and ephemerisTime() do.
    [[nodiscard]] StateVector state(Body target, Body center, const Time& t) const override;
    /// State by NAIF id at an ephemeris time (TDB seconds past J2000); meters, m/s, ICRF axes.
    /// Throws astro::Error if a body on the way from `target` to `center` is not in the kernel or
    /// has no segment covering `et`, or if that segment is of an unsupported type or not in the
    /// J2000 frame. Also throws if a record cannot be read from the file or is corrupt.
    [[nodiscard]] StateVector stateById(int target, int center, double et) const;

    /// Summaries of all segments, in file order.
    [[nodiscard]] const std::vector<SpkSegmentInfo>& segments() const noexcept;
    /// True if the file is big-endian (BIG-IEEE).
    [[nodiscard]] bool bigEndian() const noexcept;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

/// TDB seconds past J2000.0 for a time in any scale. Throws astro::Error if `t` cannot be
/// converted to TDB: a UTC time before 1960, or a UT1 time (there is no Earth orientation data here).
[[nodiscard]] double ephemerisTime(const Time& t);

} // namespace astro
