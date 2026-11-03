/**********************************************************************
 * tle.h
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

#include <string>
#include <string_view>
#include <vector>

#include "astro/time/time.h"

/// \file
/// Two-line element sets: parsing, formatting, checksums and Alpha-5 catalogue numbers.
namespace astro {

/// NORAD two-line element set. Angles in radians; mean motion in revolutions/day
/// as in the TLE; derivatives as printed (n-dot/2 in rev/day^2, n-ddot/6 in rev/day^3).
struct Tle {
    std::string name;                ///< satellite name (title line, leading "0 " removed); may be empty
    int         satnum{0};           ///< catalog number (Alpha-5 decoded, e.g. "A0001" -> 100001)
    char        classification{'U'}; ///< 'U' unclassified, 'C' classified, 'S' secret
    std::string intlDesignator;      ///< e.g. "98067A"
    Time        epoch;               ///< UTC
    double      meanMotionDot{0.0};  ///< first derivative / 2 (rev/day^2)
    double      meanMotionDdot{0.0}; ///< second derivative / 6 (rev/day^3)
    double      bstar{0.0};          ///< drag term (1 / Earth radii)
    int         ephemerisType{0};    ///< ephemeris type (line 1, column 63; normally 0)
    int         elementSetNumber{0}; ///< element set number (written modulo 10000)
    double      inclination{0.0};    ///< mean inclination (rad)
    double      raan{0.0};           ///< mean right ascension of the ascending node (rad)
    double      eccentricity{0.0};   ///< mean eccentricity
    double      argPerigee{0.0};     ///< mean argument of perigee (rad)
    double      meanAnomaly{0.0};    ///< mean anomaly at epoch (rad)
    double      meanMotion{0.0};     ///< rev/day (Kozai mean motion)
    int         revolutionNumber{0}; ///< revolution number at epoch (written modulo 100000)

    /// Parse the two data lines. Throws astro::Error on malformed lines or (when
    /// `verifyChecksum`) checksum mismatch.
    [[nodiscard]] static Tle
    parse(std::string_view line1, std::string_view line2, std::string_view name = "", bool verifyChecksum = true);
    /// Parse a text block of 2-line or 3-line (name + 2 lines) element sets.
    [[nodiscard]] static std::vector<Tle> parseMany(std::string_view text, bool verifyChecksum = true);

    /// @{

    /// Format back into the two 69-column lines (with checksums).
    [[nodiscard]] std::string line1() const;
    [[nodiscard]] std::string line2() const;
    /// @}
};

/// TLE line checksum (sum of digits, '-' counts as 1, modulo 10) over columns 1-68.
[[nodiscard]] int tleChecksum(std::string_view line) noexcept;
/// Decode an Alpha-5 / numeric 5-character catalog number.
[[nodiscard]] int decodeAlpha5(std::string_view s);
/// Encode a catalog number as 5 characters: zero-padded digits below 100000, Alpha-5 (a letter,
/// skipping I and O, plus 4 digits) up to 339999. Throws astro::Error above 339999.
[[nodiscard]] std::string encodeAlpha5(int satnum);

} // namespace astro
