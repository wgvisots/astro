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
    ///
    /// With `verifyChecksum` each line must have all 69 columns and a digit in column 69: a line without a
    /// checksum is refused, not passed. Without it, shorter lines are read: line 1 needs at least 64
    /// columns and line 2 at least 63, not counting blanks at the end, and a field that is missing is
    /// zero. Anything after column 69 is ignored.
    ///
    /// The catalog number is five digits or Alpha-5 (see decodeAlpha5()); leading blanks in place of
    /// zeros are accepted. It must be the same on both lines.
    ///
    /// The fields are read by column, as a fixed-format read would take them. A number field may hold
    /// digits, one sign, one decimal point and blank padding; anything else ("nan", "inf", an exponent,
    /// hexadecimal, a sign in a whole-number field) throws astro::Error naming the field. An all-blank
    /// number field is zero. In the fields with an implied decimal point (the eccentricity, and the
    /// mantissas of B* and of the second derivative of the mean motion) a blank column is a zero digit,
    /// as in the SGP4 reference code: an eccentricity field of " 006703" is 0.0006703.
    ///
    /// The epoch field is the day of the year (from 1) and a fraction of an 86400 s day counted from 0h UTC.
    /// That holds on a day that ends with a leap second too: element sets are produced that way, and the
    /// SGP4 reference code reads them that way. So "16366.50000000" is 2016-12-31 12:00:00 UTC exactly. A
    /// day that the year does not have throws astro::Error. Two-digit years 57-99 are 1957-1999 and
    /// 00-56 are 2000-2056; an epoch before UTC begins in the leap-second table (1960-01-01 in the
    /// built-in one) throws astro::Error, because astro::Time cannot hold it as UTC.
    [[nodiscard]] static Tle
    parse(std::string_view line1, std::string_view line2, std::string_view name = "", bool verifyChecksum = true);
    /// Parse a text block of 2-line or 3-line (name + 2 lines) element sets, with either line ending.
    ///
    /// Blank lines, and lines whose first character is '#' (comments), are skipped wherever they are. A
    /// data line starts with "1 " or "2 " and is at least 64 (line 1) or 63 (line 2) columns long; any
    /// other line is a name, and the last name before a line 1 belongs to that set. Throws astro::Error
    /// with the line number in `text` (from 1) for a line 1 that is not followed by its line 2 and for a
    /// line 2 without a line 1 before it, and whatever parse() throws for a set.
    [[nodiscard]] static std::vector<Tle> parseMany(std::string_view text, bool verifyChecksum = true);

    /// @{

    /// Format back into the two 69-column lines (with checksums). Each value is rounded to the nearest one
    /// its field can hold, and a value that no rounding can fit throws astro::Error naming the field, so a
    /// line that is returned is well formed and parse() reads it back.
    ///
    /// The epoch is written as parse() reads it: the time of day as a fraction of 86400 s, in steps of
    /// 1e-8 day, also on a day that ends with a leap second. An epoch that rounds up to the end of its
    /// day is written as 0h of the next day (of the next year after 31 December), except just before a
    /// leap second, where the last step of the day is the nearest. Throws for an epoch inside a leap
    /// second, which no day fraction names, and for a year outside 1957 to 2056, the range of the
    /// two-digit field.
    ///
    /// B* and the second derivative of the mean motion are written with five digits and a power of ten
    /// from -9 to +9. Below 1e-10 the digits get leading zeros (1e-11 is " 01000-9"), in steps of 1e-14,
    /// and below 5e-15 the field is zero. Throws for a value that rounds above 0.99999e9.
    ///
    /// The right ascension of the node, the argument of perigee and the mean anomaly are wrapped into
    /// [0, 360) deg; one that rounds up to 360.0000 is written as 0.0000. The inclination is written as it
    /// is.
    ///
    /// Also thrown for: any number that is not finite; a catalog number outside 0 to 339999; a
    /// classification that is not a printable character; an international designator of more than 8
    /// characters; a first derivative of the mean motion that rounds to 1 rev/day^2 or more in magnitude
    /// (the field has no digit before the decimal point); an ephemeris type outside 0 to 9; a negative
    /// element set or revolution number; an eccentricity that is negative or rounds to 1 or more; an
    /// inclination or a mean motion with more digits before the decimal point than the field has columns
    /// (100 rev/day or more).
    [[nodiscard]] std::string line1() const;
    [[nodiscard]] std::string line2() const;
    /// @}
};

/// TLE line checksum (sum of digits, '-' counts as 1, modulo 10) over columns 1-68.
[[nodiscard]] int tleChecksum(std::string_view line) noexcept;
/// Decode an Alpha-5 / numeric 5-character catalog number: five digits, or a capital letter other than
/// I and O followed by four digits ("A0001" is 100001, "Z9999" is 339999). Throws astro::Error for
/// anything else: another length, blanks, a sign, a lower-case letter.
[[nodiscard]] int decodeAlpha5(std::string_view s);
/// Encode a catalog number as 5 characters: zero-padded digits below 100000, Alpha-5 (a letter,
/// skipping I and O, plus 4 digits) up to 339999. Throws astro::Error below 0 and above 339999.
[[nodiscard]] std::string encodeAlpha5(int satnum);

} // namespace astro
