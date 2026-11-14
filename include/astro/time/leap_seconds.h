/**********************************************************************
 * leap_seconds.h
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
#include <optional>
#include <string>
#include <string_view>
#include <vector>

/// \file
/// The UTC leap-second table: built in, or loaded from IERS / IANA files.
namespace astro {

/// One step of the TAI-UTC history. From the given UTC date (0h) onward,
///   TAI - UTC = offset + (MJD_UTC - mjdRef) * rate   [seconds]
/// The rate term is non-zero only for 1961-1971 ("rubber second" UTC).
struct LeapSecondEntry {
    int    year;        ///< UTC year in which the step takes effect
    int    month;       ///< steps always occur on the 1st of a month
    double offset;      ///< seconds
    double mjdRef{0.0}; ///< reference MJD for the drift term
    double rate{0.0};   ///< seconds per day
};

/// Table of TAI-UTC offsets.
///
/// The built-in table covers 1960-01-01 through the 2017-01-01 leap second
/// (TAI-UTC = 37 s), the latest announced at the time of writing. A newer
/// table can be loaded from any of the standard files:
///  - IERS `Leap_Second.dat`
///  - IETF/NIST `leap-seconds.list`
///  - USNO `tai-utc.dat`
/// and installed globally with setLeapSecondTable().
class LeapSecondTable {
public:
    /// The built-in table.
    LeapSecondTable();
    /// From explicit entries (sorted by date here) and an optional expiration
    /// MJD. Throws astro::Error if `entries` is empty or an entry is invalid (year outside
    /// 1960-9999, month outside 1-12, non-finite values).
    explicit LeapSecondTable(std::vector<LeapSecondEntry> entries, std::optional<double> expiresMjd = std::nullopt);

    /// The built-in table (same as the default constructor).
    [[nodiscard]] static LeapSecondTable builtin();
    /// Load from a file. The format is recognised from the first data line, not
    /// from the comments. Every entry in the file is used; the history before
    /// the file's first entry (Leap_Second.dat and leap-seconds.list start in
    /// 1972, tai-utc.dat in 1961) is taken from the built-in table.
    ///
    /// Throws astro::Error if the file cannot be opened, if it is in none of the
    /// three formats or has no entries, and if a data line is malformed: fields
    /// missing, not numbers or followed by other text, a date that is not the
    /// first of a month, or a date that disagrees with the (M)JD on the same
    /// line. The message gives the line number and the line. A malformed line
    /// is never skipped, because TAI-UTC would then be wrong without notice.
    [[nodiscard]] static LeapSecondTable load(const std::string& path);
    /// As load(), from file contents already in memory.
    [[nodiscard]] static LeapSecondTable parse(std::string_view content);

    /// TAI-UTC in seconds for the UTC calendar date and fraction of day.
    /// Throws astro::Error before 1960-01-01, where UTC is undefined.
    [[nodiscard]] double taiMinusUtc(int year, int month, int day, double fractionOfDay = 0.0) const;
    /// TAI-UTC in seconds for an MJD (UTC).
    [[nodiscard]] double taiMinusUtc(double mjdUtc) const;

    /// The entries, sorted by date.
    [[nodiscard]] const std::vector<LeapSecondEntry>& entries() const noexcept { return entries_; }
    /// Expiration (MJD) declared by the loaded file, if any.
    [[nodiscard]] std::optional<double> expiresMjd() const noexcept { return expiresMjd_; }

private:
    std::vector<LeapSecondEntry> entries_;
    std::optional<double>        expiresMjd_;
};

/// The table used by all UTC conversions (thread-safe snapshot).
[[nodiscard]] std::shared_ptr<const LeapSecondTable> leapSecondTable();
/// Replace the global table used by UTC conversions (thread-safe).
void setLeapSecondTable(LeapSecondTable table);

} // namespace astro
