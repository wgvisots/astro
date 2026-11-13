/**********************************************************************
 * eop.h
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
#include <string>
#include <string_view>
#include <vector>

#include <boost/describe/enum.hpp>

/// \file
/// Earth orientation parameters (UT1 - UTC, polar motion, celestial pole offsets): providers and IERS file readers.
namespace astro {

/// Earth orientation parameters at one epoch. Angles in radians, times in seconds.
struct EopValues {
    double xp{0.0};   ///< polar motion x (rad)
    double yp{0.0};   ///< polar motion y (rad)
    double dut1{0.0}; ///< UT1 - UTC (s)
    double lod{0.0};  ///< excess length of day (s)
    double dX{0.0};   ///< celestial pole offset dX w.r.t. IAU 2006/2000A (rad)
    double dY{0.0};   ///< celestial pole offset dY w.r.t. IAU 2006/2000A (rad)
};

/// Source of Earth orientation parameters, evaluated at a UTC epoch.
class EopProvider {
public:
    virtual ~EopProvider() = default;
    /// EOP at the given MJD (UTC). May throw astro::Error if out of range.
    [[nodiscard]] virtual EopValues at(double mjdUtc) const = 0;
};

/// Returns the same values for every epoch (e.g. from a bulletin, or all zero
/// for an approximate UT1 = UTC, no polar motion transformation).
class ConstantEop final : public EopProvider {
public:
    /// All parameters zero (UT1 = UTC, no polar motion, no pole offsets).
    constexpr ConstantEop() noexcept = default;
    /// Always returns `v`.
    constexpr explicit ConstantEop(const EopValues& v) noexcept
        : v_(v) {}
    [[nodiscard]] EopValues at(double /*mjdUtc*/) const override { return v_; }

private:
    EopValues v_{};
};

/// EOP at one tabulated epoch (MJD, UTC, normally 0h).
struct EopRecord {
    double    mjd{0.0}; ///< epoch, MJD (UTC)
    EopValues values;   ///< parameters at that epoch
};

/// Interpolation between tabulated EOP records.
enum class EopInterpolation {
    Linear,    ///< linear between the two neighbouring records
    Lagrange4, ///< 4-point Lagrange (cubic), recommended for daily data
};

BOOST_DESCRIBE_ENUM(EopInterpolation, Linear, Lagrange4)

/// Writes "EopInterpolation::<enumerator>", or "EopInterpolation::INVALID_VALUE" for any other value.
std::ostream& operator<<(std::ostream& os, EopInterpolation value);

/// Tabulated Earth orientation parameters with interpolation.
///
/// UT1-UTC is interpolated as UT1-TAI (using the current leap-second table)
/// so leap-second steps do not corrupt the interpolation. Sub-daily ocean-tide
/// and libration corrections (IERS 2010 Ch. 5/8) are not applied.
class EopTable final : public EopProvider {
public:
    /// From records in any order (sorted here). Throws astro::Error if there
    /// are fewer than two records, if an epoch is not finite, or if two records
    /// have the same epoch (the message names it): the interpolation needs
    /// distinct epochs.
    explicit EopTable(std::vector<EopRecord> records, EopInterpolation interp = EopInterpolation::Lagrange4);

    /// IERS Rapid Service `finals2000A.all` / `finals2000A.data` / `finals2000A.daily`
    /// (IAU 2000A version, with dX/dY). Uses the Bulletin A columns; rows without
    /// polar motion are skipped, missing LOD/dX/dY are taken as zero.
    [[nodiscard]] static EopTable loadFinals2000A(const std::string& path,
                                                  EopInterpolation   interp = EopInterpolation::Lagrange4);
    /// As loadFinals2000A(), from file contents already in memory. Throws
    /// astro::Error if no usable rows are found, or as the constructor does
    /// (fewer than two rows, a repeated epoch).
    [[nodiscard]] static EopTable parseFinals2000A(std::string_view content,
                                                   EopInterpolation interp = EopInterpolation::Lagrange4);
    /// IERS EOP C04 series, either the classic "14" format
    /// (YYYY MM DD MJD x y UT1-UTC LOD dX dY ...) or the newer "20" format
    /// (YYYY MM DD HH MJD x y UT1-UTC dX dY xrt yrt LOD ...). Angles in arcsec.
    [[nodiscard]] static EopTable loadC04(const std::string& path,
                                          EopInterpolation   interp = EopInterpolation::Lagrange4);
    /// As loadC04(), from file contents already in memory. Throws astro::Error
    /// if no usable rows are found, or as the constructor does (fewer than two
    /// rows, a repeated epoch).
    [[nodiscard]] static EopTable parseC04(std::string_view content,
                                           EopInterpolation interp = EopInterpolation::Lagrange4);

    /// Interpolated EOP. Throws astro::Error outside [firstMjd(), lastMjd()].
    [[nodiscard]] EopValues at(double mjdUtc) const override;

    /// MJD (UTC) of the first / last record.
    /// @{
    [[nodiscard]] double firstMjd() const noexcept { return records_.front().mjd; }
    [[nodiscard]] double lastMjd() const noexcept { return records_.back().mjd; }
    /// @}
    /// The records, sorted by MJD.
    [[nodiscard]] const std::vector<EopRecord>& records() const noexcept { return records_; }

private:
    std::vector<EopRecord> records_;
    std::vector<double>    ut1MinusTai_; // per record
    EopInterpolation       interp_;
};

} // namespace astro
