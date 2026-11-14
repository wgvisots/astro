/**********************************************************************
 * sgp4.h
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
#include <memory>

#include <boost/describe/enum.hpp>

#include "astro/core/error.h"
#include "astro/core/vec.h"
#include "astro/orbit/tle.h"

/// \file
/// SGP4/SDP4 analytic propagator for two-line element sets, ported from the reference
/// implementation of Vallado, Crawford, Hujsak & Kelso (2006), "Revisiting Spacetrack Report #3"
/// (AIAA 2006-6753). Output is in the TEME frame (see astro::Frame::TEME).
namespace astro {

/// Earth constants used by SGP4. TLEs are generated with WGS-72.
enum class Sgp4Gravity {
    WGS72Old,
    WGS72,
    WGS84
};

BOOST_DESCRIBE_ENUM(Sgp4Gravity, WGS72Old, WGS72, WGS84)

/// Writes "Sgp4Gravity::<enumerator>", or "Sgp4Gravity::INVALID_VALUE" for any other value.
std::ostream& operator<<(std::ostream& os, Sgp4Gravity value);

/// 'Improved' (Vallado) or 'AFSPC' compatibility for the GMST at epoch and node wrapping.
enum class Sgp4OpsMode {
    Improved,
    Afspc
};

BOOST_DESCRIBE_ENUM(Sgp4OpsMode, Improved, Afspc)

/// Writes "Sgp4OpsMode::<enumerator>", or "Sgp4OpsMode::INVALID_VALUE" for any other value.
std::ostream& operator<<(std::ostream& os, Sgp4OpsMode value);

/// Outcome of one SGP4 propagation. The underlying values are the error codes of the reference code
/// (Vallado's sgp4, 1 to 6, with no 5), so `static_cast<int>(status)` compares with its output. One more
/// is this library's own: 7, TimeNotFinite.
enum class Sgp4Status {
    Ok                              = 0, ///< success
    MeanEccentricityOutOfRange      = 1, ///< mean eccentricity e >= 1 or e < -0.001
    NegativeMeanMotion              = 2, ///< mean motion less than zero
    PerturbedEccentricityOutOfRange = 3, ///< perturbed eccentricity outside [0, 1]
    NegativeSemiLatusRectum         = 4, ///< semi-latus rectum less than zero
    Decayed                         = 6, ///< the satellite has decayed (radius below one Earth radius)
    TimeNotFinite                   = 7, ///< the time asked for is NaN or infinite (not a reference code)
};

BOOST_DESCRIBE_ENUM(Sgp4Status,
                    Ok,
                    MeanEccentricityOutOfRange,
                    NegativeMeanMotion,
                    PerturbedEccentricityOutOfRange,
                    NegativeSemiLatusRectum,
                    Decayed,
                    TimeNotFinite)

/// Writes "Sgp4Status::<enumerator>", or "Sgp4Status::INVALID_VALUE" for any other value.
std::ostream& operator<<(std::ostream& os, Sgp4Status value);

/// Thrown by Sgp4 when the theory fails: an astro::Error that also carries the Sgp4Status, so that a
/// caller can tell a decayed satellite from bad elements.
class Sgp4Error : public Error {
public:
    /// `code` is the status that failed (never Sgp4Status::Ok), `what` the message.
    Sgp4Error(Sgp4Status code, const std::string& what)
        : Error(what)
        , code_(code) {}
    /// The status that failed.
    [[nodiscard]] Sgp4Status code() const noexcept { return code_; }

private:
    Sgp4Status code_;
};

/// The result of Sgp4::tryPropagate and Sgp4::tryPropagateMinutes: the state and how it was obtained.
struct Sgp4Result {
    StateVector state;                  ///< TEME position (m) and velocity (m/s); see `status` for its meaning
    Sgp4Status  status{Sgp4Status::Ok}; ///< Ok, or why the state is not a valid one
};

/// SGP4/SDP4 propagator for one TLE. The deep-space (SDP4) branch is selected automatically.
///
/// Times are seconds, as for every propagator in the library: propagate(double) takes seconds since
/// the element set's epoch. The reference code (Vallado's sgp4) and most other SGP4 implementations
/// take minutes instead; propagateMinutes takes minutes, for comparison with them and with their
/// verification output.
class Sgp4 {
public:
    /// Initialize from `tle`. Throws Sgp4Error (an astro::Error) if the elements fail at the epoch (statuses
    /// 1-4 of Sgp4Status; a decayed state, Sgp4Status::Decayed, is accepted). A negative mean motion is
    /// Sgp4Status::NegativeMeanMotion.
    ///
    /// Throws astro::Error if the epoch or an element the theory uses is not finite (B*, the inclination, the
    /// right ascension of the node, the eccentricity, the argument of perigee, the mean anomaly, the mean motion),
    /// or if finite elements give a non-finite state at the epoch. So a propagator that was built never starts
    /// from NaN. Also throws if the epoch cannot be expressed as UTC: one before 1960, or one given in UT1.
    ///
    /// The theory counts the epoch in days of 86400 s, as the element set gives it. On a day that ends with a
    /// leap second the count is taken from the time of day of `tle.epoch` over 86400 s, as the reference code
    /// takes it from the epoch field.
    explicit Sgp4(const Tle& tle, Sgp4Gravity gravity = Sgp4Gravity::WGS72, Sgp4OpsMode mode = Sgp4OpsMode::Improved);
    ~Sgp4();
    /// Deep copy.
    Sgp4(const Sgp4&);
    /// Deep copy.
    Sgp4& operator=(const Sgp4&);
    /// Move; the moved-from object may only be destroyed or assigned to.
    Sgp4(Sgp4&&) noexcept;
    /// Move; the moved-from object may only be destroyed or assigned to.
    Sgp4& operator=(Sgp4&&) noexcept;

    /// TEME position (m) and velocity (m/s) `dt` seconds after the TLE epoch: propagateMinutes(dt / 60).
    /// Throws Sgp4Error on propagation failure: for any status but Sgp4Status::Ok, Decayed and
    /// TimeNotFinite (`dt` is not finite) included.
    [[nodiscard]] StateVector propagate(double dt) const;
    /// TEME state at a time (any scale). The seconds since the epoch are the elapsed time, leap
    /// seconds included; software that subtracts Julian dates instead differs by one second (7.7 km
    /// in low Earth orbit) across a leap second. Throws Sgp4Error as propagate(dt) does, and
    /// astro::Error if t cannot be subtracted from the epoch (a UT1 time, or a UTC time before 1960;
    /// see Time).
    [[nodiscard]] StateVector propagate(const Time& t) const;
    /// Non-throwing variant of propagate(dt), for a loop that must handle a failure at every sample:
    /// the state `dt` seconds after the TLE epoch and its status. For Sgp4Status::Decayed the state is
    /// still filled in, as in Vallado's code. For Sgp4Status::TimeNotFinite (`dt` is NaN or infinite)
    /// every component of the state is NaN.
    [[nodiscard]] Sgp4Result tryPropagate(double dt) const noexcept;
    /// TEME position (m) and velocity (m/s) `minutes` after the TLE epoch, the time argument of the
    /// reference code. Throws Sgp4Error as propagate(dt) does.
    [[nodiscard]] StateVector propagateMinutes(double minutes) const;
    /// Non-throwing variant of propagateMinutes: the state `minutes` after the TLE epoch and its status,
    /// as tryPropagate(dt) gives them.
    [[nodiscard]] Sgp4Result tryPropagateMinutes(double minutes) const noexcept;

    /// The element set this propagator was built from.
    [[nodiscard]] const Tle& tle() const noexcept;
    /// The TLE epoch (UTC).
    [[nodiscard]] const Time& epoch() const noexcept;
    /// True if the deep-space (SDP4) branch is used (period >= 225 minutes).
    [[nodiscard]] bool isDeepSpace() const noexcept;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace astro
