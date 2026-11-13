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

/// SGP4 error codes (Vallado): 1 mean eccentricity out of range, 2 mean motion
/// negative, 3 perturbed eccentricity out of range, 4 semi-latus rectum negative,
/// 6 satellite has decayed.
class Sgp4Error : public Error {
public:
    /// `code` is the SGP4 error code, `what` the message.
    Sgp4Error(int code, const std::string& what)
        : Error(what)
        , code_(code) {}
    /// The SGP4 error code (see the class description).
    [[nodiscard]] int code() const noexcept { return code_; }

private:
    int code_;
};

/// SGP4/SDP4 propagator for one TLE. The deep-space (SDP4) branch is selected automatically.
class Sgp4 {
public:
    /// Initialize from `tle`. Throws Sgp4Error (an astro::Error) if the elements fail at the epoch (codes 1-4; a
    /// decayed state, code 6, is accepted). A negative mean motion is code 2.
    ///
    /// Throws astro::Error if the epoch or an element the theory uses is not finite (B*, the inclination, the
    /// right ascension of the node, the eccentricity, the argument of perigee, the mean anomaly, the mean motion),
    /// or if finite elements give a non-finite state at the epoch. So a propagator that was built never starts
    /// from NaN.
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

    /// TEME position (m) and velocity (m/s) `minutes` after the TLE epoch.
    /// Throws Sgp4Error on propagation failure.
    [[nodiscard]] StateVector propagate(double minutesSinceEpoch) const;
    /// TEME state at a time (any scale).
    [[nodiscard]] StateVector propagate(const Time& t) const;
    /// Non-throwing variant: returns the SGP4 error code (0 on success). For
    /// code 6 (decayed) the state is still filled in, as in Vallado's code.
    [[nodiscard]] int propagate(double minutesSinceEpoch, StateVector& out) const noexcept;

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
