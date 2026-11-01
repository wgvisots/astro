/**********************************************************************
 * angle.h
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

#include "astro/core/constants.h"

namespace astro {

[[nodiscard]] constexpr double deg2rad(double deg) noexcept { return deg * constants::degToRad; }
[[nodiscard]] constexpr double rad2deg(double rad) noexcept { return rad * constants::radToDeg; }
[[nodiscard]] constexpr double arcsec2rad(double as) noexcept { return as * constants::arcsecToRad; }
[[nodiscard]] constexpr double rad2arcsec(double rad) noexcept { return rad / constants::arcsecToRad; }
[[nodiscard]] constexpr double mas2rad(double mas) noexcept { return mas * constants::masToRad; }
[[nodiscard]] constexpr double rad2mas(double rad) noexcept { return rad / constants::masToRad; }

/// Wrap an angle to [0, 2*pi).
[[nodiscard]] double wrapTwoPi(double rad) noexcept;
/// Wrap an angle to [-pi, pi).
[[nodiscard]] double wrapPi(double rad) noexcept;
/// Wrap an angle to [0, 360) degrees.
[[nodiscard]] double wrap360(double deg) noexcept;
/// Wrap an angle to [-180, 180) degrees.
[[nodiscard]] double wrap180(double deg) noexcept;

/// Sexagesimal components. `negative` carries the sign so that "-00 30 00" is
/// representable; deg/hours and minutes are non-negative.
struct Sexagesimal {
    bool negative{false};
    int whole{0};   ///< degrees or hours
    int minutes{0};
    double seconds{0.0};
};

/// Hours, minutes, seconds of time to radians (15 deg per hour).
[[nodiscard]] double hmsToRad(int h, int m, double s, bool negative = false) noexcept;
/// Degrees, arcminutes, arcseconds to radians.
[[nodiscard]] double dmsToRad(int d, int m, double s, bool negative = false) noexcept;

/// Decompose radians into hours/minutes/seconds, with seconds rounded to `decimals`
/// places (carries propagate, so 59.9999 never prints as 60).
[[nodiscard]] Sexagesimal radToHms(double rad, int decimals = 3) noexcept;
/// Decompose radians into degrees/arcminutes/arcseconds.
[[nodiscard]] Sexagesimal radToDms(double rad, int decimals = 2) noexcept;

/// Parse a sexagesimal string into radians. Accepts separators made of spaces,
/// ':', 'h', 'd', 'm', 's', '\'', '"', and the degree sign; one to three fields;
/// a leading '+' or '-'. E.g. "12:34:56.7", "-05d30m00s", "12h 30m".
/// Throws astro::Error on malformed input.
[[nodiscard]] double parseHms(std::string_view text);
[[nodiscard]] double parseDms(std::string_view text);

/// Format as "HH:MM:SS.sss" (sep=':').
[[nodiscard]] std::string formatHms(double rad, int decimals = 3, char sep = ':');
/// Format as "+DD:MM:SS.ss" (sep=':'), always signed.
[[nodiscard]] std::string formatDms(double rad, int decimals = 2, char sep = ':');

} // namespace astro
