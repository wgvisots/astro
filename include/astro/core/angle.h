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

/// \file
/// Angle units and wrapping, and sexagesimal (HMS/DMS) parsing and formatting.
namespace astro {

// The unit conversions below are functions of an angle. The factors of the same names in
// astro::constants (constants::degToRad, ...) are the numbers they multiply or divide by.

/// Degrees to radians (a function): `deg` (degrees) times the factor constants::degToRad.
[[nodiscard]] constexpr double degToRad(double deg) noexcept {
    return deg * constants::degToRad;
}
/// Radians to degrees (a function): `rad` (radians) times the factor constants::radToDeg.
[[nodiscard]] constexpr double radToDeg(double rad) noexcept {
    return rad * constants::radToDeg;
}
/// Arcseconds to radians (a function): `as` (arcseconds) times the factor constants::arcsecToRad.
[[nodiscard]] constexpr double arcsecToRad(double as) noexcept {
    return as * constants::arcsecToRad;
}
/// Radians to arcseconds (a function): `rad` (radians) divided by the factor constants::arcsecToRad.
[[nodiscard]] constexpr double radToArcsec(double rad) noexcept {
    return rad / constants::arcsecToRad;
}
/// Milliarcseconds to radians (a function): `mas` (milliarcseconds) times the factor constants::masToRad.
[[nodiscard]] constexpr double masToRad(double mas) noexcept {
    return mas * constants::masToRad;
}
/// Radians to milliarcseconds (a function): `rad` (radians) divided by the factor constants::masToRad.
[[nodiscard]] constexpr double radToMas(double rad) noexcept {
    return rad / constants::masToRad;
}

/// Wrap an angle to [0, 2*pi). A value already in that range is returned unchanged, bit for bit.
/// Not-a-number (NaN) and infinity give NaN.
[[nodiscard]] double wrapTwoPi(double rad) noexcept;
/// Wrap an angle to [-pi, pi). A value already in that range is returned unchanged, bit for bit.
/// NaN and infinity give NaN.
[[nodiscard]] double wrapPi(double rad) noexcept;
/// Wrap an angle to [0, 360) degrees. A value already in that range is returned unchanged, bit for bit.
/// NaN and infinity give NaN.
[[nodiscard]] double wrap360(double deg) noexcept;
/// Wrap an angle to [-180, 180) degrees. A value already in that range is returned unchanged, bit for bit.
/// NaN and infinity give NaN.
[[nodiscard]] double wrap180(double deg) noexcept;

/// Sexagesimal components. `negative` carries the sign so that "-00 30 00" is
/// representable; deg/hours and minutes are non-negative.
struct Sexagesimal {
    bool   negative{false}; ///< true if the angle is negative (never set for zero)
    int    whole{0};        ///< degrees or hours
    int    minutes{0};      ///< minutes (of arc or of time), 0-59
    double seconds{0.0};    ///< seconds (of arc or of time), [0, 60)
};

/// Hours, minutes, seconds of time to radians (15 deg per hour).
///
/// The sign belongs to the whole angle: the result is negative if `h` is negative or `negative` is
/// true, and `h`, `m` and `s` are then used as magnitudes (the signs of `m` and `s` are ignored). So
/// hmsToRad(-1, 30, 0) is -1.5 hours. An angle between 0 and -1 hour needs `negative`, because an `int`
/// has no negative zero: hmsToRad(0, 30, 0, true) is -0.5 hours.
[[nodiscard]] double hmsToRad(int h, int m, double s, bool negative = false) noexcept;
/// Degrees, arcminutes, arcseconds to radians.
///
/// The sign is handled as in hmsToRad: dmsToRad(-5, 30, 0) is -5.5 degrees, and
/// dmsToRad(0, 30, 0, true) is -0.5 degrees.
[[nodiscard]] double dmsToRad(int d, int m, double s, bool negative = false) noexcept;

/// Decompose radians into hours/minutes/seconds, with seconds rounded to `decimals`
/// places (carries propagate, so 59.9999 never prints as 60).
///
/// An angle with no such decomposition is reported in `seconds`, with `whole` and `minutes` zero:
/// `seconds` is NaN if `rad` is not-a-number (NaN), and infinity if `rad` is infinite or so large that
/// the hours do not fit in an `int` (`negative` still gives the sign).
[[nodiscard]] Sexagesimal radToHms(double rad, int decimals = 3) noexcept;
/// Decompose radians into degrees/arcminutes/arcseconds. NaN, infinity and an angle whose degrees do
/// not fit in an `int` are reported as in radToHms.
[[nodiscard]] Sexagesimal radToDms(double rad, int decimals = 2) noexcept;

/// Parse an hours/minutes/seconds string into radians, e.g. "12:34:56.7", "12h 34m 56.7s", "12 34.5".
///
/// The text holds one, two or three numbers: hours, then minutes, then seconds.
/// - A '+' or '-' may come first. It applies to the whole angle, so "-0:30" is minus half an hour.
/// - A number is plain decimal digits. Only the last number may have a fraction ("12:34.5" is 12 h
///   34.5 min; "12.5:30" is an error). Minutes and seconds must be below 60; the first number has no limit.
/// - Between two numbers there is exactly one separator: a ':', a ',', a unit mark, or blank space.
/// - A number may be followed by its unit mark: 'h' after hours, 'm' after minutes, 's' after seconds,
///   in either case.
/// - Spaces and tabs may begin or end the text, follow the sign, and surround a separator or a unit mark.
///
/// Throws astro::Error for anything else: a value too large to be an angle (hundreds of digits),
/// an empty string or an empty field ("5::30"), a fourth
/// number, an exponent ("1e2"), "nan" or "inf", a sign that is not in front, a unit mark on the wrong
/// number ("30m" alone) or of the wrong quantity ("12d", or an arc mark such as "12h30'"), and any
/// other character (a line end included).
[[nodiscard]] double parseHms(std::string_view text);
/// Parse a degrees/arcminutes/arcseconds string into radians, e.g. "-05:30:00", "-05d30m00s",
/// "+45d 15' 36.2\"".
///
/// The grammar and the errors are those of parseHms, with these unit marks: 'd' or the degree sign
/// (UTF-8 or ISO 8859-1) after degrees, 'm' or an apostrophe after arcminutes, 's' or a double quote after
/// arcseconds. An hours mark ("12h") is an error. Throws astro::Error as parseHms does.
[[nodiscard]] double parseDms(std::string_view text);

/// Format as "HH:MM:SS.sss" (sep=':'). The angle is not wrapped: 2 pi prints as 24:00:00.000 and
/// larger angles as more than 24 hours. A NaN angle gives "nan"; an angle that is infinite, or too
/// large for radToHms, gives "inf" or "-inf".
[[nodiscard]] std::string formatHms(double rad, int decimals = 3, char sep = ':');
/// Format as "+DD:MM:SS.ss" (sep=':'), always signed. A NaN angle gives "nan"; an angle that is
/// infinite, or too large for radToDms, gives "+inf" or "-inf".
[[nodiscard]] std::string formatDms(double rad, int decimals = 2, char sep = ':');

} // namespace astro
