/**********************************************************************
 * tle_detail.h
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

#include <cmath>

#include "astro/time/julian_date.h"
#include "astro/time/time.h"

namespace astro::detail {

/// (86400 s) / (length of the UTC day that starts at the Julian date `midnight`).
///
/// The epoch of a two-line element set is a day of the year and a fraction of an 86400 s day, counted
/// from 0h UTC whatever the length of that day: the sets are produced from, and read back into, a
/// Julian date that knows nothing of leap seconds (in the reference code, twoline2rv turns the
/// fraction into hours, minutes and seconds of a 24 h day, and SGP4 counts its epoch in such days).
/// astro::Time holds UTC differently: as a quasi-Julian date whose day fraction is the time of day
/// over the length of *that* day, which is 86401 s when the day ends with a leap second (see
/// time.h). The two fractions differ by up to 1 s on such a day: 7.7 km along the track in low Earth
/// orbit.
///
/// This is the factor between them. 12:00:00 is 43200 s into any day, so it is twice the fraction at
/// which astro::Time puts noon. It is exactly 1 on an ordinary day. The element-set reader and
/// writer (tle.cpp) and Sgp4 (sgp4.cpp) all use it, so that they agree.
[[nodiscard]] inline double tleDayPerUtcDay(double midnight) {
    int year = 0, month = 0, day = 0;
    julianDayNumberToDate(std::llround(midnight + 0.5), year, month, day);
    return 2.0 * Time::fromCalendar(TimeScale::UTC, year, month, day, 12).jd().jd2;
}

} // namespace astro::detail
