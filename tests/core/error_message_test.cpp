/**********************************************************************
 * error_message_test.cpp
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

// The message format of core/error.h: "Function: what is wrong (got <value>)". A sample of failing calls
// from several modules, each of which must name its function first. It guards the format, not every message.

#include <array>
#include <functional>
#include <string>

#include <gtest/gtest.h>

#include "astro/core/angle.h"
#include "astro/core/error.h"
#include "astro/core/mat6.h"
#include "astro/ephemeris/low_precision.h"
#include "astro/frames/transform.h"
#include "astro/orbit/elements.h"
#include "astro/orbit/sgp4.h"
#include "astro/orbit/tle.h"
#include "astro/time/julian_date.h"
#include "astro/time/leap_seconds.h"
#include "astro/time/time.h"

using namespace astro;

namespace {

// The message of the astro::Error that `call` throws; empty if it throws none.
std::string messageOf(const std::function<void()>& call) {
    try {
        call();
    } catch (const Error& e) {
        return e.what();
    }
    return {};
}

} // namespace

TEST(ErrorMessage, StartsWithTheFunctionName) {
    const Time t   = Time::fromCalendar(TimeScale::UTC, 2024, 3, 1);
    Tle        tle = Tle::parse("1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927",
                         "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537");
    tle.meanMotion = -2.0;
    struct Case {
        const char*           function;
        std::function<void()> call;
    };
    const std::array<double, 6> d{};
    const Case                  cases[] = {
        {"gcrsToCelestialMatrix",        [&] { (void)gcrsToCelestialMatrix(static_cast<Frame>(99), t); }          },
        {"poleOfDateMatrix",             [&] { (void)poleOfDateMatrix(Frame::ITRS, t); }                          },
        {"dayOfYearToDate",              [] { (void)dayOfYearToDate(2023, 366); }                                 },
        {"daysInMonth",                  [] { (void)daysInMonth(2023, 13); }                                      },
        {"parseHms",                     [] { (void)parseHms("12:3x"); }                                          },
        {"parseDms",                     [] { (void)parseDms("1:2:3:4"); }                                        },
        {"parseTimeScale",               [] { (void)parseTimeScale("XYZ"); }                                      },
        {"Time::fromIso8601",            [] { (void)Time::fromIso8601("not a time"); }                            },
        {"Time::to",                     [&] { (void)Time(t.jd(), TimeScale::UT1).to(TimeScale::TT); }            },
        {"Time::to",                     [&] { (void)Time(t.jd(), static_cast<TimeScale>(42)).to(TimeScale::TT); }},
        {"LeapSecondTable::parse",       [] { (void)LeapSecondTable::parse("not a leap-second file\n"); }         },
        {"LeapSecondTable::taiMinusUtc", [] { (void)LeapSecondTable::builtin().taiMinusUtc(1e10); }               },
        {"LowPrecisionEphemeris::state", [&] { (void)LowPrecisionEphemeris().state(Body::Mars, Body::Earth, t); } },
        {"keplerianToEquinoctial",       [] { (void)keplerianToEquinoctial({-2.0e7, 1.5, 0.5, 1.0, 2.0, 1.5}); }  },
        {"equinoctialToKeplerian",       [] { (void)equinoctialToKeplerian({-2.0e7, 1.2, 0.9, 0.1, 0.2, 1.0}); }  },
        {"mahalanobis2",                 [&] { (void)mahalanobis2(d, Mat6::identity(), 7); }                      },
        {"Sgp4",                         [&] { (void)Sgp4(tle); }                                                 },
    };
    for (const Case& c : cases) {
        const std::string message = messageOf(c.call);
        EXPECT_EQ(message.rfind(std::string(c.function) + ": ", 0), 0u)
            << "expected \"" << c.function << ": ...\", got \"" << message << "\"";
    }
}

TEST(ErrorMessage, GivesTheValueItGot) {
    // "(got <value>)" ends the message where there is one value to show.
    EXPECT_EQ(messageOf([] { (void)dayOfYearToDate(2023, 366); }),
              "dayOfYearToDate: day of year out of range (got 366)");
    const Time t = Time::fromCalendar(TimeScale::UTC, 2024, 3, 1);
    EXPECT_EQ(messageOf([&] { (void)gcrsToCelestialMatrix(static_cast<Frame>(99), t); }),
              "gcrsToCelestialMatrix: invalid Frame value (got 99)");
    EXPECT_EQ(messageOf([] { (void)parseTimeScale("XYZ"); }), "parseTimeScale: unknown time scale (got 'XYZ')");
}
