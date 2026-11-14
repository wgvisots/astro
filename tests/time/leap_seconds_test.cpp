/**********************************************************************
 * leap_seconds_test.cpp
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

#include "astro/time/leap_seconds.h"

#include <string>

#include "astro/core/error.h"
#include "test_util.h"

using namespace astro;

namespace {

// The three published files as they stood in November 2026: every data line, and the header lines that
// matter here (shortened where they are only prose).

// International Earth Rotation and Reference Systems Service (IERS) Leap_Second.dat,
// https://hpiers.obspm.fr/iers/bul/bulc/Leap_Second.dat
const std::string kLeapSecondDat = "#  Value of TAI-UTC in second valid beetween the initial value until\n"
                                   "#  the epoch given on the next line. The last line reads that NO\n"
                                   "#  leap second was introduced since the corresponding date\n"
                                   "#  Updated through IERS Bulletin 72 issued in July 2026\n"
                                   "#\n"
                                   "#  File expires on 28 June 2027\n"
                                   "#\n"
                                   "#    MJD        Date        TAI-UTC (s)\n"
                                   "#           day month year\n"
                                   "#    ---    --------------   ------\n"
                                   "#\n"
                                   "    41317.0    1  1 1972       10\n"
                                   "    41499.0    1  7 1972       11\n"
                                   "    41683.0    1  1 1973       12\n"
                                   "    42048.0    1  1 1974       13\n"
                                   "    42413.0    1  1 1975       14\n"
                                   "    42778.0    1  1 1976       15\n"
                                   "    43144.0    1  1 1977       16\n"
                                   "    43509.0    1  1 1978       17\n"
                                   "    43874.0    1  1 1979       18\n"
                                   "    44239.0    1  1 1980       19\n"
                                   "    44786.0    1  7 1981       20\n"
                                   "    45151.0    1  7 1982       21\n"
                                   "    45516.0    1  7 1983       22\n"
                                   "    46247.0    1  7 1985       23\n"
                                   "    47161.0    1  1 1988       24\n"
                                   "    47892.0    1  1 1990       25\n"
                                   "    48257.0    1  1 1991       26\n"
                                   "    48804.0    1  7 1992       27\n"
                                   "    49169.0    1  7 1993       28\n"
                                   "    49534.0    1  7 1994       29\n"
                                   "    50083.0    1  1 1996       30\n"
                                   "    50630.0    1  7 1997       31\n"
                                   "    51179.0    1  1 1999       32\n"
                                   "    53736.0    1  1 2006       33\n"
                                   "    54832.0    1  1 2009       34\n"
                                   "    56109.0    1  7 2012       35\n"
                                   "    57204.0    1  7 2015       36\n"
                                   "    57754.0    1  1 2017       37\n";

// IERS / Internet Assigned Numbers Authority (IANA) leap-seconds.list,
// https://data.iana.org/time-zones/data/leap-seconds.list. Its comments use words that belong to the other two
// formats ("File expires on", "TAI-UTC", "MJD"); times are seconds of the Network Time Protocol (NTP) era.
const std::string kLeapSecondsList = "#\tThe following line shows the last update of this file in NTP timestamp:\n"
                                     "#\n"
                                     "#$\t3992312697\n"
                                     "#\n"
                                     "#\t2) Expiration date of the file given on a semi-annual basis\n"
                                     "#\n"
                                     "#\tFile expires on 28 June 2027\n"
                                     "#\n"
                                     "#\tExpire date in NTP timestamp:\n"
                                     "#\n"
                                     "#@\t4023129600\n"
                                     "#\n"
                                     "#\tLIST OF LEAP SECONDS\n"
                                     "#\tNTP timestamp (X parameter) is the number of seconds since 1900.0\n"
                                     "#\n"
                                     "#\tMJD: The Modified Julian Day number. MJD = X/86400 + 15020\n"
                                     "#\n"
                                     "#\tDTAI: The difference DTAI= TAI-UTC in units of seconds\n"
                                     "#\tIt is the quantity to add to UTC to get the time in TAI\n"
                                     "#\n"
                                     "#NTP Time      DTAI    Day Month Year\n"
                                     "#\n"
                                     "2272060800      10      # 1 Jan 1972\n"
                                     "2287785600      11      # 1 Jul 1972\n"
                                     "2303683200      12      # 1 Jan 1973\n"
                                     "2335219200      13      # 1 Jan 1974\n"
                                     "2366755200      14      # 1 Jan 1975\n"
                                     "2398291200      15      # 1 Jan 1976\n"
                                     "2429913600      16      # 1 Jan 1977\n"
                                     "2461449600      17      # 1 Jan 1978\n"
                                     "2492985600      18      # 1 Jan 1979\n"
                                     "2524521600      19      # 1 Jan 1980\n"
                                     "2571782400      20      # 1 Jul 1981\n"
                                     "2603318400      21      # 1 Jul 1982\n"
                                     "2634854400      22      # 1 Jul 1983\n"
                                     "2698012800      23      # 1 Jul 1985\n"
                                     "2776982400      24      # 1 Jan 1988\n"
                                     "2840140800      25      # 1 Jan 1990\n"
                                     "2871676800      26      # 1 Jan 1991\n"
                                     "2918937600      27      # 1 Jul 1992\n"
                                     "2950473600      28      # 1 Jul 1993\n"
                                     "2982009600      29      # 1 Jul 1994\n"
                                     "3029443200      30      # 1 Jan 1996\n"
                                     "3076704000      31      # 1 Jul 1997\n"
                                     "3124137600      32      # 1 Jan 1999\n"
                                     "3345062400      33      # 1 Jan 2006\n"
                                     "3439756800      34      # 1 Jan 2009\n"
                                     "3550089600      35      # 1 Jul 2012\n"
                                     "3644697600      36      # 1 Jul 2015\n"
                                     "3692217600      37      # 1 Jan 2017\n"
                                     "#\n"
                                     "#h\ta9bad145 84c31c70 758402aa b37bfd54 5923836a\n";

// United States Naval Observatory (USNO) tai-utc.dat. It has no header, and it starts in 1961.
const std::string kTaiUtcDat = " 1961 JAN  1 =JD 2437300.5  TAI-UTC=   1.4228180 S + (MJD - 37300.) X 0.001296 S\n"
                               " 1961 AUG  1 =JD 2437512.5  TAI-UTC=   1.3728180 S + (MJD - 37300.) X 0.001296 S\n"
                               " 1962 JAN  1 =JD 2437665.5  TAI-UTC=   1.8458580 S + (MJD - 37665.) X 0.0011232S\n"
                               " 1963 NOV  1 =JD 2438334.5  TAI-UTC=   1.9458580 S + (MJD - 37665.) X 0.0011232S\n"
                               " 1964 JAN  1 =JD 2438395.5  TAI-UTC=   3.2401300 S + (MJD - 38761.) X 0.001296 S\n"
                               " 1964 APR  1 =JD 2438486.5  TAI-UTC=   3.3401300 S + (MJD - 38761.) X 0.001296 S\n"
                               " 1964 SEP  1 =JD 2438639.5  TAI-UTC=   3.4401300 S + (MJD - 38761.) X 0.001296 S\n"
                               " 1965 JAN  1 =JD 2438761.5  TAI-UTC=   3.5401300 S + (MJD - 38761.) X 0.001296 S\n"
                               " 1965 MAR  1 =JD 2438820.5  TAI-UTC=   3.6401300 S + (MJD - 38761.) X 0.001296 S\n"
                               " 1965 JUL  1 =JD 2438942.5  TAI-UTC=   3.7401300 S + (MJD - 38761.) X 0.001296 S\n"
                               " 1965 SEP  1 =JD 2439004.5  TAI-UTC=   3.8401300 S + (MJD - 38761.) X 0.001296 S\n"
                               " 1966 JAN  1 =JD 2439126.5  TAI-UTC=   4.3131700 S + (MJD - 39126.) X 0.002592 S\n"
                               " 1968 FEB  1 =JD 2439887.5  TAI-UTC=   4.2131700 S + (MJD - 39126.) X 0.002592 S\n"
                               " 1972 JAN  1 =JD 2441317.5  TAI-UTC=  10.0       S + (MJD - 41317.) X 0.0      S\n"
                               " 1972 JUL  1 =JD 2441499.5  TAI-UTC=  11.0       S + (MJD - 41317.) X 0.0      S\n"
                               " 1973 JAN  1 =JD 2441683.5  TAI-UTC=  12.0       S + (MJD - 41317.) X 0.0      S\n"
                               " 1974 JAN  1 =JD 2442048.5  TAI-UTC=  13.0       S + (MJD - 41317.) X 0.0      S\n"
                               " 1975 JAN  1 =JD 2442413.5  TAI-UTC=  14.0       S + (MJD - 41317.) X 0.0      S\n"
                               " 1976 JAN  1 =JD 2442778.5  TAI-UTC=  15.0       S + (MJD - 41317.) X 0.0      S\n"
                               " 1977 JAN  1 =JD 2443144.5  TAI-UTC=  16.0       S + (MJD - 41317.) X 0.0      S\n"
                               " 1978 JAN  1 =JD 2443509.5  TAI-UTC=  17.0       S + (MJD - 41317.) X 0.0      S\n"
                               " 1979 JAN  1 =JD 2443874.5  TAI-UTC=  18.0       S + (MJD - 41317.) X 0.0      S\n"
                               " 1980 JAN  1 =JD 2444239.5  TAI-UTC=  19.0       S + (MJD - 41317.) X 0.0      S\n"
                               " 1981 JUL  1 =JD 2444786.5  TAI-UTC=  20.0       S + (MJD - 41317.) X 0.0      S\n"
                               " 1982 JUL  1 =JD 2445151.5  TAI-UTC=  21.0       S + (MJD - 41317.) X 0.0      S\n"
                               " 1983 JUL  1 =JD 2445516.5  TAI-UTC=  22.0       S + (MJD - 41317.) X 0.0      S\n"
                               " 1985 JUL  1 =JD 2446247.5  TAI-UTC=  23.0       S + (MJD - 41317.) X 0.0      S\n"
                               " 1988 JAN  1 =JD 2447161.5  TAI-UTC=  24.0       S + (MJD - 41317.) X 0.0      S\n"
                               " 1990 JAN  1 =JD 2447892.5  TAI-UTC=  25.0       S + (MJD - 41317.) X 0.0      S\n"
                               " 1991 JAN  1 =JD 2448257.5  TAI-UTC=  26.0       S + (MJD - 41317.) X 0.0      S\n"
                               " 1992 JUL  1 =JD 2448804.5  TAI-UTC=  27.0       S + (MJD - 41317.) X 0.0      S\n"
                               " 1993 JUL  1 =JD 2449169.5  TAI-UTC=  28.0       S + (MJD - 41317.) X 0.0      S\n"
                               " 1994 JUL  1 =JD 2449534.5  TAI-UTC=  29.0       S + (MJD - 41317.) X 0.0      S\n"
                               " 1996 JAN  1 =JD 2450083.5  TAI-UTC=  30.0       S + (MJD - 41317.) X 0.0      S\n"
                               " 1997 JUL  1 =JD 2450630.5  TAI-UTC=  31.0       S + (MJD - 41317.) X 0.0      S\n"
                               " 1999 JAN  1 =JD 2451179.5  TAI-UTC=  32.0       S + (MJD - 41317.) X 0.0      S\n"
                               " 2006 JAN  1 =JD 2453736.5  TAI-UTC=  33.0       S + (MJD - 41317.) X 0.0      S\n"
                               " 2009 JAN  1 =JD 2454832.5  TAI-UTC=  34.0       S + (MJD - 41317.) X 0.0      S\n"
                               " 2012 JUL  1 =JD 2456109.5  TAI-UTC=  35.0       S + (MJD - 41317.) X 0.0      S\n"
                               " 2015 JUL  1 =JD 2457204.5  TAI-UTC=  36.0       S + (MJD - 41317.) X 0.0      S\n"
                               " 2017 JAN  1 =JD 2457754.5  TAI-UTC=  37.0       S + (MJD - 41317.) X 0.0      S\n";

// Good first lines in each format; the tests below add one damaged line to them.
const std::string kDatHead  = "    41317.0    1  1 1972       10\n"
                              "    57204.0    1  7 2015       36\n";
const std::string kListHead = "#@ 4023129600\n"
                              "2272060800 10\n"
                              "3644697600 36\n";
const std::string kUsnoHead = " 1972 JAN  1 =JD 2441317.5  TAI-UTC=  10.0       S + (MJD - 41317.) X 0.0      S\n"
                              " 2015 JUL  1 =JD 2457204.5  TAI-UTC=  36.0       S + (MJD - 41317.) X 0.0      S\n";

// The same text with Windows line endings.
std::string withCrLf(const std::string& text) {
    std::string out;
    for (const char ch : text) {
        if (ch == '\n') {
            out += '\r';
        }
        out += ch;
    }
    return out;
}

// Parsing `content` must throw astro::Error, and the message must contain `what`.
void expectParseError(const std::string& content, const std::string& what) {
    try {
        const LeapSecondTable t = LeapSecondTable::parse(content);
        ADD_FAILURE() << "parsed without error (TAI-UTC on 2024-01-01 = " << t.taiMinusUtc(2024, 1, 1) << " s):\n"
                      << content;
    } catch (const Error& e) {
        EXPECT_NE(std::string(e.what()).find(what), std::string::npos)
            << "expected \"" << what << "\" in the message, which was: " << e.what();
    }
}

} // namespace

TEST(LeapSecondFile, AllThreeFormatsGiveTheBuiltInTable) {
    // The built-in table holds the same history (it was checked against these files and against the Essential
    // Routines for Fundamental Astronomy (ERFA) function eraDat), so each file must reproduce it entry by entry.
    const LeapSecondTable builtin = LeapSecondTable::builtin();
    for (const std::string* file : {&kLeapSecondDat, &kLeapSecondsList, &kTaiUtcDat}) {
        for (const bool crlf : {false, true}) {
            const LeapSecondTable t = LeapSecondTable::parse(crlf ? withCrLf(*file) : *file);
            ASSERT_EQ(t.entries().size(), builtin.entries().size()) << file->substr(0, 40);
            for (std::size_t i = 0; i < t.entries().size(); ++i) {
                const LeapSecondEntry& a = t.entries()[i];
                const LeapSecondEntry& b = builtin.entries()[i];
                EXPECT_EQ(a.year, b.year) << i;
                EXPECT_EQ(a.month, b.month) << i;
                EXPECT_EQ(a.offset, b.offset) << i;
                if (a.rate != 0.0 || b.rate != 0.0) { // the reference MJD only matters with a rate
                    EXPECT_EQ(a.mjdRef, b.mjdRef) << i;
                    EXPECT_EQ(a.rate, b.rate) << i;
                }
            }
            // Published values: IERS Bulletin C (37 s from 2017-01-01, 36 s before it, 10 s at the start of 1972).
            EXPECT_EQ(t.taiMinusUtc(2026, 11, 12), 37.0);
            EXPECT_EQ(t.taiMinusUtc(2017, 1, 1), 37.0);
            EXPECT_EQ(t.taiMinusUtc(2016, 12, 31), 36.0);
            EXPECT_EQ(t.taiMinusUtc(1972, 1, 1), 10.0);
        }
    }
}

TEST(LeapSecondFile, KeepsTheBuiltInHistoryBeforeTheFilesFirstEntry) {
    // tai-utc.dat starts on 1961-01-01. The year 1960 comes from the built-in table (ERFA eraDat: 1.4178180 s on
    // 1960-01-01, drifting 0.001296 s/day about MJD 37300); 1960-06-01 is MJD 37086.
    const LeapSecondTable usno = LeapSecondTable::parse(kTaiUtcDat);
    EXPECT_EQ(usno.entries().front().year, 1960);
    EXPECT_NEAR(usno.taiMinusUtc(1960, 6, 1), 1.4178180 + (37086.0 - 37300.0) * 0.001296, 1e-12);
    // From its own first line on, the file is used: 1961-01-01 is MJD 37300, 1965-03-01 12h is MJD 38820.5.
    EXPECT_NEAR(usno.taiMinusUtc(1961, 1, 1), 1.4228180, 1e-12);
    EXPECT_NEAR(usno.taiMinusUtc(1965, 3, 1, 0.5), 3.6401300 + (38820.5 - 38761.0) * 0.001296, 1e-12);

    // A file cut down to its last lines still has the earlier steps (19 s through 1980, Bulletin C).
    const LeapSecondTable recent = LeapSecondTable::parse("    57204.0    1  7 2015       36\n"
                                                          "    57754.0    1  1 2017       37\n");
    EXPECT_EQ(recent.taiMinusUtc(1980, 6, 1), 19.0);
    EXPECT_EQ(recent.taiMinusUtc(2015, 6, 30), 35.0);
    EXPECT_EQ(recent.taiMinusUtc(2015, 7, 1), 36.0);
    EXPECT_EQ(recent.taiMinusUtc(2024, 1, 1), 37.0);

    // No entry of a file is dropped, whatever its date: here a (made up) step on 1970-01-01, NTP 2208988800.
    const LeapSecondTable early = LeapSecondTable::parse("2208988800 5\n2272060800 10\n");
    EXPECT_EQ(early.taiMinusUtc(1970, 6, 1), 5.0);
    EXPECT_EQ(early.taiMinusUtc(1972, 1, 1), 10.0);
    EXPECT_NEAR(early.taiMinusUtc(1960, 6, 1), 1.4178180 + (37086.0 - 37300.0) * 0.001296, 1e-12);
}

TEST(LeapSecondFile, FormatIsRecognisedFromTheDataNotTheComments) {
    // Comments that use the words of another format must not change how the data lines are read.
    for (const char* comment : {"# times here are MJD, not NTP seconds\n",
                                "# see the #@ line of leap-seconds.list for the expiry\n",
                                "# TAI-UTC= value in seconds\n"}) {
        const LeapSecondTable t = LeapSecondTable::parse(comment + kLeapSecondDat);
        EXPECT_EQ(t.taiMinusUtc(2024, 1, 1), 37.0) << comment;
        EXPECT_EQ(t.taiMinusUtc(1990, 6, 1), 25.0) << comment;
    }
    EXPECT_EQ(LeapSecondTable::parse("# TAI-UTC= value in seconds\n" + kLeapSecondsList).taiMinusUtc(2024, 1, 1), 37.0);

    // A leap-seconds.list without its "#@" line and without the word NTP is still a leap-seconds.list.
    const LeapSecondTable bare = LeapSecondTable::parse("2272060800 10\n3692217600 37\n");
    EXPECT_EQ(bare.taiMinusUtc(2024, 1, 1), 37.0);
    EXPECT_EQ(bare.taiMinusUtc(1980, 1, 1), 10.0);
    EXPECT_FALSE(bare.expiresMjd().has_value());

    // Blank lines and indented comments are not data.
    EXPECT_EQ(LeapSecondTable::parse("\n   \n  # indented comment\n" + kDatHead).taiMinusUtc(2024, 1, 1), 36.0);

    // Text in none of the three formats is refused, with the line that could not be recognised.
    expectParseError("<html><body>404 Not Found</body></html>\n", "line 1");
    expectParseError("# only a comment\n\n", "no entries");
    expectParseError("", "no entries");
}

TEST(LeapSecondFile, MalformedDataLineThrowsAndNamesTheLine) {
    // Leap_Second.dat: "MJD day month year TAI-UTC". Skipping the damaged third line instead would leave
    // TAI-UTC at 36 s after 2017, one second wrong, with no error.
    expectParseError(kDatHead + "    57754.0    1  1 20", "line 3");                // cut off
    expectParseError(kDatHead + "    57754.0    1  1 2017\n", "line 3");            // no value
    expectParseError(kDatHead + "    57754.O    1  1 2017       37\n", "line 3");   // letter O for zero
    expectParseError(kDatHead + "    57754.0    1  1 2017       37 x\n", "line 3"); // text after the value
    expectParseError(kDatHead + "    57754.0    1  1 2017       3 7\n", "line 3");  // six fields
    expectParseError(kDatHead + "    57754.0    1 13 2017       37\n", "line 3");   // month 13
    expectParseError(kDatHead + "    57754.0    1  1 2071       37\n", "line 3");   // the MJD is another date
    expectParseError(kDatHead + "    62697.0   15  7 2030       38\n", "line 3");   // not the first of a month
    expectParseError("# header\n\n" + kDatHead + "garbage\n", "line 5");            // numbering counts every line
    expectParseError("    41317.0    1  1 1972       10\nhello world\n", "'hello world'"); // the line is quoted

    // leap-seconds.list: "<NTP seconds> <TAI-UTC>", optionally followed by a comment.
    expectParseError(kListHead + "x3692217600 37\n", "line 4");
    expectParseError(kListHead + "3692217600\n", "line 4");       // no value
    expectParseError(kListHead + "3692217600 37 38\n", "line 4"); // three fields
    expectParseError(kListHead + "3692217600 3x # 1 Jan 2017\n", "line 4");
    expectParseError(kListHead + "3692260800 37\n", "line 4"); // 12h, not 0h, on 2017-01-01
    expectParseError(kListHead + "3693427200 37\n", "line 4"); // 2017-01-15
    expectParseError(kListHead + "1e300 37\n", "line 4");      // not a time
    expectParseError("#@ soon\n2272060800 10\n", "line 1");    // expiry is not a number

    // tai-utc.dat.
    expectParseError(kUsnoHead + " 2017 JAN  1 =JD 2457754.5  TAI-UTC\n", "line 3"); // cut off
    expectParseError(kUsnoHead + " 2017 JAN  1 =JD 2457754.5  TAI-UTC=  3O.0       S + (MJD - 41317.) X 0.0      S\n",
                     "line 3"); // letter O for zero
    expectParseError(kUsnoHead + " 2017 JNA  1 =JD 2457754.5  TAI-UTC=  37.0       S + (MJD - 41317.) X 0.0      S\n",
                     "line 3"); // no such month
    expectParseError(kUsnoHead + " 2017 JAN 15 =JD 2457768.5  TAI-UTC=  37.0       S + (MJD - 41317.) X 0.0      S\n",
                     "line 3"); // not the first of a month
    expectParseError(kUsnoHead + " 2071 JAN  1 =JD 2457754.5  TAI-UTC=  37.0       S + (MJD - 41317.) X 0.0      S\n",
                     "line 3"); // the JD is another date
    expectParseError(kUsnoHead + " 2017 JAN  1 =JD 2457754.5  TAI-UTC=  37.0       S\n", "line 3"); // no rate term
}

TEST(LeapSecondFile, ExpiryDateIsTheOneInTheFile) {
    // Both IERS files expire on 2027-06-28, which is MJD 61584 (MJD 51544 is 2000-01-01; 27 years with 7 leap
    // days, then 178 days to 28 June). In leap-seconds.list it is NTP second 4023129600 = 86400 * (61584 - 15020).
    const LeapSecondTable dat = LeapSecondTable::parse(kLeapSecondDat);
    ASSERT_TRUE(dat.expiresMjd().has_value());
    EXPECT_EQ(*dat.expiresMjd(), 61584.0);
    const LeapSecondTable list = LeapSecondTable::parse(kLeapSecondsList);
    ASSERT_TRUE(list.expiresMjd().has_value());
    EXPECT_EQ(*list.expiresMjd(), 61584.0);
    // tai-utc.dat and the built-in table do not state one.
    EXPECT_FALSE(LeapSecondTable::parse(kTaiUtcDat).expiresMjd().has_value());
    EXPECT_FALSE(LeapSecondTable::builtin().expiresMjd().has_value());
    // A table built from entries reports the date it was given.
    EXPECT_EQ(LeapSecondTable(
                  {
                      {1972, 1, 10.0}
    },
                  61584.0)
                  .expiresMjd(),
              61584.0);
}
