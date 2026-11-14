/**********************************************************************
 * frames_test.cpp
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

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

#include "astro/core/angle.h"
#include "astro/core/error.h"
#include "astro/core/rotation.h"
#include "astro/frames/earth_rotation.h"
#include "astro/frames/eop.h"
#include "astro/frames/precession_nutation.h"
#include "astro/frames/transform.h"
#include "test_util.h"

using namespace astro;
using astro::test::dataPath;
using astro::test::expectMatNear;
using astro::test::expectOrthonormal;
using astro::test::expectVecNear;
using astro::test::readCsv;

namespace {
constexpr double kUas = constants::uasToRad;

Mat3 matFrom(const test::Csv& csv, const std::vector<std::string>& r, const std::string& prefix) {
    Mat3 m;
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            m(i, j) = std::stod(r[csv.col(prefix + std::to_string(i) + std::to_string(j))]);
        }
    }
    return m;
}
} // namespace

TEST(PrecessionNutation, MatchesErfa) {
    const auto csv = readCsv(dataPath("frames.csv"));
    ASSERT_GE(csv.rows.size(), 50u);
    for (const auto& r : csv.rows) {
        auto             v = [&](const char* n) { return std::stod(r[csv.col(n)]); };
        const JulianDate tt{v("tt1"), v("tt2")};
        const JulianDate ut1{v("tt1"), v("ut12")};
        const double     t  = tt.centuriesSinceJ2000();
        const CipXY      xy = cipXY06(t);
        EXPECT_NEAR(xy.x, v("x"), 0.5 * kUas) << "X at t=" << t;
        EXPECT_NEAR(xy.y, v("y"), 0.5 * kUas) << "Y at t=" << t;
        EXPECT_NEAR(cioLocatorS06(t, {v("x"), v("y")}), v("s"), 0.05 * kUas) << "s";
        EXPECT_NEAR(tioLocatorSp(t), v("sp"), 1e-18);
        EXPECT_NEAR(earthRotationAngle(ut1), v("era"), 1e-12);
        EXPECT_NEAR(gmst06(ut1, tt), v("gmst06"), 1e-12);
        EXPECT_NEAR(gmst82(ut1), v("gmst82"), 1e-12);
        // The IERS R06 nutation table and SOFA's nut06a apply the IAU 2006 adjustment differently: they agree
        // to ~0.1 uas near J2000 and drift apart by a few uas at +/-100 years (GAST inherits this via EE).
        EXPECT_NEAR(wrapPi(gast06(ut1, tt) - v("gst06a")), 0.0, 5.0 * kUas);
        EXPECT_NEAR(meanObliquity06(t), v("obl06"), 1e-15);
        const Nutation na = nutationIau2000A(t);
        EXPECT_NEAR(na.dpsi, v("dpsi06a"), 5.0 * kUas) << "dpsi";
        EXPECT_NEAR(na.deps, v("deps06a"), 1.0 * kUas) << "deps";
        const Nutation nb = nutationIau2000B(t);
        EXPECT_NEAR(nb.dpsi, v("dpsi00b"), 1e-15);
        EXPECT_NEAR(nb.deps, v("deps00b"), 1e-15);
        // 2000B is a ~1 mas approximation of 2000A (over 1995-2050; looser over 1900-2100).
        EXPECT_NEAR(nb.dpsi, na.dpsi, 5.0 * constants::masToRad);

        expectMatNear(frameBiasMatrix06(), matFrom(csv, r, "rb"), 1e-15, "bias");
        expectMatNear(precessionMatrix06(t), matFrom(csv, r, "rp"), 1e-14, "precession");
        expectMatNear(biasPrecessionMatrix06(t), matFrom(csv, r, "rbp"), 1e-14, "bias-precession");
        expectMatNear(npbMatrix06(t), matFrom(csv, r, "npb"), 5.0 * kUas, "npb");

        // Full CIO-based GCRS -> ITRS with polar motion and pole offsets.
        const CipXY xyc{xy.x + v("dx"), xy.y + v("dy")};
        const Mat3  c2t = polarMotionMatrix(v("xp"), v("yp"), tioLocatorSp(t)) * R3(earthRotationAngle(ut1))
                       * celestialToIntermediateMatrix(xyc, cioLocatorS06(t, xy));
        expectMatNear(c2t, matFrom(csv, r, "c2t"), 1.0 * kUas, "c2t");
        expectOrthonormal(c2t);
    }
}

TEST(PrecessionNutation, PoleOffsetsToNutationRoundTrip) {
    const double   t  = 0.2;
    const double   dX = 0.3 * constants::masToRad, dY = -0.2 * constants::masToRad;
    const Nutation n = poleOffsetsToNutation(t, dX, dY);
    // Applying the nutation correction should move the CIP by (dX, dY).
    const CipXY a = cipXYFromMatrix(npbMatrix06(t));
    const CipXY b = cipXYFromMatrix(npbMatrix06(t, PrecessionNutationModel::IAU2006_2000A, n.dpsi, n.deps));
    EXPECT_NEAR(b.x - a.x, dX, 1e-3 * constants::masToRad);
    EXPECT_NEAR(b.y - a.y, dY, 1e-3 * constants::masToRad);
}

// Vallado, Fundamentals of Astrodynamics and Applications (4th ed.), Example 3-15,
// and Vallado et al. (2006) "Revisiting Spacetrack Report #3", Appendix C.
class ValladoExample : public ::testing::Test {
protected:
    const Time        utc = Time::fromCalendar(TimeScale::UTC, 2004, 4, 6, 7, 51, 28.386009);
    const ConstantEop eop{
        EopValues{.xp   = arcsec2rad(-0.140682),
                  .yp   = arcsec2rad(0.333309),
                  .dut1 = -0.4399619,
                  .lod  = 0.0015563,
                  .dX   = arcsec2rad(-0.000205),
                  .dY   = arcsec2rad(-0.000136)}
    };
    const StateVector itrf{
        {-1033.4793830e3, 7901.2952754e3, 6380.3565958e3},
        {-3.225636520e3,  -2.872451450e3, 5.531924446e3 }
    };
};

TEST_F(ValladoExample, ItrfToGcrf) {
    const StateVector gcrf = itrsToGcrs(itrf, utc, eop);
    // Vallado's values use IAU 2000; IAU 2006 differs at the centimeter level.
    expectVecNear(gcrf.r, {5102.508959e3, 6123.011403e3, 6378.136925e3}, 0.05, "r GCRF");
    expectVecNear(gcrf.v, {-4.743220156e3, 0.790536497e3, 5.533755724e3}, 5e-5, "v GCRF");
    const StateVector back = gcrsToItrs(gcrf, utc, eop);
    expectVecNear(back.r, itrf.r, 1e-6, "round trip r");
    expectVecNear(back.v, itrf.v, 1e-9, "round trip v");
}

TEST_F(ValladoExample, TemeAndPef) {
    const FrameTransformer ft(utc, eop);
    const StateVector      teme = ft.transform(itrf, Frame::ITRS, Frame::TEME);
    expectVecNear(teme.r, {5094.18016210e3, 6127.64465950e3, 6380.34453270e3}, 0.05, "r TEME");
    expectVecNear(teme.v, {-4.746131487e3, 0.785818041e3, 5.531931288e3}, 5e-5, "v TEME");
    const StateVector pef = ft.transform(itrf, Frame::ITRS, Frame::PEF);
    expectVecNear(pef.r, {-1033.47503130e3, 7901.30558560e3, 6380.34453270e3}, 0.05, "r PEF");
    expectVecNear(pef.v, {-3.2256327470e3, -2.8724425110e3, 5.5319312880e3}, 5e-5, "v PEF");
    // TIRS and PEF agree closely for IAU 2006.
    expectVecNear(ft.transform(itrf.r, Frame::ITRS, Frame::TIRS), pef.r, 0.01, "TIRS vs PEF");
    // J2000 (EME2000) differs from GCRS by the frame bias (~23 mas -> ~0.7 m here).
    const Vec3 j2000 = ft.transform(itrf.r, Frame::ITRS, Frame::J2000);
    const Vec3 gcrs  = ft.transform(itrf.r, Frame::ITRS, Frame::GCRS);
    EXPECT_GT(norm(j2000 - gcrs), 0.3);
    EXPECT_LT(norm(j2000 - gcrs), 1.5);
}

TEST(FrameTransformer, ConsistencyAndRoundTrips) {
    const Time             t = Time::fromIso8601("2021-07-04T12:34:56Z");
    const ConstantEop      eop(EopValues{.xp = arcsec2rad(0.2), .yp = arcsec2rad(0.4), .dut1 = -0.15, .lod = 0.001});
    const FrameTransformer ft(t, eop);
    const StateVector      s{
             {7000e3, -1200e3, 300e3},
             {1.0e3,  7.2e3,   0.5e3}
    };
    const Frame all[] = {Frame::GCRS,
                         Frame::J2000,
                         Frame::MOD,
                         Frame::TOD,
                         Frame::TEME,
                         Frame::CIRS,
                         Frame::TIRS,
                         Frame::PEF,
                         Frame::ITRS};
    for (Frame a : all) {
        expectOrthonormal(ft.gcrsTo(a), 1e-14);
        for (Frame b : all) {
            const StateVector x = ft.transform(s, Frame::GCRS, a);
            const StateVector y = ft.transform(ft.transform(x, a, b), b, Frame::GCRS);
            expectVecNear(y.r, s.r, 1e-6);
            expectVecNear(y.v, s.v, 1e-9);
        }
    }
    // An ITRS-fixed point has inertial velocity omega x r.
    const StateVector fixed{
        {6378137.0, 0, 0},
        {}
    };
    const StateVector inertial = ft.transform(fixed, Frame::ITRS, Frame::GCRS);
    EXPECT_NEAR(norm(inertial.v), kEarthRotationRateEra * 6378137.0, 0.01);
    // The CIO and equinox routes agree (TIRS vs PEF) to well under a milliarcsecond.
    const Mat3 d = ft.gcrsTo(Frame::TIRS) * ft.gcrsTo(Frame::PEF).transpose();
    EXPECT_LT(std::abs(d(0, 1)), 0.1 * constants::masToRad);
    // GAST - GMST is the equation of the equinoxes (< 1.2 s of time).
    EXPECT_LT(std::abs(wrapPi(ft.gast() - ft.gmst())), deg2rad(1.2 * 15.0 / 3600.0));
    // Without Earth-orientation data (ConstantEop{}: UT1 taken as UTC, no polar motion) the transformer
    // is still built.
    EXPECT_NO_THROW(FrameTransformer(t, ConstantEop{}));
}

TEST(FrameTransformer, Model2000BWithinMilliarcsecond) {
    const Time        t = Time::fromIso8601("2030-01-01T00:00:00Z");
    const ConstantEop eop;
    const Mat3        a     = gcrsToItrsMatrix(t, eop);
    const Mat3        b     = gcrsToItrsMatrix(t, eop, {.model = PrecessionNutationModel::IAU2006_2000B});
    const Mat3        d     = a * b.transpose();
    const double      angle = std::sqrt(d(0, 1) * d(0, 1) + d(0, 2) * d(0, 2) + d(1, 2) * d(1, 2));
    EXPECT_LT(angle, 2.0 * constants::masToRad);
}

namespace {
std::string finalsLine(int    y,
                       int    m,
                       int    d,
                       double mjd,
                       double xp,
                       double yp,
                       double dut1,
                       double lodMs,
                       double dXmas,
                       double dYmas) {
    char buf[200];
    std::snprintf(buf,
                  sizeof(buf),
                  "%2d%02d%02d %8.2f %c %9.6f%9.6f %9.6f%9.6f  %c%10.7f%10.7f %7.4f%7.4f  %c %9.3f%9.3f %9.3f%9.3f",
                  y % 100,
                  m,
                  d,
                  mjd,
                  'I',
                  xp,
                  0.00003,
                  yp,
                  0.00003,
                  'I',
                  dut1,
                  0.00001,
                  lodMs,
                  0.007,
                  'I',
                  dXmas,
                  0.1,
                  dYmas,
                  0.1);
    return buf;
}
} // namespace

TEST(EopTable, Finals2000AParsingAndLeapSecondInterpolation) {
    // Around the 2016-12-31 leap second: UT1-UTC jumps by +1 s.
    std::string content;
    content += finalsLine(2016, 12, 29, 57751.0, 0.030, 0.280, -0.4050000, 1.0, 0.10, -0.10) + "\n";
    content += finalsLine(2016, 12, 30, 57752.0, 0.029, 0.279, -0.4060000, 1.0, 0.11, -0.11) + "\n";
    content += finalsLine(2016, 12, 31, 57753.0, 0.028, 0.278, -0.4070000, 1.0, 0.12, -0.12) + "\n";
    content += finalsLine(2017, 1, 1, 57754.0, 0.027, 0.277, 0.5920000, 1.0, 0.13, -0.13) + "\n";
    content += finalsLine(2017, 1, 2, 57755.0, 0.026, 0.276, 0.5910000, 1.0, 0.14, -0.14) + "\n";
    content += "170103 57756.00 P  \n"; // prediction rows without data are skipped
    for (auto interp : {EopInterpolation::Linear, EopInterpolation::Lagrange4}) {
        const EopTable table = EopTable::parseFinals2000A(content, interp);
        ASSERT_EQ(table.records().size(), 5u);
        const EopValues a = table.at(57752.0);
        EXPECT_NEAR(rad2arcsec(a.xp), 0.029, 1e-12);
        EXPECT_NEAR(a.dut1, -0.406, 1e-12);
        EXPECT_NEAR(a.lod, 0.001, 1e-15);
        EXPECT_NEAR(rad2mas(a.dX), 0.11, 1e-12);
        // Mid-day before the leap: smooth, no 0.5 s artifact.
        EXPECT_NEAR(table.at(57753.5).dut1, -0.4075, 1e-4);
        // After the leap.
        EXPECT_NEAR(table.at(57754.5).dut1, 0.5915, 1e-4);
        EXPECT_THROW((void)table.at(57760.0), Error);
    }
}

TEST(EopTable, C04Formats) {
    const char*    c04_14 = "# header\n"
                            "2020   1   1  58849   0.076614   0.282383  -0.1772561   0.0003791   0.000098  -0.000233\n"
                            "2020   1   2  58850   0.074996   0.282597  -0.1775922   0.0002929   0.000099  -0.000227\n";
    const EopTable a      = EopTable::parseC04(c04_14);
    EXPECT_NEAR(rad2arcsec(a.at(58849.0).xp), 0.076614, 1e-12);
    EXPECT_NEAR(a.at(58849.0).dut1, -0.1772561, 1e-12);
    EXPECT_NEAR(a.at(58849.0).lod, 0.0003791, 1e-15);
    const char* c04_20 =
        "# YR MM DD HH MJD x y UT1-UTC dX dY xrt yrt LOD\n"
        "2020  1  1  0 58849.00   0.076614   0.282383  -0.1772561   0.000098  -0.000233   0.0  0.0   0.0003791\n"
        "2020  1  2  0 58850.00   0.074996   0.282597  -0.1775922   0.000099  -0.000227   0.0  0.0   0.0002929\n";
    const EopTable b = EopTable::parseC04(c04_20);
    EXPECT_NEAR(rad2arcsec(b.at(58850.0).yp), 0.282597, 1e-12);
    EXPECT_NEAR(b.at(58850.0).lod, 0.0002929, 1e-15);
    EXPECT_NEAR(rad2arcsec(b.at(58850.0).dX), 0.000099, 1e-12);
}

TEST(EopTable, EpochThatIsNotFiniteIsRejected) {
    // Such a record used to reach std::sort, whose comparison must be a strict order (it is not for
    // NaN: undefined behaviour), and was refused only later, with a message about leap seconds.
    for (double bad : {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity()}) {
        EopRecord a, b, c;
        a.mjd = 59000.0;
        b.mjd = bad;
        c.mjd = 59002.0;
        try {
            const EopTable table({a, b, c});
            FAIL() << "accepted an epoch of " << bad;
        } catch (const Error& e) {
            EXPECT_NE(std::string(e.what()).find("not finite"), std::string::npos) << e.what();
        }
    }
}

TEST(EopTable, DuplicateEpochIsRejected) {
    // Two records at one epoch put a zero in the denominators of the interpolation weights: the table
    // used to be accepted and at() then returned infinity or NaN. Real IERS files have one row per
    // day, so a repeated epoch is an error in the input and is refused when the table is built.
    const auto rec = [](double mjd, double xpArcsec) {
        EopRecord r;
        r.mjd       = mjd;
        r.values.xp = arcsec2rad(xpArcsec);
        return r;
    };
    const auto message = [](const std::vector<EopRecord>& records, EopInterpolation interp) -> std::string {
        try {
            const EopTable table(records, interp);
            return "accepted; xp(59000.5) = " + std::to_string(table.at(59000.5).xp);
        } catch (const Error& e) {
            return e.what();
        }
    };
    // Records may come in any order, so the two at MJD 59001 need not be next to each other.
    const std::vector<EopRecord> repeated{rec(59001.0, 2.0),
                                          rec(59003.0, 4.0),
                                          rec(59000.0, 1.0),
                                          rec(59001.0, 2.0),
                                          rec(59002.0, 3.0)};
    for (auto interp : {EopInterpolation::Linear, EopInterpolation::Lagrange4}) {
        const std::string what = message(repeated, interp);
        EXPECT_NE(what.find("duplicate epoch"), std::string::npos) << what;
        EXPECT_NE(what.find("59001"), std::string::npos) << what; // names the epoch
    }
    // The same through the file readers.
    std::string finals;
    finals += finalsLine(2020, 5, 31, 59000.0, 0.10, 0.40, -0.2400000, 0.5, 0.10, -0.10) + "\n";
    finals += finalsLine(2020, 6, 1, 59001.0, 0.11, 0.41, -0.2405000, 0.5, 0.11, -0.11) + "\n";
    finals += finalsLine(2020, 6, 1, 59001.0, 0.11, 0.41, -0.2405000, 0.5, 0.11, -0.11) + "\n";
    finals += finalsLine(2020, 6, 2, 59002.0, 0.12, 0.42, -0.2410000, 0.5, 0.12, -0.12) + "\n";
    EXPECT_THROW((void)EopTable::parseFinals2000A(finals), Error);
    const char* c04 = "2020   5  31  59000   0.100000   0.400000  -0.2400000   0.0005000   0.000100  -0.000100\n"
                      "2020   6   1  59001   0.110000   0.410000  -0.2405000   0.0005000   0.000110  -0.000110\n"
                      "2020   6   1  59001   0.110000   0.410000  -0.2405000   0.0005000   0.000110  -0.000110\n"
                      "2020   6   2  59002   0.120000   0.420000  -0.2410000   0.0005000   0.000120  -0.000120\n";
    EXPECT_THROW((void)EopTable::parseC04(c04), Error);

    // Without the repeat, the same records in the same (unsorted) order are accepted and sorted. The
    // values lie on a straight line (xp = MJD - 58999 arcsec), which both interpolations reproduce.
    std::vector<EopRecord> distinct = repeated;
    distinct.erase(distinct.begin());
    for (auto interp : {EopInterpolation::Linear, EopInterpolation::Lagrange4}) {
        const EopTable table(distinct, interp);
        ASSERT_EQ(table.records().size(), 4u);
        EXPECT_EQ(table.firstMjd(), 59000.0);
        EXPECT_EQ(table.lastMjd(), 59003.0);
        for (double mjd : {59000.0, 59000.5, 59001.0, 59001.75, 59002.25, 59003.0}) {
            EXPECT_NEAR(rad2arcsec(table.at(mjd).xp), mjd - 58999.0, 1e-12) << mjd; // rounding only
        }
    }
}

namespace {

/// Daily EOP records from MJD 60000 to 60007 (2023, no leap second) whose values are polynomials of
/// u = MJD - 60003.4: polar motion x, UT1-UTC, LOD, dX and dY cubic, polar motion y quartic.
struct PolynomialEop {
    static constexpr double c4 = 1e-3; ///< coefficient of u^4 in y, arcsec

    static double u(double mjd) { return mjd - 60003.4; }
    static double xp(double mjd) {
        return arcsec2rad(0.1 + 0.01 * u(mjd) - 0.002 * std::pow(u(mjd), 2) + 3e-4 * std::pow(u(mjd), 3));
    }
    static double yp(double mjd) { return arcsec2rad(0.3 + 0.02 * u(mjd) + c4 * std::pow(u(mjd), 4)); }
    static double dut1(double mjd) { return -0.1 + 1e-3 * u(mjd) - 2e-4 * std::pow(u(mjd), 3); }
    static double lod(double mjd) { return 1e-3 + 1e-4 * std::pow(u(mjd), 2) - 1e-5 * std::pow(u(mjd), 3); }
    static double dX(double mjd) { return mas2rad(0.2 - 0.03 * u(mjd) + 0.004 * std::pow(u(mjd), 3)); }
    static double dY(double mjd) { return mas2rad(-0.1 + 0.01 * std::pow(u(mjd), 2) - 0.002 * std::pow(u(mjd), 3)); }

    static std::vector<EopRecord> records() {
        std::vector<EopRecord> out;
        for (int d = 0; d <= 7; ++d) {
            const double m = 60000.0 + d;
            out.push_back({
                m,
                {xp(m), yp(m), dut1(m), lod(m), dX(m), dY(m)}
            });
        }
        return out;
    }
};

} // namespace

TEST(EopTable, Lagrange4ReproducesACubicExactly) {
    // Four-point Lagrange interpolation is exact for a polynomial of degree three: at nodes, in the
    // first and last intervals (where the four points cannot straddle the epoch) and inside.
    const EopTable table(PolynomialEop::records(), EopInterpolation::Lagrange4);
    for (double mjd : {60000.0, 60000.3, 60001.5, 60003.7, 60005.0, 60006.8, 60007.0}) {
        SCOPED_TRACE(testing::Message() << "MJD " << mjd);
        const EopValues v = table.at(mjd);
        // Rounding only: a few units in the last place of each value, summed with weights below 2 (and
        // UT1-UTC passes through UT1-TAI, 37 s). Measured 1.1e-22 rad, 5.3e-15 s, 4.3e-19 s, 2.1e-25 rad.
        EXPECT_NEAR(v.xp, PolynomialEop::xp(mjd), 1e-21);
        EXPECT_NEAR(v.dut1, PolynomialEop::dut1(mjd), 5e-14);
        EXPECT_NEAR(v.lod, PolynomialEop::lod(mjd), 5e-18);
        EXPECT_NEAR(v.dX, PolynomialEop::dX(mjd), 2e-24);
        EXPECT_NEAR(v.dY, PolynomialEop::dY(mjd), 2e-24);
    }
}

TEST(EopTable, Lagrange4ErrorOnAQuarticIsTheNodePolynomial) {
    // For f(u) = ... + c4 u^4 the interpolation error of a cubic through four nodes x_k is exactly
    // f''''/4! prod(x - x_k) = c4 prod(x - x_k). The nodes are the two records on each side of the
    // epoch, or the first (last) four records in the first (last) interval.
    const EopTable table(PolynomialEop::records(), EopInterpolation::Lagrange4);
    for (double mjd : {60000.3, 60001.5, 60003.7, 60004.2, 60006.8}) {
        SCOPED_TRACE(testing::Message() << "MJD " << mjd);
        const double lower          = std::clamp(std::floor(mjd) - 1.0, 60000.0, 60004.0); // first of the four nodes
        double       nodePolynomial = 1.0;
        for (int k = 0; k < 4; ++k) {
            nodePolynomial *= mjd - (lower + k);
        }
        const double expected = PolynomialEop::yp(mjd) - arcsec2rad(PolynomialEop::c4 * nodePolynomial);
        // The error itself is 1.7e-9 to 4.7e-9 rad (0.35 to 0.96 mas) at these epochs; rounding is
        // 4.2e-22 rad (measured).
        EXPECT_NEAR(table.at(mjd).yp, expected, 4e-21);
        EXPECT_GT(std::abs(PolynomialEop::yp(mjd) - expected), 1e-9); // the error term is not negligible
    }
}

TEST(EopTable, LinearInterpolationIsTheChord) {
    const EopTable table(PolynomialEop::records(), EopInterpolation::Linear);
    for (double mjd : {60000.3, 60003.7, 60006.8}) {
        const double x0 = std::floor(mjd), t = mjd - x0;
        EXPECT_NEAR(table.at(mjd).yp, (1.0 - t) * PolynomialEop::yp(x0) + t * PolynomialEop::yp(x0 + 1.0), 4e-21)
            << mjd;
        EXPECT_NEAR(table.at(mjd).dut1, (1.0 - t) * PolynomialEop::dut1(x0) + t * PolynomialEop::dut1(x0 + 1.0), 5e-14)
            << mjd;
    }
}

TEST(FrameTransformer, InvalidFrameValueGivesNaN) {
    // A Frame that is not one of the enumerators (a bad cast, an uninitialised value) used to index
    // past the end of the transformer's table of matrices. gcrsTo() and matrix() are noexcept, so they
    // answer with NaN, which cannot be taken for a rotation; the functions that may throw, throw.
    const Time             t = Time::fromIso8601("2021-07-04T12:34:56Z");
    const FrameTransformer ft(t, ConstantEop{});
    const auto             allNaN = [](const Mat3& m) {
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 3; ++j) {
                if (!std::isnan(m(i, j))) {
                    return false;
                }
            }
        }
        return true;
    };
    const auto        allNaNVec = [](const Vec3& v) { return std::isnan(v.x) && std::isnan(v.y) && std::isnan(v.z); };
    const StateVector s{
        {7000e3, -1200e3, 300e3},
        {1.0e3,  7.2e3,   0.5e3}
    };
    for (int raw : {9, 99, -1}) { // 9 is one past the last enumerator (ITRS = 8)
        const Frame bad = static_cast<Frame>(raw);
        EXPECT_TRUE(allNaN(ft.gcrsTo(bad))) << raw;
        EXPECT_TRUE(allNaN(ft.matrix(bad, Frame::GCRS))) << raw;
        EXPECT_TRUE(allNaN(ft.matrix(Frame::ITRS, bad))) << raw;
        EXPECT_TRUE(allNaNVec(ft.transform(s.r, bad, Frame::ITRS))) << raw;
        EXPECT_TRUE(allNaNVec(ft.transform(s.r, Frame::GCRS, bad))) << raw;
        for (const StateVector& out : {ft.transform(s, bad, Frame::ITRS), ft.transform(s, Frame::GCRS, bad)}) {
            EXPECT_TRUE(allNaNVec(out.r)) << raw;
            EXPECT_TRUE(allNaNVec(out.v)) << raw;
        }
        EXPECT_THROW((void)gcrsToCelestialMatrix(bad, t), Error) << raw;
        EXPECT_THROW((void)poleOfDateMatrix(bad, t), Error) << raw;
    }
}
