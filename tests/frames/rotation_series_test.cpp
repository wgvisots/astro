/**********************************************************************
 * rotation_series_test.cpp
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

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

#include "astro/core/constants.h"
#include "astro/core/error.h"
#include "astro/frames/rotation_series.h"
#include "astro/frames/transform.h"
#include "test_util.h"

using namespace astro;

namespace {

const Time kStart = Time::fromIso8601("2026-03-01T00:00:00Z");

/// Self-consistent EOP: UT1 - UTC drifts at the rate the length of day implies.
class DriftingEop final : public EopProvider {
public:
    [[nodiscard]] EopValues at(double mjdUtc) const override {
        EopValues v;
        v.xp   = 0.12 * constants::arcsecToRad;
        v.yp   = 0.38 * constants::arcsecToRad;
        v.lod  = kLod;
        v.dut1 = 0.21 - kLod * (mjdUtc - 61100.0); // seconds per day of excess length
        v.dX   = 0.2 * constants::masToRad;
        v.dY   = -0.1 * constants::masToRad;
        return v;
    }
    static constexpr double kLod = 1.3e-3;
};

double maxAbsDiff(const Mat3& a, const Mat3& b) {
    double m = 0.0;
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            m = std::max(m, std::abs(a(i, j) - b(i, j)));
        }
    }
    return m;
}

} // namespace

TEST(ItrsRotationSeries, ExactAtTheNodes) {
    // At the nodes (every nodeSpacing = 1 h from the start) the parts are those of a FrameTransformer at
    // that time, so the matrices agree to rounding: measured 4.3e-15 (the Earth rotation angle is
    // carried from node to node), tolerance 2e-14.
    const DriftingEop        eop;
    const ItrsRotationSeries series(kStart, 6 * 3600.0, eop);
    double                   worst = 0.0;
    for (int i = 0; i <= 6; ++i) {
        const double dt = 3600.0 * i;
        worst =
            std::max(worst,
                     maxAbsDiff(series.gcrsToItrsMatrix(dt), FrameTransformer(kStart + dt, eop).gcrsTo(Frame::ITRS)));
    }
    EXPECT_LT(worst, 2e-14);
}

TEST(ItrsRotationSeries, MatchesFrameTransformerBetweenTheNodes) {
    // Between the nodes precession-nutation and polar motion are linear in time: within a few 1e-11
    // rad of the transformer (the header's 0.2 mm at the surface), at random times.
    const DriftingEop                      eop;
    const ItrsRotationSeries               series(kStart, 6 * 3600.0, eop);
    std::mt19937_64                        rng(3);
    std::uniform_real_distribution<double> u(0.0, 6 * 3600.0);
    const StateVector                      s{
                             {6.0e6, -2.0e6, 1.5e6 },
                             {1.2e3, 5.1e3,  -4.7e3}
    };
    double worstM = 0.0, worstV = 0.0;
    for (int i = 0; i < 200; ++i) {
        const double           dt = u(rng);
        const FrameTransformer ft(kStart + dt, eop);
        worstM                  = std::max(worstM, maxAbsDiff(series.gcrsToItrsMatrix(dt), ft.gcrsTo(Frame::ITRS)));
        const StateVector exact = ft.transform(s, Frame::GCRS, Frame::ITRS);
        worstV                  = std::max(worstV, norm(series.gcrsToItrs(s, dt).v - exact.v));
    }
    // Measured 2.8e-11 (0.18 mm at the surface) and 1.4e-7 m/s, the matrix error times the 7 km/s
    // speed and the Earth-rotation term.
    EXPECT_LT(worstM, 5e-11);
    EXPECT_LT(worstV, 3e-7);
}

TEST(ItrsRotationSeries, RoundTripRateAndErrors) {
    const DriftingEop        eop;
    const ItrsRotationSeries series(kStart, 1800.0, eop, {PrecessionNutationModel::IAU2006_2000B, true});
    const StateVector        s{
               {-4.1e6, 3.3e6, 4.4e6 },
               {2.0e3,  6.0e3, -1.0e3}
    };
    const StateVector back = series.itrsToGcrs(series.gcrsToItrs(s, 900.0), 900.0);
    test::expectVecNear(back.r, s.r, 1e-8);
    test::expectVecNear(back.v, s.v, 1e-11);
    test::expectVecNear(series.gcrsToItrs(s.r, 900.0), series.gcrsToItrs(s, 900.0).r, 1e-8);
    // The rotation rate includes the length of day: omega = 7.292115e-5 (1 - LOD / 86400).
    EXPECT_NEAR(series.rotationRate(100.0) / (7.2921151467064e-5 * (1.0 - DriftingEop::kLod / 86400.0)), 1.0, 1e-11);
    EXPECT_DOUBLE_EQ(series.duration(), 1800.0);
    EXPECT_EQ(series.start(), kStart);

    EXPECT_THROW((void)series.gcrsToItrsMatrix(1802.0), Error);
    EXPECT_THROW((void)series.gcrsToItrsMatrix(-2.0), Error);
    EXPECT_THROW(ItrsRotationSeries(kStart, -1.0, eop), Error);
    EXPECT_THROW(ItrsRotationSeries(kStart, 100.0, eop, {}, 0.0), Error);
    // A zero-length window still works (one second wide).
    const ItrsRotationSeries point(kStart, 0.0, eop);
    EXPECT_LT(maxAbsDiff(point.gcrsToItrsMatrix(0.0), FrameTransformer(kStart, eop).gcrsTo(Frame::ITRS)), 1e-15);
}

TEST(ItrsRotationSeries, EopNotCoveringTheWindowOrAStartBefore1960Throws) {
    // Documented: the EOP must cover the window, and the window must start in 1960 or later.
    std::vector<EopRecord> records;
    for (int d = 0; d < 3; ++d) {
        records.push_back({kStart.mjd() + d, {}}); // three days from the start
    }
    const EopTable table(records);
    EXPECT_NO_THROW(ItrsRotationSeries(kStart, 86400.0, table));
    EXPECT_THROW(ItrsRotationSeries(kStart, 4 * 86400.0, table), Error);
    EXPECT_THROW(ItrsRotationSeries(Time::fromCalendar(TimeScale::TT, 1955, 1, 1, 0), 3600.0, ConstantEop{}), Error);
}
