/**********************************************************************
 * enum_io_test.cpp
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

#include <iomanip>
#include <sstream>
#include <string>

#include "astro/astrometry/rise_set.h"
#include "astro/core/rotation.h"
#include "astro/ephemeris/ephemeris.h"
#include "astro/frames/eop.h"
#include "astro/frames/precession_nutation.h"
#include "astro/frames/transform.h"
#include "astro/geodesy/helmert.h"
#include "astro/orbit/secular.h"
#include "astro/orbit/sgp4.h"
#include "astro/orbit/state12.h"
#include "astro/satellite/eclipse.h"
#include "astro/satellite/estimate_tracker.h"
#include "astro/time/time.h"

#include <boost/describe/enumerators.hpp>
#include <boost/mp11/algorithm.hpp>

using namespace astro;

namespace {

template <class E>
std::string str(E e) {
    std::ostringstream os;
    os << e;
    return os.str();
}

template <class E>
E invalid(int v) {
    return static_cast<E>(v);
}

} // namespace

TEST(EnumIo, WritesQualifiedEnumerator) {
    EXPECT_EQ(str(PrecessionNutationModel::IAU2006_2000A), "PrecessionNutationModel::IAU2006_2000A");
    EXPECT_EQ(str(PrecessionNutationModel::IAU2006_2000B), "PrecessionNutationModel::IAU2006_2000B");
    EXPECT_EQ(str(Frame::GCRS), "Frame::GCRS");
    EXPECT_EQ(str(Frame::ITRS), "Frame::ITRS");
    EXPECT_EQ(str(EulerSequence::XYZ), "EulerSequence::XYZ");
    EXPECT_EQ(str(EulerSequence::ZYZ), "EulerSequence::ZYZ");
    EXPECT_EQ(str(EopInterpolation::Linear), "EopInterpolation::Linear");
    EXPECT_EQ(str(EopInterpolation::Lagrange4), "EopInterpolation::Lagrange4");
    EXPECT_EQ(str(HelmertConvention::PositionVector), "HelmertConvention::PositionVector");
    EXPECT_EQ(str(HelmertConvention::CoordinateFrame), "HelmertConvention::CoordinateFrame");
    EXPECT_EQ(str(Body::SolarSystemBarycenter), "Body::SolarSystemBarycenter");
    EXPECT_EQ(str(Body::Moon), "Body::Moon");
    EXPECT_EQ(str(Body::Mars), "Body::Mars");
    EXPECT_EQ(str(SkyEventType::Rise), "SkyEventType::Rise");
    EXPECT_EQ(str(SkyEventType::Transit), "SkyEventType::Transit");
    EXPECT_EQ(str(Twilight::Civil), "Twilight::Civil");
    EXPECT_EQ(str(Twilight::Astronomical), "Twilight::Astronomical");
    EXPECT_EQ(str(ShadowModel::Cylindrical), "ShadowModel::Cylindrical");
    EXPECT_EQ(str(ShadowModel::Conical), "ShadowModel::Conical");
    EXPECT_EQ(str(ShadowState::Sunlit), "ShadowState::Sunlit");
    EXPECT_EQ(str(ShadowState::Annular), "ShadowState::Annular");
    EXPECT_EQ(str(TimeScale::UTC), "TimeScale::UTC");
    EXPECT_EQ(str(TimeScale::BDT), "TimeScale::BDT");
    EXPECT_EQ(str(GnssSystem::GPS), "GnssSystem::GPS");
    EXPECT_EQ(str(GnssSystem::BeiDou), "GnssSystem::BeiDou");
    EXPECT_EQ(str(SecularOrder::J2), "SecularOrder::J2");
    EXPECT_EQ(str(SecularOrder::J2SquaredJ4), "SecularOrder::J2SquaredJ4");
    EXPECT_EQ(str(Sgp4Gravity::WGS72Old), "Sgp4Gravity::WGS72Old");
    EXPECT_EQ(str(Sgp4Gravity::WGS84), "Sgp4Gravity::WGS84");
    EXPECT_EQ(str(Sgp4OpsMode::Improved), "Sgp4OpsMode::Improved");
    EXPECT_EQ(str(Sgp4OpsMode::Afspc), "Sgp4OpsMode::Afspc");
    EXPECT_EQ(str(HigherOrderMode::ModelOnly), "HigherOrderMode::ModelOnly");
    EXPECT_EQ(str(HigherOrderMode::UnmodeledForce), "HigherOrderMode::UnmodeledForce");
    EXPECT_EQ(str(ConsistencyTest::Position), "ConsistencyTest::Position");
    EXPECT_EQ(str(ConsistencyTest::PositionVelocity), "ConsistencyTest::PositionVelocity");
    EXPECT_EQ(str(ConsistencyTest::Elements), "ConsistencyTest::Elements");
}

TEST(EnumIo, WritesInvalidValue) {
    EXPECT_EQ(str(invalid<PrecessionNutationModel>(2)), "PrecessionNutationModel::INVALID_VALUE");
    EXPECT_EQ(str(invalid<Frame>(-1)), "Frame::INVALID_VALUE");
    EXPECT_EQ(str(invalid<EulerSequence>(12)), "EulerSequence::INVALID_VALUE");
    EXPECT_EQ(str(invalid<EopInterpolation>(7)), "EopInterpolation::INVALID_VALUE");
    EXPECT_EQ(str(invalid<HelmertConvention>(2)), "HelmertConvention::INVALID_VALUE");
    EXPECT_EQ(str(invalid<Body>(11)), "Body::INVALID_VALUE"); // gap in the NAIF ids
    EXPECT_EQ(str(invalid<SkyEventType>(3)), "SkyEventType::INVALID_VALUE");
    EXPECT_EQ(str(invalid<Twilight>(3)), "Twilight::INVALID_VALUE");
    EXPECT_EQ(str(invalid<ShadowModel>(2)), "ShadowModel::INVALID_VALUE");
    EXPECT_EQ(str(invalid<ShadowState>(4)), "ShadowState::INVALID_VALUE");
    EXPECT_EQ(str(invalid<TimeScale>(10)), "TimeScale::INVALID_VALUE");
    EXPECT_EQ(str(invalid<GnssSystem>(3)), "GnssSystem::INVALID_VALUE");
    EXPECT_EQ(str(invalid<SecularOrder>(2)), "SecularOrder::INVALID_VALUE");
    EXPECT_EQ(str(invalid<Sgp4Gravity>(3)), "Sgp4Gravity::INVALID_VALUE");
    EXPECT_EQ(str(invalid<Sgp4OpsMode>(2)), "Sgp4OpsMode::INVALID_VALUE");
    EXPECT_EQ(str(invalid<HigherOrderMode>(2)), "HigherOrderMode::INVALID_VALUE");
    EXPECT_EQ(str(invalid<ConsistencyTest>(3)), "ConsistencyTest::INVALID_VALUE");
}

TEST(EnumIo, FieldWidthAppliesToWholeToken) {
    std::ostringstream os;
    os << std::setw(16) << std::left << TimeScale::TT << '|' << std::setw(16) << std::right << Frame::TEME;
    EXPECT_EQ(os.str(), "TimeScale::TT   |     Frame::TEME");
}

TEST(EnumIo, ToStringStillUnqualified) {
    EXPECT_EQ(toString(TimeScale::TAI), "TAI");
    EXPECT_EQ(toString(Frame::CIRS), "CIRS");
    EXPECT_EQ(toString(Body::Earth), "Earth");
    EXPECT_EQ(toString(invalid<TimeScale>(42)), "?");
}

namespace {

/// operator<< must agree with Boost.Describe for every described enumerator.
template <class E>
void checkAgainstDescribe(const char* enumName, std::size_t expectedCount) {
    using Enumerators = boost::describe::describe_enumerators<E>;
    EXPECT_EQ(boost::mp11::mp_size<Enumerators>::value, expectedCount) << enumName;
    boost::mp11::mp_for_each<Enumerators>(
        [&](auto d) { EXPECT_EQ(str(d.value), std::string(enumName) + "::" + d.name); });
}

} // namespace

TEST(EnumIo, AgreesWithBoostDescribe) {
    static_assert(boost::describe::has_describe_enumerators<TimeScale>::value);
    checkAgainstDescribe<PrecessionNutationModel>("PrecessionNutationModel", 2);
    checkAgainstDescribe<Frame>("Frame", 9);
    checkAgainstDescribe<EulerSequence>("EulerSequence", 12);
    checkAgainstDescribe<EopInterpolation>("EopInterpolation", 2);
    checkAgainstDescribe<HelmertConvention>("HelmertConvention", 2);
    checkAgainstDescribe<Body>("Body", 16);
    checkAgainstDescribe<SkyEventType>("SkyEventType", 3);
    checkAgainstDescribe<Twilight>("Twilight", 3);
    checkAgainstDescribe<ShadowModel>("ShadowModel", 2);
    checkAgainstDescribe<ShadowState>("ShadowState", 4);
    checkAgainstDescribe<TimeScale>("TimeScale", 10);
    checkAgainstDescribe<GnssSystem>("GnssSystem", 3);
    checkAgainstDescribe<SecularOrder>("SecularOrder", 2);
    checkAgainstDescribe<Sgp4Gravity>("Sgp4Gravity", 3);
    checkAgainstDescribe<Sgp4OpsMode>("Sgp4OpsMode", 2);
    checkAgainstDescribe<HigherOrderMode>("HigherOrderMode", 2);
    checkAgainstDescribe<ConsistencyTest>("ConsistencyTest", 3);
}
