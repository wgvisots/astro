/**********************************************************************
 * event_test.cpp
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

#include <cmath>
#include <memory>
#include <sstream>

#include "astro/core/angle.h"
#include "astro/core/error.h"
#include "astro/orbit/elements.h"
#include "astro/orbit/gravity.h"
#include "astro/orbit/numerical.h"
#include "test_util.h"

using namespace astro;

namespace {

constexpr double km = 1000.0;

/// Point-mass propagator for a circular inclined orbit starting at its ascending node, so the
/// nodes (z = 0) are at exact multiples of half a period.
struct Circular {
    double              a{7000 * km};
    double              period{constants::twoPi * std::sqrt(a * a * a / constants::gmEarth)};
    NumericalPropagator prop{Time(),
                             elementsToState({a, 0.0, deg2rad(51.6), 0.3, 0.0, 0.0}),
                             {std::make_shared<GravityField>(std::vector<double>{}, constants::gmEarth)}};
};

double z(double, const StateVector& s) {
    return s.r.z;
}

} // namespace

TEST(PropagationEvent, FindsNodesInEachDirection) {
    const Circular c;
    // The epoch is on the ascending node, which is not an event itself.
    const auto rising = c.prop.propagateToEvent(z, 1.5 * c.period, Crossing::Rising);
    ASSERT_TRUE(rising.has_value());
    EXPECT_NEAR(rising->dt, c.period, 1e-5);
    const auto falling = c.prop.propagateToEvent(z, 1.5 * c.period, Crossing::Falling);
    ASSERT_TRUE(falling.has_value());
    EXPECT_NEAR(falling->dt, 0.5 * c.period, 1e-5);
    const auto any = c.prop.propagateToEvent(z, 1.5 * c.period);
    ASSERT_TRUE(any.has_value());
    EXPECT_NEAR(any->dt, 0.5 * c.period, 1e-5);
    // The event state is the trajectory at the event time.
    const StateVector direct = c.prop.propagate(falling->dt);
    test::expectVecNear(falling->state.r, direct.r, 1e-5);
    test::expectVecNear(falling->state.v, direct.v, 1e-8);
    EXPECT_NEAR(falling->state.r.z, 0.0, 1e-2); // |vz| ~ 6 km/s times the 1e-6 s tolerance

    // Backward: the descending node half a period ago, the ascending node a period ago. Directions
    // refer to increasing time.
    const auto backRising = c.prop.propagateToEvent(z, -1.5 * c.period, Crossing::Rising);
    ASSERT_TRUE(backRising.has_value());
    EXPECT_NEAR(backRising->dt, -c.period, 1e-5);
    const auto backFalling = c.prop.propagateToEvent(z, -1.5 * c.period, Crossing::Falling);
    ASSERT_TRUE(backFalling.has_value());
    EXPECT_NEAR(backFalling->dt, -0.5 * c.period, 1e-5);
}

TEST(PropagationEvent, NoEventAndErrors) {
    const Circular c;
    AdaptiveStats  st{9, 9, 9};
    EXPECT_FALSE(c.prop.propagateToEvent(z, 0.4 * c.period, Crossing::Any, 1e-6, &st).has_value());
    EXPECT_GT(st.accepted, 0);
    EXPECT_FALSE(c.prop.propagateToEvent(z, 0.0).has_value());
    EXPECT_THROW((void)c.prop.propagateToEvent(EventFunction{}, 100.0), Error);
    EXPECT_THROW((void)c.prop.propagateToEvent(z, std::nan("")), Error);
    EXPECT_THROW((void)c.prop.propagateToEvent(z, 100.0, Crossing::Any, 0.0), Error);
    EXPECT_THROW((void)c.prop.propagateToEvent([](double, const StateVector&) { return std::nan(""); }, 100.0), Error);
    std::ostringstream os;
    os << Crossing::Falling << ' ' << static_cast<Crossing>(42);
    EXPECT_EQ(os.str(), "Crossing::Falling Crossing::INVALID_VALUE");
}
