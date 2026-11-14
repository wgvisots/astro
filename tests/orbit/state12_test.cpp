/**********************************************************************
 * state12_test.cpp
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

#include <memory>
#include <sstream>
#include <type_traits>
#include <vector>

#include "astro/core/angle.h"
#include "astro/orbit/elements.h"
#include "astro/orbit/gravity.h"
#include "astro/orbit/numerical.h"
#include "astro/orbit/secular.h"
#include "astro/orbit/state12.h"
#include "astro/orbit/two_body.h"
#include "test_util.h"

using namespace astro;
using astro::test::expectVecNear;

namespace {

constexpr double  km    = 1000.0;
const CentralBody kJ2J4 = CentralBody::earthJ2J4();

// An eccentric, inclined LEO orbit.
StateVector leo() {
    return elementsToState({7200 * km, 0.05, deg2rad(51.6), 0.3, 1.2, 0.4});
}

StateVector12 pointMassState(const StateVector& s) {
    return {s.r, s.v, pointMassAcceleration(s.r), pointMassJerk(s.r, s.v)};
}

StateVector12 fieldState(const ForceModel& f, const StateVector& s) {
    const Vec3 a = f.acceleration(0.0, s.r, s.v);
    return {s.r, s.v, a, f.jerk(0.0, s.r, s.v, a)};
}

// Taylor series through jerk; the error is O(dt^4) times the snap.
StateVector taylor(const StateVector12& s, double dt) {
    return {s.r + s.v * dt + s.a * (0.5 * dt * dt) + s.j * (dt * dt * dt / 6.0),
            s.v + s.a * dt + s.j * (0.5 * dt * dt)};
}

// A few mm/s^2 of slowly changing thrust on top of the model.
StateVector12 withThrust(StateVector12 s) {
    s.a += Vec3{2e-3, -1e-3, 3e-3};
    s.j += Vec3{1e-6, 2e-6, -1e-6};
    return s;
}

} // namespace

// {r, v} must still mean a 6-component state, so existing calls stay unambiguous.
static_assert(std::is_same_v<decltype(propagateKepler(
                                 {
                                     Vec3{7e6, 0,     0},
                                     Vec3{0,   7.5e3, 0}
},
                                 60.0)),
                             StateVector>);

TEST(State12, PointMassJerkIsDerivativeOfAcceleration) {
    const StateVector s = leo();
    const double      h = 0.5;
    const Vec3        fd =
        (pointMassAcceleration(propagateKepler(s, h).r) - pointMassAcceleration(propagateKepler(s, -h).r)) / (2 * h);
    expectVecNear(pointMassJerk(s.r, s.v), fd, 1e-9); // jerk ~1e-2 m/s^3
}

TEST(State12, ForceModelJerkMatchesTrajectoryDerivative) {
    // J2/J4 field jerk against a central difference of acceleration along the integrated orbit.
    const auto                field = std::make_shared<GravityField>(kJ2J4);
    const NumericalPropagator prop(Time(), leo(), {field});
    const double              h = 1.0;
    const auto                s = prop.propagate12(std::vector<double>{-h, 0.0, h});
    expectVecNear(s[1].j, (s[2].a - s[0].a) / (2 * h), 1e-8);
    // The generic central-difference default agrees with GravityField's partly exact override.
    expectVecNear(field->ForceModel::jerk(0.0, s[1].r, s[1].v, s[1].a), s[1].j, 1e-11);
}

TEST(State12, KeplerModelOnly) {
    const StateVector s0       = leo();
    const double      dt       = 2000.0;
    const StateVector expected = propagateKepler(s0, dt);

    // Consistent input: no residual; output a and j are the point mass at the new state.
    const Propagation12 p = propagateKepler(pointMassState(s0), dt);
    EXPECT_EQ(p.state, pointMassState(expected));
    EXPECT_EQ(p.accelerationResidual, Vec3{});
    EXPECT_EQ(p.jerkResidual, Vec3{});

    // Inconsistent input is reported but does not change the trajectory.
    const Propagation12 q = propagateKepler(withThrust(pointMassState(s0)), dt);
    EXPECT_EQ(q.state, p.state);
    expectVecNear(q.accelerationResidual, {2e-3, -1e-3, 3e-3}, 1e-15);
    expectVecNear(q.jerkResidual, {1e-6, 2e-6, -1e-6}, 1e-15);
}

TEST(State12, KeplerUnmodeledForce) {
    const StateVector   s0   = leo();
    const StateVector12 in   = withThrust(pointMassState(s0));
    constexpr auto      mode = HigherOrderMode::UnmodeledForce;

    // At dt = 0 the input is reproduced, acceleration and jerk included.
    const Propagation12 p0 = propagateKepler(in, 0.0, constants::gmEarth, mode);
    EXPECT_EQ(p0.state.rv(), in.rv());
    expectVecNear(p0.state.a, in.a, 1e-15);
    expectVecNear(p0.state.j, in.j, 1e-15);

    // Over a short arc the trajectory follows the input's own Taylor series (error ~ dt^4 snap).
    const double        dt = 10.0;
    const StateVector   t  = taylor(in, dt);
    const Propagation12 p  = propagateKepler(in, dt, constants::gmEarth, mode);
    expectVecNear(p.state.r, t.r, 5e-3);
    expectVecNear(p.state.v, t.v, 2e-3);
    // ModelOnly ignores the thrust and misses by about da dt^2 / 2 (~0.2 m).
    EXPECT_GT(norm(propagateKepler(in, dt).state.r - t.r), 0.1);

    // With a consistent input there is no unmodeled force and the result is the analytic one.
    const StateVector12 consistent = pointMassState(s0);
    EXPECT_EQ(propagateKepler(consistent, 3000.0, constants::gmEarth, mode).state,
              propagateKepler(consistent, 3000.0).state);
}

TEST(State12, KeplerUnmodeledForceDoesWorkOnOrbit) {
    // Constant along-track thrust over half an orbit: the change in specific orbital energy equals
    // the work done, the integral of da . v dt.
    const StateVector s0          = leo();
    StateVector12     in          = pointMassState(s0);
    const Vec3        da          = normalized(s0.v) * 1e-3;
    in.a                         += da;
    const double              T   = 3000.0;
    const int                 n   = 300; // Simpson's rule
    const std::vector<double> ts  = [&] {
        std::vector<double> v;
        for (int k = 0; k <= n; ++k) {
            v.push_back(T * k / n);
        }
        return v;
    }();
    const NumericalPropagator prop(Time(),
                                   in,
                                   {std::make_shared<GravityField>(CentralBody::pointMass())},
                                   HigherOrderMode::UnmodeledForce);
    const auto                states = prop.propagate12(ts);
    double                    work   = 0.0;
    for (int k = 0; k <= n; ++k) {
        expectVecNear(states[k].a - pointMassAcceleration(states[k].r), da, 1e-13); // ~1 ulp of 8 m/s^2
        work += ((k == 0 || k == n) ? 1.0 : (k % 2 ? 4.0 : 2.0)) * dot(da, states[k].v);
    }
    work        *= T / n / 3.0;
    auto energy  = [](const StateVector12& s) { return 0.5 * norm2(s.v) - constants::gmEarth / norm(s.r); };
    EXPECT_NEAR((energy(states.back()) - energy(in)) / work, 1.0, 1e-8);
    // The free function gives the same end state.
    const StateVector12 end = propagateKepler(in, T, constants::gmEarth, HigherOrderMode::UnmodeledForce).state;
    expectVecNear(end.r, states.back().r, 1e-4); // integrator tolerance, different steps
}

TEST(State12, SecularStateVersionMatchesElementVersion) {
    const KeplerianElements el{7000 * km, 0.01, deg2rad(98.0), 0.5, 1.0, 2.0};
    const double            dt       = 5 * 86400.0;
    const StateVector       expected = elementsToState(propagateSecular(el, dt));
    const StateVector       s        = propagateSecular(elementsToState(el), dt);
    expectVecNear(s.r, expected.r, 1e-6);
    expectVecNear(s.v, expected.v, 1e-9);

    // SecularOrder::J2 uses J2 only, in both the drift and the output acceleration.
    const CentralBody body{
        .J = {0.0, 0.0, constants::J2}
    };
    const Propagation12 p =
        propagateSecular(fieldState(GravityField(body), elementsToState(el)), dt, body, SecularOrder::J2);
    expectVecNear(p.state.r, elementsToState(propagateSecular(el, dt, body, SecularOrder::J2)).r, 1e-6);
    expectVecNear(p.accelerationResidual, {}, 1e-18);
    expectVecNear(p.state.a, GravityField(body).acceleration(0.0, p.state.r, p.state.v), 1e-18);
}

TEST(State12, SecularModelOnlyAndUnmodeledForce) {
    const GravityField field(kJ2J4);
    const StateVector  s0 = leo();
    const double       dt = 600.0;

    // Consistent input: zero residuals; output a and j are the J2+J4 field at the new state.
    const Propagation12 p = propagateSecular(fieldState(field, s0), dt);
    expectVecNear(p.accelerationResidual, {}, 1e-18);
    expectVecNear(p.jerkResidual, {}, 1e-18);
    const StateVector sec = propagateSecular(s0, dt);
    EXPECT_EQ(p.state.rv(), sec);
    const StateVector12 expected = fieldState(field, sec);
    expectVecNear(p.state.a, expected.a, 1e-18);
    expectVecNear(p.state.j, expected.j, 1e-18);

    // Thrust: residuals reported; UnmodeledForce reproduces the input at dt = 0 and superposes
    // da dt^2/2 + dj dt^3/6 on the secular trajectory.
    const StateVector12 in = withThrust(fieldState(field, s0));
    const Propagation12 q0 = propagateSecular(in, 0.0, {}, SecularOrder::J2SquaredJ4, HigherOrderMode::UnmodeledForce);
    expectVecNear(q0.state.r, in.r, 1e-6);
    expectVecNear(q0.state.v, in.v, 1e-9);
    expectVecNear(q0.state.a, in.a, 1e-13); // r went through elements and back
    expectVecNear(q0.state.j, in.j, 1e-13);
    const Propagation12 q = propagateSecular(in, dt, {}, SecularOrder::J2SquaredJ4, HigherOrderMode::UnmodeledForce);
    expectVecNear(q.accelerationResidual, {2e-3, -1e-3, 3e-3}, 1e-15);
    expectVecNear(q.state.r - sec.r,
                  Vec3{2e-3, -1e-3, 3e-3} * (dt * dt / 2) + Vec3{1e-6, 2e-6, -1e-6} * (dt * dt * dt / 6),
                  1e-6);
    // ModelOnly ignores the thrust.
    EXPECT_EQ(propagateSecular(in, dt).state.rv(), sec);
}

TEST(State12, NumericalPropagatorModes) {
    const StateVector   s0         = leo();
    const auto          field      = std::make_shared<GravityField>(kJ2J4);
    const StateVector12 consistent = fieldState(*field, s0);
    const double        dt         = 5400.0;

    // ModelOnly from a consistent 12-state matches the 6-state propagator.
    const NumericalPropagator six(Time(), s0, {field});
    const NumericalPropagator twelve = NumericalPropagator::j2j4(Time(), consistent, Frame::TOD);
    expectVecNear(twelve.accelerationResidual(), {}, 1e-18);
    expectVecNear(twelve.jerkResidual(), {}, 1e-18);
    const StateVector12 a = twelve.propagate12(dt);
    EXPECT_EQ(a.rv(), six.propagate(dt));
    EXPECT_EQ(a, six.withDerivatives(dt, six.propagate(dt)));
    EXPECT_EQ(six.accelerationResidual(), Vec3{});

    // With thrust, ModelOnly reports it and ignores it; UnmodeledForce follows it.
    const StateVector12       in = withThrust(consistent);
    const NumericalPropagator ignore(Time(), in, {field});
    expectVecNear(ignore.accelerationResidual(), {2e-3, -1e-3, 3e-3}, 1e-15);
    EXPECT_EQ(ignore.propagate12(dt), a);
    const NumericalPropagator thrust =
        NumericalPropagator::j2j4(Time(), in, Frame::TOD, HigherOrderMode::UnmodeledForce);
    const StateVector12 b0 = thrust.propagate12(0.0);
    EXPECT_EQ(b0.rv(), in.rv());
    expectVecNear(b0.a, in.a, 1e-15);
    expectVecNear(b0.j, in.j, 1e-15);
    const double      h = 10.0;
    const StateVector t = taylor(in, h);
    expectVecNear(thrust.propagate(h).r, t.r, 5e-3);
    EXPECT_GT(norm(thrust.propagate(dt).r - a.r), 1e3); // thrust moves it km over 1.5 h
}

TEST(State12, ModeWritesQualifiedName) {
    std::ostringstream os;
    os << HigherOrderMode::UnmodeledForce;
    EXPECT_EQ(os.str(), "HigherOrderMode::UnmodeledForce");
}
