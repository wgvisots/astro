/**********************************************************************
 * integrator_test.cpp
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

// The adaptive integrator at the edges of its input: step bounds, values that are not numbers, and
// the one thing its error estimate cannot see.

#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "astro/core/angle.h"
#include "astro/core/constants.h"
#include "astro/core/error.h"
#include "astro/orbit/elements.h"
#include "astro/orbit/integrators.h"
#include "astro/orbit/numerical.h"
#include "astro/orbit/two_body.h"

using namespace astro;

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kInf = std::numeric_limits<double>::infinity();

/// y' = y, counting the evaluations.
struct Growth {
    int*        calls;
    OdeState<1> operator()(double /*t*/, const OdeState<1>& y) const {
        ++*calls;
        return {y[0]};
    }
};

} // namespace

TEST(Rkf78, MaxStepBoundsEveryStepIncludingTheFirst) {
    // y' = 1: any step passes the error test, so only maxStep limits it. The first step used to be
    // min(initialStep, |t1 - t0|) whatever maxStep said: the whole 10 s here.
    AdaptiveOptions opt;
    opt.maxStep        = 1.0;
    double     longest = 0.0;
    const auto y       = integrateRkf78Steps<1>(
        [](double, const OdeState<1>&) { return OdeState<1>{1.0}; },
        0.0,
        OdeState<1>{0.0},
        10.0,
        [&](double ta, const OdeState<1>&, double tb, const OdeState<1>&) { longest = std::max(longest, tb - ta); },
        opt);
    EXPECT_LE(longest, 1.0);
    EXPECT_NEAR(y[0], 10.0, 1e-13);
}

TEST(Rkf78, TimesThatAreNotFiniteAreRejectedAtOnce) {
    // These used to run the full ten million steps (seconds of work) before giving up.
    for (const auto& [t0, t1] : {
             std::pair{0.0,  kNaN },
             std::pair{0.0,  kInf },
             std::pair{0.0,  -kInf},
             std::pair{kNaN, 1.0  }
    }) {
        int calls = 0;
        EXPECT_THROW((void)integrateRkf78<1>(Growth{&calls}, t0, OdeState<1>{1.0}, t1), Error) << t0 << " " << t1;
        EXPECT_EQ(calls, 0) << t0 << " " << t1;
    }
}

TEST(Rkf78, UnusableOptionsAreRejectedAtOnce) {
    const auto rejected = [](const AdaptiveOptions& opt) {
        int calls = 0;
        EXPECT_THROW((void)integrateRkf78<1>(Growth{&calls}, 0.0, OdeState<1>{1.0}, 1.0, opt), Error);
        return calls == 0;
    };
    AdaptiveOptions opt;
    opt.initialStep = 0.0; // every step would have zero length
    EXPECT_TRUE(rejected(opt));
    opt             = {};
    opt.initialStep = kNaN;
    EXPECT_TRUE(rejected(opt));
    opt         = {};
    opt.maxStep = 0.0;
    EXPECT_TRUE(rejected(opt));
    opt         = {};
    opt.minStep = 10.0;
    opt.maxStep = 1.0;
    EXPECT_TRUE(rejected(opt));
    opt        = {};
    opt.relTol = 0.0;
    opt.absTol = 0.0; // the error would be divided by zero
    EXPECT_TRUE(rejected(opt));
    opt        = {};
    opt.relTol = -1e-12;
    EXPECT_TRUE(rejected(opt));
}

TEST(Rkf78, ADerivativeThatIsNotANumberIsAnError) {
    // The error norm was formed with std::max, which ignores a NaN operand: the step was accepted
    // and the integration returned NaN as if it had succeeded.
    const auto f = [](double t, const OdeState<2>& y) { return OdeState<2>{y[1], t > 5.0 ? kNaN : -y[0]}; };
    EXPECT_THROW((void)integrateRkf78<2>(f, 0.0, OdeState<2>{1.0, 0.0}, 10.0), Error);
    // Also in a component left out of the error control (a state transition matrix, say).
    AdaptiveOptions opt;
    opt.errorComponents = 1;
    const auto g        = [](double t, const OdeState<2>& y) { return OdeState<2>{-y[0], t > 5.0 ? kInf : 1.0}; };
    EXPECT_THROW((void)integrateRkf78<2>(g, 0.0, OdeState<2>{1.0, 0.0}, 10.0, opt), Error);
    // A start that is not a number.
    const auto h = [](double, const OdeState<1>& y) { return OdeState<1>{-y[0]}; };
    EXPECT_THROW((void)integrateRkf78<1>(h, 0.0, OdeState<1>{kNaN}, 10.0), Error);
}

TEST(Rkf78, ATrialStepThatLeavesTheDomainIsRetriedShorter) {
    // y' = -(sqrt y)^2 is y' = -y where y >= 0, and not a number below. The solution exp(-t) stays
    // positive, but the stages of the first trial step (60 s) overshoot far below zero. That is not
    // the end of the integration: a shorter step is tried. (It used to be accepted as it stood.)
    const auto f = [](double, const OdeState<1>& y) {
        const double root = std::sqrt(y[0]);
        return OdeState<1>{-root * root};
    };
    AdaptiveStats stats;
    const auto    y = integrateRkf78<1>(f, 0.0, OdeState<1>{1.0}, 10.0, {}, &stats);
    EXPECT_NEAR(y[0], std::exp(-10.0), 1e-8); // absTol is 1e-9 per step
    EXPECT_GT(stats.rejected, 0);
}

TEST(Rkf78, ATermInTimeAloneNeedsMaxStep) {
    // This pins a documented limitation, and shows the documented remedy works. If a later integrator
    // removes the limitation the first half of this test should go, with the paragraph in the header.
    // Fehlberg's error estimate cancels exactly for a part of f that depends
    // on t alone, so a short pulse can be stepped over with no step rejected. y' = pulse(t), a
    // Gaussian of width 0.5 s at t = 50 s with unit area: y(100) = 1.
    const auto pulse = [](double t, const OdeState<1>&) {
        const double u = (t - 50.0) / 0.5;
        return OdeState<1>{std::exp(-0.5 * u * u) / (0.5 * std::sqrt(constants::twoPi))};
    };
    AdaptiveStats stats;
    const auto    coarse = integrateRkf78<1>(pulse, 0.0, OdeState<1>{0.0}, 100.0, {}, &stats);
    EXPECT_EQ(stats.rejected, 0);               // the estimate never objects...
    EXPECT_GT(std::abs(coarse[0] - 1.0), 1e-3); // ...and the answer is wrong
    // The documented remedy: bound the step by a fraction of the time over which the term changes.
    AdaptiveOptions opt;
    opt.maxStep     = 0.5;
    const auto fine = integrateRkf78<1>(pulse, 0.0, OdeState<1>{0.0}, 100.0, opt);
    EXPECT_NEAR(fine[0], 1.0, 1e-10); // measured: 6e-12
}

TEST(NumericalPropagator, RejectsNullForceModelsAndEphemeris) {
    // Both used to be accepted and dereferenced at the first propagation.
    const StateVector s0 = elementsToState({7000e3, 0.001, deg2rad(51.6), 0.3, 0.6, 0.9});
    EXPECT_THROW(NumericalPropagator(Time(), s0, {std::shared_ptr<const ForceModel>{}}), Error);
    EXPECT_THROW(NumericalPropagator(Time(), s0, {std::make_shared<GravityField>(std::vector<double>{}), nullptr}),
                 Error);
    EXPECT_THROW(ThirdBodyForce(Body::Moon, 4.9e12, nullptr, Time()), Error);
}

TEST(NumericalPropagator, RejectsTimesThatAreNotFinite) {
    // These used to integrate for ten million steps before throwing. Count the force evaluations.
    struct CountingPointMass final : ForceModel {
        mutable int calls{0};
        Vec3        acceleration(double, const Vec3& r, const Vec3&) const override {
            ++calls;
            const double d = norm(r);
            return r * (-constants::gmEarth / (d * d * d));
        }
    };
    const auto                force = std::make_shared<CountingPointMass>();
    const NumericalPropagator prop(Time(), elementsToState({7000e3, 0.001, deg2rad(51.6), 0.3, 0.6, 0.9}), {force});
    EXPECT_THROW((void)prop.propagate(kNaN), Error);
    EXPECT_THROW((void)prop.propagate(kInf), Error);
    EXPECT_THROW((void)prop.propagateGrid(kNaN, 10.0, 5), Error);
    EXPECT_EQ(force->calls, 0);
}

TEST(Rk4, ErrorOnTheHarmonicOscillatorIsTheFifthOrderTermOfTheExponential) {
    // On y'' = -y one RK4 step multiplies (cos t, -sin t), seen as exp(-i t), by the Taylor
    // polynomial of exp(-i h) to fourth order, so the error of the step is the first term left out,
    // h^5 / 120 in magnitude (the next is h^6 / 720, 0.2 % at h = 0.01). It is also the error of
    // rk4Step itself, so this checks the coefficients of the method.
    const auto f = [](double, const OdeState<2>& y) { return OdeState<2>{y[1], -y[0]}; };
    for (const double h : {0.01, 0.02, -0.01}) {
        const double      t0    = 0.3;
        const OdeState<2> y     = rk4Step<2>(f, t0, {std::cos(t0), -std::sin(t0)}, h);
        const double      error = std::hypot(y[0] - std::cos(t0 + h), y[1] + std::sin(t0 + h));
        EXPECT_NEAR(error / (std::pow(std::abs(h), 5) / 120.0), 1.0, 0.005) << "h = " << h;
    }
}

TEST(Rk4, GlobalErrorIsFourthOrder) {
    // Over T the step errors add up along the orbit to T h^4 / 120 (no damping on the oscillator),
    // so halving the step divides the error by 16. Measured: 15.999, and the error within 1e-4 of
    // T h^4 / 120; a third-order method would halve it only to an eighth.
    const auto   f        = [](double, const OdeState<2>& y) { return OdeState<2>{y[1], -y[0]}; };
    const double T        = 10.0;
    double       previous = 0.0;
    for (const double h : {0.1, 0.05, 0.025}) {
        const OdeState<2> y     = integrateRk4<2>(f, 0.0, {1.0, 0.0}, T, h);
        const double      error = std::hypot(y[0] - std::cos(T), y[1] + std::sin(T));
        EXPECT_NEAR(error / (T * std::pow(h, 4) / 120.0), 1.0, 1e-3) << "h = " << h;
        if (previous > 0.0) {
            EXPECT_NEAR(previous / error, 16.0, 0.05) << "h = " << h;
        }
        previous = error;
    }
    // Backwards, with a step that does not divide the interval (the last step is shortened).
    const OdeState<2> back = integrateRk4<2>(f, T, {std::cos(T), -std::sin(T)}, 0.0, 0.03);
    EXPECT_LT(std::hypot(back[0] - 1.0, back[1]), 1.01 * T * std::pow(0.03, 4) / 120.0);
}

TEST(Rk4, RejectsAStepThatIsNotPositive) {
    const auto f = [](double, const OdeState<1>& y) { return OdeState<1>{y[0]}; };
    for (const double step : {0.0, -0.1}) {
        try {
            (void)integrateRk4<1>(f, 0.0, OdeState<1>{1.0}, 1.0, step);
            ADD_FAILURE() << "no exception for step " << step;
        } catch (const Error& e) {
            EXPECT_NE(std::string(e.what()).find("step must be positive"), std::string::npos) << e.what();
        }
    }
}

TEST(Rkf78, FixedStepErrorIsEighthOrder) {
    // The integrator propagates Fehlberg's eighth-order solution. With minStep = maxStep = h and
    // tolerances so loose that every step passes, it is a fixed-step method: halving h must divide
    // the error over T = 10 by 2^8 = 256. Measured: 231, 250, 254 as h goes 1 -> 0.125.
    const auto f        = [](double, const OdeState<2>& y) { return OdeState<2>{y[1], -y[0]}; };
    double     previous = 0.0;
    for (const double h : {0.5, 0.25, 0.125}) {
        AdaptiveOptions opt;
        opt.initialStep = h;
        opt.minStep     = h;
        opt.maxStep     = h;
        opt.relTol      = 1e3;
        opt.absTol      = 1e3;
        AdaptiveStats     stats;
        const OdeState<2> y     = integrateRkf78<2>(f, 0.0, {1.0, 0.0}, 10.0, opt, &stats);
        const double      error = std::hypot(y[0] - std::cos(10.0), y[1] + std::sin(10.0));
        EXPECT_EQ(stats.accepted, static_cast<int>(10.0 / h));
        EXPECT_EQ(stats.rejected, 0);
        if (previous > 0.0) {
            EXPECT_GT(previous / error, 225.0) << "h = " << h;
            EXPECT_LT(previous / error, 270.0) << "h = " << h;
        }
        previous = error;
    }
}

TEST(Rkf78, GlobalErrorFollowsTheTolerance) {
    // With relTol = absTol = tol the error after 10 units on the oscillator is 2.3 tol at every
    // tolerance from 1e-7 to 1e-13 (measured), forward and backward. Allowed: 0.3 tol to 10 tol,
    // and at least a factor 50 for two decades of tolerance.
    const auto f        = [](double, const OdeState<2>& y) { return OdeState<2>{y[1], -y[0]}; };
    double     previous = 0.0;
    for (const double tol : {1e-7, 1e-9, 1e-11, 1e-13}) {
        AdaptiveOptions opt;
        opt.relTol              = tol;
        opt.absTol              = tol;
        const OdeState<2> y     = integrateRkf78<2>(f, 0.0, {1.0, 0.0}, 10.0, opt);
        const double      error = std::hypot(y[0] - std::cos(10.0), y[1] + std::sin(10.0));
        EXPECT_GT(error, 0.3 * tol) << "tol = " << tol;
        EXPECT_LT(error, 10.0 * tol) << "tol = " << tol;
        if (previous > 0.0) {
            EXPECT_GT(previous / error, 50.0) << "tol = " << tol;
        }
        previous = error;
        // Backward from the exact state at t = 10 to the start.
        const OdeState<2> back = integrateRkf78<2>(f, 10.0, {std::cos(10.0), -std::sin(10.0)}, 0.0, opt);
        EXPECT_LT(std::hypot(back[0] - 1.0, back[1]), 10.0 * tol) << "tol = " << tol;
    }
}

TEST(Rkf78, GivesUpWhenTheStepUnderflowsOrTheStepsRunOut) {
    const auto f         = [](double, const OdeState<2>& y) { return OdeState<2>{y[1], -y[0]}; };
    const auto failsWith = [&](const AdaptiveOptions& opt, const std::string& message) {
        try {
            (void)integrateRkf78<2>(f, 0.0, OdeState<2>{1.0, 0.0}, 10.0, opt);
            ADD_FAILURE() << "no exception, expected: " << message;
        } catch (const Error& e) {
            EXPECT_NE(std::string(e.what()).find(message), std::string::npos) << e.what();
        }
    };
    // A step of 0.5 has an error of ~1e-10, far above a tolerance of 1e-15, and may not shrink.
    AdaptiveOptions opt;
    opt.relTol  = 1e-15;
    opt.absTol  = 1e-15;
    opt.minStep = 0.5;
    failsWith(opt, "step size underflow");
    // The default tolerance needs about 60 steps for this interval.
    opt          = {};
    opt.maxSteps = 5;
    failsWith(opt, "maximum number of steps exceeded");
}

TEST(NumericalPropagator, PointMassOrbitMatchesKepler) {
    // A point-mass field is the two-body problem, whose solution (propagateKepler) is checked against
    // an extended-precision reference in two_body_test.cpp. One day on a 1600 x 4000 km orbit,
    // forward, backward and at the epoch, several outputs in one call. Measured after a day: 0.57 mm
    // and 6e-7 m/s, the integrator's own error at relTol = 1e-13.
    const StateVector         s0 = elementsToState({8000e3, 0.3, 0.9, 0.2, 0.4, 0.0});
    const NumericalPropagator prop(Time(),
                                   s0,
                                   {std::make_shared<GravityField>(std::vector<double>{})},
                                   {.relTol = 1e-13, .absTol = 1e-9});
    EXPECT_EQ(prop.epoch(), Time());
    const std::vector<double> dts{86400.0, 3600.0, -3600.0, 0.0, 7200.0, -86400.0};
    const auto                many = prop.propagate(dts);
    for (std::size_t k = 0; k < dts.size(); ++k) {
        const StateVector kepler = propagateKepler(s0, dts[k]);
        EXPECT_LT(norm(many[k].r - kepler.r), 1e-3) << "dt = " << dts[k];
        EXPECT_LT(norm(many[k].v - kepler.v), 1e-6) << "dt = " << dts[k];
    }
    EXPECT_EQ(many[3], s0);
}
