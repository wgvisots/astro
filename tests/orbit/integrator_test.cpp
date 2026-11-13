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
#include <vector>

#include <gtest/gtest.h>

#include "astro/core/angle.h"
#include "astro/core/constants.h"
#include "astro/core/error.h"
#include "astro/orbit/elements.h"
#include "astro/orbit/integrators.h"
#include "astro/orbit/numerical.h"

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
