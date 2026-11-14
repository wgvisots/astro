/**********************************************************************
 * orbit_test.cpp
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

#include <cmath>
#include <limits>
#include <memory>
#include <random>
#include <vector>

#include "astro/core/angle.h"
#include "astro/core/error.h"
#include "astro/orbit/elements.h"
#include "astro/orbit/gravity.h"
#include "astro/orbit/integrators.h"
#include "astro/orbit/kepler.h"
#include "astro/orbit/numerical.h"
#include "astro/orbit/secular.h"
#include "astro/orbit/two_body.h"
#include "test_util.h"

using namespace astro;
using astro::test::expectVecNear;

namespace {

constexpr double km = 1000.0;

/// Brouwer's (1959) secular rates of the mean anomaly, argument of perigee and node, to first order
/// in J2 and with the J2^2 and J4 terms, written out in extended precision. Brouwer's small
/// parameters are gamma2' = (J2 / 2) (Re / p)^2 and gamma4' = -(3/8) J4 (Re / p)^4 (p the
/// semi-latus rectum), eta = sqrt(1 - e^2) and theta = cos i; the coefficients are those of the
/// paper's secular terms, as also used by Orekit's BrouwerLyddanePropagator.
SecularRates brouwerRates(double a, double e, double i, const ZonalParameters& zp, bool secondOrder) {
    using L      = long double;
    const L n    = std::sqrt(static_cast<L>(zp.mu) / (static_cast<L>(a) * a * a));
    const L eta2 = 1.0L - static_cast<L>(e) * e;
    const L eta  = std::sqrt(eta2);
    const L q    = static_cast<L>(zp.re) / (static_cast<L>(a) * eta2); // Re / p
    const L g2   = 0.5L * zp.J2 * q * q;
    const L g4   = -0.375L * zp.J4 * q * q * q * q;
    const L th   = std::cos(static_cast<L>(i));
    const L th2  = th * th;
    L       lDot = 1.0L + 1.5L * g2 * eta * (3.0L * th2 - 1.0L);
    L       gDot = 1.5L * g2 * (5.0L * th2 - 1.0L);
    L       hDot = -3.0L * g2 * th;
    if (secondOrder) {
        const L th4  = th2 * th2;
        lDot        += 3.0L / 32.0L * g2 * g2 * eta
                  * (-15.0L + 16.0L * eta + 25.0L * eta2 + (30.0L - 96.0L * eta - 90.0L * eta2) * th2
                     + (105.0L + 144.0L * eta + 25.0L * eta2) * th4)
              + 15.0L / 16.0L * g4 * eta * e * e * (3.0L - 30.0L * th2 + 35.0L * th4);
        gDot += 3.0L / 32.0L * g2 * g2
                  * (-35.0L + 24.0L * eta + 25.0L * eta2 + (90.0L - 192.0L * eta - 126.0L * eta2) * th2
                     + (385.0L + 360.0L * eta + 45.0L * eta2) * th4)
              + 5.0L / 16.0L * g4
                    * (21.0L - 9.0L * eta2 + (-270.0L + 126.0L * eta2) * th2 + (385.0L - 189.0L * eta2) * th4);
        hDot += 3.0L / 8.0L * g2 * g2
                  * ((-5.0L + 12.0L * eta + 9.0L * eta2) * th - (35.0L + 36.0L * eta + 5.0L * eta2) * th2 * th)
              + 5.0L / 4.0L * g4 * (5.0L - 3.0L * eta2) * th * (3.0L - 7.0L * th2);
    }
    return {static_cast<double>(n * hDot), static_cast<double>(n * gDot), static_cast<double>(n * lDot)};
}

} // namespace

TEST(Kepler, AnomalyRoundTrips) {
    std::mt19937                           rng(5);
    std::uniform_real_distribution<double> ang(-3.1, 3.1);
    for (double e : {0.0, 1e-8, 0.1, 0.5, 0.9, 0.99, 0.999}) {
        for (int k = 0; k < 200; ++k) {
            const double M = ang(rng);
            const double E = meanToEccentricAnomaly(M, e);
            EXPECT_NEAR(eccentricToMeanAnomaly(E, e), M, 1e-13) << e;
            const double nu = eccentricToTrueAnomaly(E, e);
            EXPECT_NEAR(wrapPi(trueToEccentricAnomaly(nu, e) - E), 0.0, 1e-12) << e;
            EXPECT_NEAR(wrapPi(trueToMeanAnomaly(meanToTrueAnomaly(M, e), e) - M), 0.0, 1e-11) << e;
        }
    }
    // Multiple revolutions are preserved.
    EXPECT_NEAR(meanToEccentricAnomaly(4.0 * constants::pi + 0.3, 0.2) - meanToEccentricAnomaly(0.3, 0.2),
                4.0 * constants::pi,
                1e-12);
    for (double e : {1.01, 1.5, 3.0, 10.0}) {
        for (double M : {-50.0, -1.0, -1e-3, 0.0, 0.5, 5.0, 200.0}) {
            const double F = meanToHyperbolicAnomaly(M, e);
            EXPECT_NEAR(hyperbolicToMeanAnomaly(F, e), M, 1e-10 * (1.0 + std::abs(M))) << e << " " << M;
            EXPECT_NEAR(trueToHyperbolicAnomaly(hyperbolicToTrueAnomaly(F, e), e), F, 1e-9 * (1.0 + std::abs(F)));
        }
    }
    for (double M : {-10.0, -0.1, 0.0, 0.7, 30.0}) {
        EXPECT_NEAR(parabolicToMeanAnomaly(meanToParabolicAnomaly(M)), M, 1e-12 * (1.0 + std::abs(M)));
    }
    EXPECT_THROW((void)meanToEccentricAnomaly(1.0, 1.0), Error);
}

TEST(Elements, ValladoExample2_5) {
    // Vallado, Fundamentals of Astrodynamics, Example 2-5 (RV2COE).
    const StateVector s{
        {6524.834 * km, 6862.875 * km, 6448.296 * km },
        {4.901327 * km, 5.533756 * km, -1.976341 * km}
    };
    const KeplerianElements el = stateToElements(s, 398600.4418e9);
    EXPECT_NEAR(el.semiLatusRectum() / km, 11067.790, 0.01);
    EXPECT_NEAR(el.a / km, 36127.343, 0.01);
    EXPECT_NEAR(el.e, 0.832853, 1e-6);
    EXPECT_NEAR(rad2deg(el.i), 87.870, 1e-3);
    EXPECT_NEAR(rad2deg(el.raan), 227.89, 0.01);
    EXPECT_NEAR(rad2deg(el.argp), 53.38, 0.01);
    EXPECT_NEAR(rad2deg(el.nu), 92.335, 1e-3);
    const StateVector back = elementsToState(el, 398600.4418e9);
    expectVecNear(back.r, s.r, 1e-6);
    expectVecNear(back.v, s.v, 1e-9);
}

TEST(Elements, RandomRoundTripsAndSpecialCases) {
    std::mt19937                           rng(9);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    for (int k = 0; k < 500; ++k) {
        KeplerianElements el;
        el.e                 = k % 5 == 0 ? 1.2 + 2.0 * u(rng) : 0.95 * u(rng);
        el.a                 = el.e < 1.0 ? (6600.0 + 40000.0 * u(rng)) * km : -(7000.0 + 20000.0 * u(rng)) * km;
        el.i                 = 0.01 + 3.1 * u(rng);
        el.raan              = constants::twoPi * u(rng);
        el.argp              = constants::twoPi * u(rng);
        const double nuMax   = el.e < 1.0 ? constants::pi : 0.9 * std::acos(-1.0 / el.e);
        el.nu                = (2.0 * u(rng) - 1.0) * nuMax;
        const StateVector s  = elementsToState(el);
        const StateVector s2 = elementsToState(stateToElements(s));
        expectVecNear(s2.r, s.r, 1e-6 * norm(s.r) * 1e-3);
        expectVecNear(s2.v, s.v, 1e-9 * norm(s.v));
        // Equinoctial round trip (elliptic, prograde-ish).
        if (el.e < 1.0) {
            const StateVector s3 = equinoctialToState(stateToEquinoctial(s));
            expectVecNear(s3.r, s.r, 1e-6);
            expectVecNear(s3.v, s.v, 1e-9);
            const EquinoctialElements q  = keplerianToEquinoctial(el);
            const StateVector         s4 = equinoctialToState(q);
            expectVecNear(s4.r, s.r, 1e-5);
        }
    }
    // Circular inclined: argp = 0, nu = argument of latitude.
    const StateVector       circ = elementsToState({7000 * km, 0.0, 0.5, 1.0, 0.0, 0.7});
    const KeplerianElements c    = stateToElements(circ);
    EXPECT_NEAR(c.e, 0.0, 1e-12);
    EXPECT_DOUBLE_EQ(c.argp, 0.0);
    EXPECT_NEAR(c.nu, 0.7, 1e-10);
    EXPECT_NEAR(c.raan, 1.0, 1e-12);
    // Circular equatorial: true longitude.
    const StateVector ce = elementsToState({7000 * km, 0.0, 0.0, 0.0, 0.0, 2.0});
    EXPECT_NEAR(stateToElements(ce).nu, 2.0, 1e-10);
    // Elliptic equatorial: longitude of periapsis.
    const StateVector       ee = elementsToState({9000 * km, 0.2, 0.0, 0.0, 1.1, 0.4});
    const KeplerianElements e2 = stateToElements(ee);
    EXPECT_NEAR(e2.argp, 1.1, 1e-10);
    EXPECT_NEAR(e2.nu, 0.4, 1e-10);
    EXPECT_NEAR(orbitalPeriod(42164.17 * km) / 86164.09, 1.0, 1e-4); // GEO period = sidereal day
    EXPECT_NEAR(semiMajorAxisFromMeanMotion(meanMotion(7000 * km)), 7000 * km, 1e-6);
}

TEST(TwoBody, ValladoExample2_4) {
    // Vallado Example 2-4 (KEPLER): 40 minutes.
    const double      mu = 398600.4418e9;
    const StateVector s0{
        {1131.340 * km, -2282.343 * km, 6672.423 * km},
        {-5.64305 * km, 4.30333 * km,   2.42879 * km }
    };
    const StateVector s = propagateKepler(s0, 40.0 * 60.0, mu);
    expectVecNear(s.r, {-4219.7527 * km, 4363.0292 * km, -3958.7666 * km}, 1.0);
    expectVecNear(s.v, {3.689866 * km, -1.916735 * km, -6.112511 * km}, 1e-3);
}

TEST(TwoBody, AgreesWithElementPropagation) {
    std::mt19937                           rng(17);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    for (int k = 0; k < 200; ++k) {
        KeplerianElements el{(7000.0 + 30000.0 * u(rng)) * km,
                             0.9 * u(rng),
                             3.0 * u(rng),
                             6.0 * u(rng),
                             6.0 * u(rng),
                             6.0 * u(rng)};
        const double      dt       = (u(rng) - 0.3) * 5.0 * orbitalPeriod(el.a);
        const StateVector s0       = elementsToState(el);
        KeplerianElements e1       = el;
        e1.nu                      = meanToTrueAnomaly(el.meanAnomaly() + meanMotion(el.a) * dt, el.e);
        const StateVector expected = elementsToState(e1);
        const StateVector s        = propagateKepler(s0, dt);
        expectVecNear(s.r, expected.r, 1e-3 + 1e-10 * norm(s.r));
        expectVecNear(s.v, expected.v, 1e-6);
    }
    // Hyperbolic: energy and angular momentum conserved; round trip.
    const StateVector h0{
        {7000 * km, 0,         0       },
        {0,         12.0 * km, 1.0 * km}
    };
    const StateVector h1 = propagateKepler(h0, 20000.0);
    EXPECT_NEAR(specificEnergy(h1), specificEnergy(h0), 1e-6 * std::abs(specificEnergy(h0)));
    expectVecNear(cross(h1.r, h1.v), cross(h0.r, h0.v), 1e-3 * norm(cross(h0.r, h0.v)) * 1e-6);
    const StateVector h2 = propagateKepler(h1, -20000.0);
    expectVecNear(h2.r, h0.r, 1e-3);
    // Near-parabolic.
    const double      vesc = std::sqrt(2.0 * constants::gmEarth / (7000 * km));
    const StateVector p0{
        {7000 * km, 0,    0},
        {0,         vesc, 0}
    };
    const StateVector p1 = propagateKepler(p0, 3600.0);
    expectVecNear(propagateKepler(p1, -3600.0).r, p0.r, 1e-2);
}

TEST(Gravity, ZonalMatchesClosedFormAndPotentialGradient) {
    const Vec3 r{5000 * km, -3000 * km, 4000 * km};
    expectVecNear(zonalAcceleration(r, {0, 0, constants::J2}), j2Acceleration(r), 1e-15);
    // Gradient of U = -(mu/r) sum Jn (R/r)^n Pn(z/r) by central differences.
    const std::vector<double> J{0, 0, constants::J2, constants::J3, constants::J4, constants::J5, constants::J6};
    auto                      U = [&](const Vec3& p) {
        const double rm = norm(p), s = p.z / rm, q = constants::earthRadiusEgm2008 / rm;
        double       P0 = 1.0, P1 = s, sum = 0.0, qn = q;
        for (int n = 2; n <= 6; ++n) {
            const double P2  = ((2.0 * n - 1.0) * s * P1 - (n - 1.0) * P0) / n;
            qn              *= q;
            sum             += J[static_cast<std::size_t>(n)] * qn * P2;
            P0               = P1;
            P1               = P2;
        }
        return -constants::gmEarth / rm * sum;
    };
    const double h = 1.0;
    const Vec3   grad{(U(r + Vec3{h, 0, 0}) - U(r - Vec3{h, 0, 0})) / (2 * h),
                    (U(r + Vec3{0, h, 0}) - U(r - Vec3{0, h, 0})) / (2 * h),
                    (U(r + Vec3{0, 0, h}) - U(r - Vec3{0, 0, h})) / (2 * h)};
    expectVecNear(zonalAcceleration(r, J), grad, 1e-10); // finite-difference roundoff ~1e-11
    // Third body: zero at the origin, tidal elsewhere.
    expectVecNear(thirdBodyAcceleration({}, {384400 * km, 0, 0}, constants::gmMoon), {}, 1e-20);
}

TEST(Integrators, OrderAndTwoBodyAccuracy) {
    // Harmonic oscillator y'' = -y.
    auto              f  = [](double, const OdeState<2>& y) { return OdeState<2>{y[1], -y[0]}; };
    const OdeState<2> r4 = integrateRk4<2>(f, 0.0, {1.0, 0.0}, 10.0, 0.01);
    EXPECT_NEAR(r4[0], std::cos(10.0), 1e-8);
    AdaptiveStats     st;
    const OdeState<2> r8 = integrateRkf78<2>(f, 0.0, {1.0, 0.0}, 10.0, {.relTol = 1e-13, .absTol = 1e-14}, &st);
    EXPECT_NEAR(r8[0], std::cos(10.0), 1e-12);
    EXPECT_NEAR(r8[1], -std::sin(10.0), 1e-12);
    EXPECT_GT(st.accepted, 0);
    const OdeState<2> back = integrateRkf78<2>(f, 10.0, r8, 0.0, {.relTol = 1e-13, .absTol = 1e-14});
    EXPECT_NEAR(back[0], 1.0, 1e-11);
    // Two-body orbit over one day vs Kepler.
    const StateVector         s0 = elementsToState({8000 * km, 0.3, 0.9, 0.2, 0.4, 0.0});
    const NumericalPropagator prop(Time(),
                                   s0,
                                   {std::make_shared<GravityField>(std::vector<double>{})},
                                   {.relTol = 1e-13, .absTol = 1e-9});
    const StateVector         s = prop.propagate(86400.0);
    const StateVector         k = propagateKepler(s0, 86400.0);
    expectVecNear(s.r, k.r, 1e-3);
    expectVecNear(s.v, k.v, 1e-5);
    // Multiple outputs in both directions in one call.
    const auto many = prop.propagate(std::vector<double>{3600.0, -3600.0, 0.0, 7200.0});
    expectVecNear(many[1].r, propagateKepler(s0, -3600.0).r, 1e-3);
    expectVecNear(many[2].r, s0.r, 1e-9);
    expectVecNear(many[3].r, propagateKepler(s0, 7200.0).r, 1e-3);
}

TEST(Secular, RatesMatchAnAveragedNumericalIntegration) {
    // The truth: a numerical integration with the zonal field alone (J2, then J2 and J4) over 20 days,
    // about 270 revolutions, sampled every minute. The secular rates are the least-squares slopes of
    // the osculating angles, and the mean a, e and i are the averages of the osculating ones (the
    // short-period terms of Brouwer's theory average to zero over the mean anomaly). What remains is
    // of order J2^3 and the long-period terms, which leave the node and the argument of latitude
    // (perigee plus mean anomaly) within 1e-6 of their rates and the perigee within 2e-4 at e = 0.05.
    // The second-order terms (J2^2, J4) are 1e-3 of the first-order ones, so this checks them too.
    const KeplerianElements el{7500 * km, 0.05, deg2rad(50.0), 0.3, 0.5, 0.0};
    const double            span = 20.0 * 86400.0;
    const std::size_t       n    = 20 * 1440;
    for (const double J4 : {0.0, constants::J4}) {
        const ZonalParameters     zp{constants::gmEarth, constants::earthRadiusEgm2008, constants::J2, J4};
        const NumericalPropagator prop(
            Time(),
            elementsToState(el),
            {std::make_shared<GravityField>(std::vector<double>{0, 0, constants::J2, 0, J4})});
        const double                   step   = span / static_cast<double>(n);
        const auto                     states = prop.propagateGrid(0.0, step, n + 1);
        std::vector<KeplerianElements> osc;
        double                         am = 0.0, em = 0.0, im = 0.0;
        for (const StateVector& s : states) {
            osc.push_back(stateToElements(s));
            am += osc.back().a;
            em += osc.back().e;
            im += osc.back().i;
        }
        const double count = static_cast<double>(osc.size());
        // Least-squares slope of an unwrapped angle against time.
        const auto slope = [&](auto angle) {
            double prev = angle(osc[0]), off = 0.0, sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;
            for (std::size_t k = 0; k < osc.size(); ++k) {
                const double a  = angle(osc[k]);
                off            += a - prev > constants::pi  ? -constants::twoPi
                                : a - prev < -constants::pi ? constants::twoPi
                                                            : 0.0;
                prev            = a;
                const double x = step * static_cast<double>(k), y = a + off;
                sx  += x;
                sy  += y;
                sxx += x * x;
                sxy += x * y;
            }
            return (sxy - sx * sy / count) / (sxx - sx * sx / count);
        };
        const double node     = slope([](const KeplerianElements& e) { return e.raan; });
        const double perigee  = slope([](const KeplerianElements& e) { return e.argp; });
        const double latitude = slope([](const KeplerianElements& e) { return e.argp + e.meanAnomaly(); });

        const SecularRates r = secularRates(am / count, em / count, im / count, zp, SecularOrder::J2SquaredJ4);
        // Measured 6e-7 (node), 1.9e-4 (perigee), 1.1e-7 (argument of latitude) for J2 alone.
        EXPECT_NEAR(node / r.raanDot, 1.0, 3e-6) << "J4 = " << J4;
        EXPECT_NEAR(perigee / r.argpDot, 1.0, 1e-3) << "J4 = " << J4;
        EXPECT_NEAR(latitude / (r.argpDot + r.meanAnomalyDot), 1.0, 1e-6) << "J4 = " << J4;
        // First order alone misses the node by 1e-3 of its rate: the check resolves the second order.
        const SecularRates r1 = secularRates(am / count, em / count, im / count, zp, SecularOrder::J2);
        EXPECT_GT(std::abs(node / r1.raanDot - 1.0), 1e-4) << "J4 = " << J4;
    }
}

TEST(Secular, BrouwerReducesToSgp4FormForCircularOrbits) {
    // SGP4's secular rates (Kozai/Brouwer) for e -> 0, written out from sgp4init.
    const double          a = 7000 * km, i = deg2rad(63.0);
    const ZonalParameters zp;
    const double          n = std::sqrt(zp.mu / (a * a * a));
    const double          p = a / zp.re, pinvsq = 1.0 / (p * p), th = std::cos(i), th2 = th * th, th4 = th2 * th2;
    const double          temp1 = 1.5 * zp.J2 * pinvsq * n, temp2 = 0.5 * temp1 * zp.J2 * pinvsq,
                 temp3   = -0.46875 * zp.J4 * pinvsq * pinvsq * n;
    const double mdot    = n + 0.5 * temp1 * (3.0 * th2 - 1.0) + 0.0625 * temp2 * (13.0 - 78.0 * th2 + 137.0 * th4);
    const double argpdot = -0.5 * temp1 * (1.0 - 5.0 * th2) + 0.0625 * temp2 * (7.0 - 114.0 * th2 + 395.0 * th4)
                         + temp3 * (3.0 - 36.0 * th2 + 49.0 * th4);
    const double nodedot = -temp1 * th + (0.5 * temp2 * (4.0 - 19.0 * th2) + 2.0 * temp3 * (3.0 - 7.0 * th2)) * th;
    const SecularRates r = secularRates(a, 0.0, i, zp);
    EXPECT_NEAR(r.meanAnomalyDot, mdot, 1e-15 * n);
    EXPECT_NEAR(r.argpDot, argpdot, 1e-15 * n);
    EXPECT_NEAR(r.raanDot, nodedot, 1e-15 * n);
    // J2^2 / J4 corrections are ~1e-3 of the J2 terms.
    const SecularRates r1 = secularRates(a, 0.0, i, zp, SecularOrder::J2);
    EXPECT_GT(std::abs(r.raanDot - r1.raanDot), 1e-4 * std::abs(r1.raanDot));
    EXPECT_LT(std::abs(r.raanDot - r1.raanDot), 1e-2 * std::abs(r1.raanDot));
}

TEST(Secular, MissionDesignHelpers) {
    const ZonalParameters zp{constants::gmEarth, constants::wgs84A, 1.08263e-3, 0.0};
    // Sun-synchronous at 700 km: ~98.19 deg.
    EXPECT_NEAR(rad2deg(sunSynchronousInclination(constants::wgs84A + 700 * km, 0.0, zp)), 98.19, 0.02);
    EXPECT_THROW((void)sunSynchronousInclination(constants::wgs84A + 7000 * km, 0.0, zp), Error);
    // Landsat-8: 233 revolutions in 16 days, i = 98.2 deg; published semi-major axis 7077.7 km
    // (the "705 km" nominal altitude is measured above the ellipsoid at mid latitudes).
    const double a = repeatGroundTrackSemiMajorAxis(233, 16, 0.0, deg2rad(98.2), zp);
    EXPECT_NEAR(a / km, 7077.7, 1.0);
}

TEST(Secular, CriticalInclinationIsWhereTheJ2ApsidalRateVanishes) {
    // acos(1 / sqrt(5)) = atan(2): the root of 5 cos^2 i - 1 = 0.
    EXPECT_NEAR(kCriticalInclination, std::atan(2.0), 2.0 * std::numeric_limits<double>::epsilon());
    const SecularRates r = secularRates(7000 * km, 0.1, kCriticalInclination, {}, SecularOrder::J2);
    EXPECT_LT(std::abs(r.argpDot), 1e-15 * r.meanAnomalyDot);
}

TEST(Secular, RatesMatchBrouwersFormulasForEccentricOrbits) {
    // Both orders, from circular to e = 0.95 and from equatorial prograde to equatorial retrograde
    // (i = 180 deg, where the node moves eastward at the rate of i = 0 reversed).
    // Tolerance: 16 rounding errors of the largest term in each rate (n for the mean anomaly, the
    // first-order term for the node and perigee). Measured: 2 eps.
    const ZonalParameters zp;
    const double          eps = std::numeric_limits<double>::epsilon();
    for (double a : {7000 * km, 26600 * km}) {
        for (double e : {0.0, 0.01, 0.3, 0.7, 0.95}) {
            if (a * (1.0 - e) < 6400 * km) {
                continue;
            }
            for (double i : {0.0, 0.3, kCriticalInclination, 1.5, 2.5, constants::pi}) {
                for (SecularOrder order : {SecularOrder::J2, SecularOrder::J2SquaredJ4}) {
                    const SecularRates ref   = brouwerRates(a, e, i, zp, order == SecularOrder::J2SquaredJ4);
                    const SecularRates first = brouwerRates(a, e, i, zp, false);
                    const SecularRates r     = secularRates(a, e, i, zp, order);
                    const double       scale = std::abs(first.raanDot) + std::abs(first.argpDot);
                    EXPECT_NEAR(r.meanAnomalyDot, ref.meanAnomalyDot, 16.0 * eps * ref.meanAnomalyDot)
                        << order << " a=" << a << " e=" << e << " i=" << i;
                    EXPECT_NEAR(r.argpDot, ref.argpDot, 16.0 * eps * scale)
                        << order << " a=" << a << " e=" << e << " i=" << i;
                    EXPECT_NEAR(r.raanDot, ref.raanDot, 16.0 * eps * scale)
                        << order << " a=" << a << " e=" << e << " i=" << i;
                }
            }
        }
    }
    // The second-order terms matter at e != 0: at e = 0.7 they move the node by about 1e-3 of its
    // rate, far above the tolerance.
    const SecularRates j2 = secularRates(26600 * km, 0.7, 1.1, zp, SecularOrder::J2);
    const SecularRates j4 = secularRates(26600 * km, 0.7, 1.1, zp, SecularOrder::J2SquaredJ4);
    EXPECT_GT(std::abs(j4.raanDot - j2.raanDot), 1e-4 * std::abs(j2.raanDot));
}

TEST(Secular, RetrogradeEquatorialNodeMovesEastwardAsFastAsThePrograde) {
    // cos(pi) = -1: the J2 node rate at i = 180 deg is exactly minus the rate at i = 0, the perigee
    // and mean anomaly rates are the same. Not a special case in the formula, but in the elements:
    // a state there has RAAN 0 and its longitude of periapsis in argp, measured clockwise. The
    // propagator must still give a continuous answer (compare i = pi - 1e-7 below).
    const ZonalParameters zp;
    const SecularRates    pro   = secularRates(7000 * km, 0.01, 0.0, zp);
    const SecularRates    retro = secularRates(7000 * km, 0.01, constants::pi, zp);
    EXPECT_GT(retro.raanDot, 0.0);
    EXPECT_NEAR(retro.raanDot, -pro.raanDot, 1e-15 * std::abs(pro.raanDot));
    EXPECT_NEAR(retro.argpDot, pro.argpDot, 1e-15 * std::abs(pro.argpDot));
    EXPECT_NEAR(retro.meanAnomalyDot, pro.meanAnomalyDot, 1e-15 * pro.meanAnomalyDot);
    // States on an orbit at exactly 180 deg and on one 1e-7 rad away separate by about 1e-7 of the
    // radius after a day (0.7 m); a special case handled the wrong way (the perigee drifting the
    // wrong way round, say) would put them kilometres apart.
    const double      day   = 86400.0;
    const StateVector exact = elementsToState({7000 * km, 0.01, constants::pi, 0.5, 1.0, 2.0});
    const StateVector near  = elementsToState({7000 * km, 0.01, constants::pi - 1e-7, 0.5, 1.0, 2.0});
    EXPECT_LT(norm(exact.r - near.r), 1.0); // m, at the epoch
    EXPECT_LT(norm(SecularPropagator(exact).propagate(day).r - SecularPropagator(near).propagate(day).r), 2.0);
}

TEST(Numerical, J2J4PropagatorConservesEnergyAndFollowsSecularRates) {
    const KeplerianElements   el{7200 * km, 0.08, deg2rad(40.0), 0.1, 1.0, 0.0};
    const StateVector         s0   = elementsToState(el);
    const NumericalPropagator prop = NumericalPropagator::j2j4(Time(), s0);
    const std::vector<double> J{0, 0, constants::J2, 0, constants::J4};
    auto                      energy = [&](const StateVector& s) {
        const double rm = norm(s.r), sn = s.r.z / rm, q = constants::earthRadiusEgm2008 / rm;
        const double P2 = 0.5 * (3 * sn * sn - 1), P4 = (35 * std::pow(sn, 4) - 30 * sn * sn + 3) / 8;
        const double U =
            -constants::gmEarth / rm * (1.0 - constants::J2 * q * q * P2 - constants::J4 * std::pow(q, 4) * P4);
        return 0.5 * norm2(s.v) + U;
    };
    const StateVector s = prop.propagate(10.0 * 86400.0);
    EXPECT_NEAR(energy(s) / energy(s0), 1.0, 1e-10);
    // RAAN drift over 10 days agrees with the J2 + J2^2 + J4 secular rate to ~0.3%.
    const double       draan = wrapPi(stateToElements(s).raan - el.raan);
    const SecularRates r     = secularRates(el.a, el.e, el.i);
    EXPECT_NEAR(draan / wrapPi(r.raanDot * 10.0 * 86400.0), 1.0, 3e-3);
}
