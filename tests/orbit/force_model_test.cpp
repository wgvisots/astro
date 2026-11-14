/**********************************************************************
 * force_model_test.cpp
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

// The force models of the numerical propagator, each against what it is defined to be: the central
// field in a rotated body frame, the tidal pull of a third body at its ephemeris position at
// epoch + t, and an unmodeled acceleration that is linear in time (whose trajectory, alone, is a
// cubic the integrator must reproduce exactly).

#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <memory>
#include <vector>

#include "astro/core/angle.h"
#include "astro/core/constants.h"
#include "astro/core/error.h"
#include "astro/ephemeris/ephemeris.h"
#include "astro/orbit/elements.h"
#include "astro/orbit/gravity.h"
#include "astro/orbit/numerical.h"
#include "astro/time/time.h"
#include "test_util.h"

using namespace astro;
using astro::test::expectMatNear;
using astro::test::expectVecNear;

namespace {

constexpr double km   = 1000.0;
constexpr double kEps = std::numeric_limits<double>::epsilon();

/// Rotation by `angle` about x (a frame turned by -angle).
Mat3 rotationX(double angle) {
    const double c = std::cos(angle), s = std::sin(angle);
    Mat3         R = Mat3::identity();
    R(1, 1)        = c;
    R(1, 2)        = s;
    R(2, 1)        = -s;
    R(2, 2)        = c;
    return R;
}

/// An ephemeris in which the Moon moves on a straight line, s(t) = s0 + w (t - epoch), so that the
/// position a force model asks for at any time is known exactly. Any other body throws.
class LinearMoon final : public Ephemeris {
public:
    LinearMoon(const Time& epoch, const Vec3& s0, const Vec3& w)
        : epoch_(epoch)
        , s0_(s0)
        , w_(w) {}
    [[nodiscard]] StateVector state(Body target, Body center, const Time& t) const override {
        if (target != Body::Moon || center != Body::Earth) {
            throw Error("LinearMoon: only the Moon relative to the Earth");
        }
        return {s0_ + w_ * (t - epoch_), w_};
    }

private:
    Time epoch_;
    Vec3 s0_, w_;
};

} // namespace

TEST(CentralBody, NamedValuesAndZonalLookup) {
    // The default is the Earth with EGM2008 J2, J3 and J4; earthJ2J4 drops J3 only; a point mass
    // has no zonal terms. zonal(n) reads past the end of the list as 0.
    const CentralBody def, earth = CentralBody::earth(), j2j4 = CentralBody::earthJ2J4();
    const CentralBody mass = CentralBody::pointMass(4.9028e12);
    for (const CentralBody* b : {&def, &earth, &j2j4}) {
        EXPECT_EQ(b->mu, constants::gmEarth);
        EXPECT_EQ(b->re, constants::earthRadiusEgm2008);
        EXPECT_EQ(b->zonal(2), constants::J2);
        EXPECT_EQ(b->zonal(4), constants::J4);
        EXPECT_EQ(b->zonal(5), 0.0);
    }
    EXPECT_EQ(def.J, earth.J);
    EXPECT_EQ(earth.zonal(3), constants::J3);
    EXPECT_EQ(j2j4.zonal(3), 0.0);
    EXPECT_EQ(mass.mu, 4.9028e12);
    EXPECT_TRUE(mass.J.empty());
    EXPECT_EQ(mass.zonal(2), 0.0);
    EXPECT_EQ(CentralBody::pointMass().mu, constants::gmEarth);
    // A point-mass body gives the point-mass field and no zonal acceleration.
    const Vec3 r{7000 * km, -1000 * km, 2000 * km};
    const Vec3 a = pointMassAcceleration(r, mass.mu);
    expectVecNear(GravityField(mass).acceleration(0.0, r, {}), a, 4.0 * kEps * norm(a), "point mass");
    EXPECT_EQ(zonalAcceleration(r, mass), Vec3{});
}

TEST(GravityField, RotatedBodyFrameTurnsTheFieldWithIt) {
    // With toBodyFrame = R the field is the body-frame field at R r, turned back: a = R^T a_b(R r),
    // jerk = R^T j_b(R r, R v), gradient = R^T G_b(R r) R. The body-frame truths are the free
    // functions (checked against closed forms and the potential elsewhere). A non-default mu and
    // reference radius are passed through, and body() reports the field's.
    const double      mu = 4.9028e12, re = 1737.4 * km; // a Moon-sized body
    const CentralBody body{
        mu,
        re,
        {0.0, 0.0, 1.1e-3, -2.5e-6, -1.6e-6}
    };
    const Mat3         R = rotationX(0.4);
    const GravityField field(body, R);
    EXPECT_EQ(field.body().mu, mu);
    EXPECT_EQ(field.body().re, re);
    EXPECT_EQ(GravityField().body().mu, constants::gmEarth);
    for (const Vec3 r : {
             Vec3{2000 * km, 300 * km,  -400 * km},
             Vec3{-500 * km, 1200 * km, 1500 * km},
             Vec3{100 * km,  -100 * km, 1900 * km}
    }) {
        const Vec3 v{0.3 * km, -1.1 * km, 0.7 * km};
        const Vec3 rb = R * r, vb = R * v;
        const Vec3 a = R.transpose() * (pointMassAcceleration(rb, mu) + zonalAcceleration(rb, body));
        const Vec3 j = R.transpose() * (pointMassJerk(rb, vb, mu) + zonalJerk(rb, vb, body));
        const Mat3 G = R.transpose() * (pointMassGradient(rb, mu) + zonalGradient(rb, body)) * R;
        // A few rounding errors of the point-mass term, which dominates (measured 1.8 eps).
        expectVecNear(field.acceleration(0.0, r, v), a, 8.0 * kEps * norm(a), "acceleration");
        expectVecNear(field.jerk(0.0, r, v, a), j, 8.0 * kEps * norm(j), "jerk");
        Mat3       dadr, dadv;
        const Vec3 ap = field.partials(0.0, r, v, dadr, dadv);
        expectVecNear(ap, a, 8.0 * kEps * norm(a), "partials' acceleration");
        expectMatNear(dadr, G, 8.0 * kEps * mu / std::pow(norm(r), 3), "gradient");
        expectMatNear(dadv, Mat3{}, 0.0, "no velocity dependence");
    }
}

TEST(ThirdBodyForce, IsTheTidalPullOfTheBodyAtEpochPlusT) {
    // a = mu ((s - r) / |s - r|^3 - s / |s|^3), with s the body's geocentric position at epoch + t
    // (the second term is the body's pull on the Earth), and da/dr = mu / |d|^3 (3 dhat dhat^T - I)
    // with d = s - r. Written out here; the ephemeris is a Moon on a straight line, so s(t) is exact.
    const Time           epoch = Time::fromIso8601("2026-11-02T00:00:00Z");
    const Vec3           s0{384400 * km, 20000 * km, -15000 * km}, w{-300.0, 900.0, 100.0};
    const double         mu = constants::gmMoon;
    const ThirdBodyForce force(Body::Moon, mu, std::make_shared<LinearMoon>(epoch, s0, w), epoch);
    for (const double t : {0.0, 3600.0, -86400.0}) {
        const Vec3 s = s0 + w * t;
        for (const Vec3 r : {
                 Vec3{7000 * km,  0,         0         },
                 Vec3{-3000 * km, 5000 * km, 4000 * km },
                 Vec3{0,          0,         42164 * km}
        }) {
            const Vec3   d  = s - r;
            const double dm = norm(d), sm = norm(s);
            const Vec3   a = d * (mu / (dm * dm * dm)) - s * (mu / (sm * sm * sm));
            // The two terms nearly cancel (the tide is 1e-4 of each in LEO): rounding of the terms
            // (measured 0.8 eps mu / |s|^2).
            const double tol = 8.0 * kEps * mu / (sm * sm);
            expectVecNear(force.acceleration(t, r, Vec3{}), a, tol, "acceleration");
            Mat3       dadr, dadv;
            const Vec3 ap = force.partials(t, r, Vec3{}, dadr, dadv);
            expectVecNear(ap, a, tol, "partials' acceleration");
            Mat3 G;
            for (std::size_t i = 0; i < 3; ++i) {
                for (std::size_t k = 0; k < 3; ++k) {
                    G(i, k) = mu / (dm * dm * dm) * (3.0 * d[i] * d[k] / (dm * dm) - (i == k ? 1.0 : 0.0));
                }
            }
            expectMatNear(dadr, G, 8.0 * kEps * mu / (dm * dm * dm), "gradient");
            expectMatNear(dadv, Mat3{}, 0.0, "no velocity dependence");
        }
    }
    // Zero at the geocentre: the pull on the satellite is the pull on the Earth.
    expectVecNear(force.acceleration(0.0, Vec3{}, Vec3{}), Vec3{}, 0.0);
    // An ephemeris error passes through.
    const ThirdBodyForce sun(Body::Sun, constants::gmSun, std::make_shared<LinearMoon>(epoch, s0, w), epoch);
    EXPECT_THROW((void)sun.acceleration(0.0, s0, Vec3{}), Error);
}

TEST(UnmodeledAcceleration, IsLinearInTimeAndItsTrajectoryAloneIsACubic) {
    const Vec3                  a0{1e-3, -2e-3, 5e-4}, j0{1e-6, 3e-6, -2e-6};
    const UnmodeledAcceleration force(a0, j0);
    const Vec3                  r{7000 * km, 0.0, 0.0}, v{0.0, 7.5 * km, 0.0};
    for (const double t : {0.0, 100.0, -250.0}) {
        expectVecNear(force.acceleration(t, r, v), a0 + j0 * t, 4.0 * kEps * norm(a0 + j0 * t));
        expectVecNear(force.jerk(t, r, v, Vec3{1.0, 2.0, 3.0}), j0, 0.0);
        Mat3 dadr, dadv;
        expectVecNear(force.partials(t, r, v, dadr, dadv), a0 + j0 * t, 4.0 * kEps * norm(a0 + j0 * t));
        expectMatNear(dadr, Mat3{}, 0.0);
        expectMatNear(dadv, Mat3{}, 0.0);
    }
    // With no other force the motion is r0 + v0 t + a0 t^2 / 2 + j0 t^3 / 6, a cubic, which the
    // eighth-order integrator reproduces to rounding (measured 1.2e-16 relative).
    const NumericalPropagator prop(Time(), {r, v}, {std::make_shared<UnmodeledAcceleration>(a0, j0)});
    for (const double t : {600.0, -1800.0, 7200.0}) {
        const StateVector s = prop.propagate(t);
        const Vec3        x = r + v * t + a0 * (0.5 * t * t) + j0 * (t * t * t / 6.0);
        const Vec3        u = v + a0 * t + j0 * (0.5 * t * t);
        expectVecNear(s.r, x, 1e-15 * norm(x), "position");
        expectVecNear(s.v, u, 1e-15 * norm(u), "velocity");
    }
}

TEST(NumericalPropagator, RejectsAnEmptyListOfForces) {
    const StateVector s0 = elementsToState({7000 * km, 0.001, degToRad(51.6), 0.3, 0.6, 0.9});
    EXPECT_THROW(NumericalPropagator(Time(), s0, {}), Error);
    EXPECT_THROW(NumericalPropagator(Time(), StateVector12{s0.r, s0.v, {}, {}}, {}), Error);
}
