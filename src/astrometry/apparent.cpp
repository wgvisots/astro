/**********************************************************************
 * apparent.cpp
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

#include "astro/astrometry/apparent.h"

#include <algorithm>
#include <cmath>

#include "astro/core/angle.h"
#include "astro/core/error.h"

namespace astro {

namespace {
// Schwarzschild radius of the Sun (au): 2 GM_sun / c^2.
constexpr double kSrs = 1.97412574336e-8;
constexpr double kSecondsPerJulianYear = constants::daysPerJulianYear * constants::secondsPerDay;

Aer horizonFromItrsDirection(const Vec3& dItrs, const Geodetic& site, double range) {
    const Aer a = enuToAer(ecefToEnuMatrix(site.lat, site.lon) * dItrs);
    return {a.az, a.el, range};
}

double applyRefraction(double el, const Atmosphere* atm) {
    return atm ? refractTrueToObserved(el, refractionConstants(*atm)) : el;
}
} // namespace

ObserverEphemeris observerEphemeris(const Time& t, const Ephemeris& eph, const StateVector& geo) {
    ObserverEphemeris o;
    StateVector earthHelio = eph.state(Body::Earth, Body::Sun, t);
    StateVector earthBary;
    try {
        earthBary = eph.state(Body::Earth, Body::SolarSystemBarycenter, t);
    } catch (const Error&) {
        earthBary = earthHelio;
    }
    o.barycentric = {earthBary.r + geo.r, earthBary.v + geo.v};
    o.heliocentric = earthHelio.r + geo.r;
    return o;
}

Vec3 lightDeflection(const Vec3& p, const Vec3& q, const Vec3& e, double emAu, double massSolar) noexcept {
    const double qdqpe = dot(q, q + e);
    const double w = massSolar * kSrs / emAu / std::max(qdqpe, 1e-9);
    return p + w * cross(p, cross(e, q));
}

Vec3 aberration(const Vec3& p, const Vec3& v, double s) noexcept {
    const double bm1 = std::sqrt(1.0 - norm2(v));
    const double pdv = dot(p, v);
    const double w1 = 1.0 + pdv / (1.0 + bm1);
    const double w2 = kSrs / s;
    return normalized(p * bm1 + w1 * v + w2 * (v - pdv * p));
}

Vec3 starBarycentricPosition(const CatalogStar& star, const Time& t) {
    const double dtYears = t.to(TimeScale::TDB).jd().daysSince(JulianDate{star.epochJdTdb, 0.0}) /
                           constants::daysPerJulianYear;
    const double sa = std::sin(star.ra), ca = std::cos(star.ra), sd = std::sin(star.dec), cd = std::cos(star.dec);
    const Vec3 u{cd * ca, cd * sa, sd};
    const Vec3 ea{-sa, ca, 0.0};             // direction of increasing RA
    const Vec3 ed{-sd * ca, -sd * sa, cd};   // direction of increasing Dec
    const Vec3 mu = star.pmRaCosDec * ea + star.pmDec * ed;  // rad / yr
    if (star.parallax <= 0.0) return normalized(u + mu * dtYears);
    const double dist = constants::au / star.parallax;  // m
    const Vec3 v = mu * dist + u * (star.radialVelocity * kSecondsPerJulianYear);  // m / yr
    return u * dist + v * dtYears;
}

Vec3 starApparentDirection(const CatalogStar& star, const Time& t, const ObserverEphemeris& obs) {
    const Vec3 sb = starBarycentricPosition(star, t);
    // Natural direction from the observer (parallax only when the distance is known).
    const Vec3 p = star.parallax > 0.0 ? normalized(sb - obs.barycentric.r) : sb;
    const double em = norm(obs.heliocentric);
    const Vec3 e = obs.heliocentric / em;
    const Vec3 pd = normalized(lightDeflection(p, p, e, em / constants::au));
    return aberration(pd, obs.barycentric.v / constants::c, em / constants::au);
}

RaDec starApparentGcrs(const CatalogStar& star, const Time& t, const Ephemeris& eph) {
    return raDecFromVector(starApparentDirection(star, t, observerEphemeris(t, eph)));
}

Aer starObserved(const CatalogStar& star, const Time& t, const Geodetic& site, const EopProvider& eop,
                 const Ephemeris& eph, const Atmosphere* atm) {
    const FrameTransformer ft(t, eop);
    const StateVector siteGcrs = ft.transform(StateVector{geodeticToEcef(site), {}}, Frame::ITRS, Frame::GCRS);
    const Vec3 dir = starApparentDirection(star, t, observerEphemeris(t, eph, siteGcrs));
    Aer a = horizonFromItrsDirection(ft.gcrsTo(Frame::ITRS) * dir, site, 0.0);
    a.el = applyRefraction(a.el, atm);
    return a;
}

Aer bodyObserved(Body body, const Time& t, const Geodetic& site, const EopProvider& eop, const Ephemeris& eph,
                 const Atmosphere* atm) {
    const FrameTransformer ft(t, eop);
    const StateVector siteGcrs = ft.transform(StateVector{geodeticToEcef(site), {}}, Frame::ITRS, Frame::GCRS);
    const ObserverEphemeris obs = observerEphemeris(t, eph, siteGcrs);
    const Vec3 vEarth = obs.barycentric.v - siteGcrs.v;
    // Light time: rho = body(t - tau) - observer(t), using geocentric body positions.
    double tau = 0.0;
    Vec3 rho;
    for (int i = 0; i < 4; ++i) {
        rho = eph.position(body, Body::Earth, t - tau) - vEarth * tau - siteGcrs.r;
        tau = norm(rho) / constants::c;
    }
    const Vec3 dir = aberration(normalized(rho), obs.barycentric.v / constants::c, norm(obs.heliocentric) / constants::au);
    Aer a = horizonFromItrsDirection(ft.gcrsTo(Frame::ITRS) * dir, site, norm(rho));
    a.el = applyRefraction(a.el, atm);
    return a;
}

Aer topocentricFromGcrs(const Vec3& rGcrs, const FrameTransformer& ft, const Geodetic& site, const Ellipsoid& ell) {
    const Vec3 rItrs = ft.gcrsTo(Frame::ITRS) * rGcrs;
    return ecefToAer(rItrs, site, ell);
}

double localMeanSiderealTime(const Time& t, const EopProvider& eop, double lon) {
    const FrameTransformer ft(t, eop);
    return wrapTwoPi(ft.gmst() + lon);
}

double localApparentSiderealTime(const Time& t, const EopProvider& eop, double lon) {
    const FrameTransformer ft(t, eop);
    return wrapTwoPi(ft.gast() + lon);
}

double hourAngle(double lst, double ra) noexcept {
    return wrapPi(lst - ra);
}

} // namespace astro
