/**********************************************************************
 * elements.cpp
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

#include "astro/orbit/elements.h"

#include <algorithm>
#include <cmath>

#include "astro/core/angle.h"
#include "astro/core/error.h"
#include "astro/core/rotation.h"
#include "astro/orbit/kepler.h"

namespace astro {

namespace {
constexpr double kSmall = 1e-11;
}

double KeplerianElements::meanAnomaly() const noexcept { return trueToMeanAnomaly(nu, e); }

double KeplerianElements::eccentricAnomaly() const noexcept {
    if (e < 1.0) return trueToEccentricAnomaly(nu, e);
    return trueToHyperbolicAnomaly(nu, e);
}

KeplerianElements stateToElements(const StateVector& s, double mu) {
    const Vec3& r = s.r;
    const Vec3& v = s.v;
    const double rm = norm(r);
    const Vec3 h = cross(r, v);
    const double hm = norm(h);
    if (rm == 0.0 || hm == 0.0) throw Error("stateToElements: degenerate (rectilinear) state");
    const Vec3 n{-h.y, h.x, 0.0};  // node vector k x h
    const double nm = norm(n);
    const Vec3 ev = ((norm2(v) - mu / rm) * r - dot(r, v) * v) / mu;
    KeplerianElements el;
    el.e = norm(ev);
    const double energy = 0.5 * norm2(v) - mu / rm;
    if (std::abs(el.e - 1.0) < 1e-14) throw Error("stateToElements: parabolic orbit");
    el.a = -mu / (2.0 * energy);
    el.i = std::acos(std::clamp(h.z / hm, -1.0, 1.0));

    const bool circular = el.e < kSmall;
    const bool equatorial = el.i < kSmall || std::abs(el.i - constants::pi) < kSmall;

    if (!equatorial) {
        el.raan = wrapTwoPi(std::atan2(n.y, n.x));
    }
    if (!circular && !equatorial) {
        el.argp = std::acos(std::clamp(dot(n, ev) / (nm * el.e), -1.0, 1.0));
        if (ev.z < 0.0) el.argp = constants::twoPi - el.argp;
        el.nu = std::acos(std::clamp(dot(ev, r) / (el.e * rm), -1.0, 1.0));
        if (dot(r, v) < 0.0) el.nu = constants::twoPi - el.nu;
    } else if (circular && !equatorial) {
        // Argument of latitude.
        el.nu = std::acos(std::clamp(dot(n, r) / (nm * rm), -1.0, 1.0));
        if (r.z < 0.0) el.nu = constants::twoPi - el.nu;
    } else if (!circular && equatorial) {
        // Longitude of periapsis.
        el.argp = wrapTwoPi(std::atan2(ev.y, ev.x));
        if (h.z < 0.0) el.argp = wrapTwoPi(-el.argp);
        el.nu = std::acos(std::clamp(dot(ev, r) / (el.e * rm), -1.0, 1.0));
        if (dot(r, v) < 0.0) el.nu = constants::twoPi - el.nu;
    } else {
        // True longitude.
        el.nu = wrapTwoPi(std::atan2(r.y, r.x));
        if (h.z < 0.0) el.nu = wrapTwoPi(-el.nu);
    }
    if (el.e >= 1.0) el.nu = wrapPi(el.nu);
    return el;
}

StateVector elementsToState(const KeplerianElements& el, double mu) {
    const double p = el.semiLatusRectum();
    if (!(p > 0.0)) throw Error("elementsToState: invalid semi-latus rectum (check a and e)");
    const double cn = std::cos(el.nu), sn = std::sin(el.nu);
    const double rPf = p / (1.0 + el.e * cn);
    const Vec3 rp{rPf * cn, rPf * sn, 0.0};
    const double k = std::sqrt(mu / p);
    const Vec3 vp{-k * sn, k * (el.e + cn), 0.0};
    // Perifocal -> inertial: R3(-raan) R1(-i) R3(-argp).
    const Mat3 Q = R3(-el.raan) * R1(-el.i) * R3(-el.argp);
    return {Q * rp, Q * vp};
}

KeplerianElements elementsFromMeanAnomaly(double a, double e, double i, double raan, double argp, double M) {
    return {a, e, i, raan, argp, meanToTrueAnomaly(M, e)};
}

EquinoctialElements keplerianToEquinoctial(const KeplerianElements& el) noexcept {
    const double lonPeri = el.raan + el.argp;
    const double t = std::tan(0.5 * el.i);
    return {el.a, el.e * std::sin(lonPeri), el.e * std::cos(lonPeri), t * std::sin(el.raan), t * std::cos(el.raan),
            wrapTwoPi(lonPeri + el.meanAnomaly())};
}

KeplerianElements equinoctialToKeplerian(const EquinoctialElements& q) {
    KeplerianElements el;
    el.a = q.a;
    el.e = std::hypot(q.h, q.k);
    el.i = 2.0 * std::atan(std::hypot(q.p, q.q));
    el.raan = (q.p == 0.0 && q.q == 0.0) ? 0.0 : wrapTwoPi(std::atan2(q.p, q.q));
    const double lonPeri = (q.h == 0.0 && q.k == 0.0) ? 0.0 : std::atan2(q.h, q.k);
    el.argp = wrapTwoPi(lonPeri - el.raan);
    const double M = wrapTwoPi(q.lambda - lonPeri);
    el.nu = wrapTwoPi(meanToTrueAnomaly(M, el.e));
    return el;
}

EquinoctialElements stateToEquinoctial(const StateVector& s, double mu) {
    // Vector formulation (Broucke & Cefola 1972), free of the classical singularities.
    const Vec3 w = normalized(cross(s.r, s.v));
    if (w.z <= -1.0 + 1e-14) throw Error("stateToEquinoctial: singular for retrograde equatorial orbits");
    EquinoctialElements q;
    q.p = w.x / (1.0 + w.z);   // tan(i/2) sin(raan)
    q.q = -w.y / (1.0 + w.z);  // tan(i/2) cos(raan)
    const double d = 1.0 + q.p * q.p + q.q * q.q;
    const Vec3 f = Vec3{1.0 - q.p * q.p + q.q * q.q, 2.0 * q.p * q.q, -2.0 * q.p} / d;
    const Vec3 g = Vec3{2.0 * q.p * q.q, 1.0 + q.p * q.p - q.q * q.q, 2.0 * q.q} / d;
    const double rm = norm(s.r);
    const Vec3 ev = ((norm2(s.v) - mu / rm) * s.r - dot(s.r, s.v) * s.v) / mu;
    q.a = 1.0 / (2.0 / rm - norm2(s.v) / mu);
    q.k = dot(ev, f);
    q.h = dot(ev, g);
    const double trueLon = std::atan2(dot(s.r, g), dot(s.r, f));
    const double e = std::hypot(q.h, q.k);
    const double lonPeri = (e == 0.0) ? 0.0 : std::atan2(q.h, q.k);
    q.lambda = wrapTwoPi(lonPeri + trueToMeanAnomaly(trueLon - lonPeri, e));
    return q;
}

StateVector equinoctialToState(const EquinoctialElements& q, double mu) {
    return elementsToState(equinoctialToKeplerian(q), mu);
}

double orbitalPeriod(double a, double mu) noexcept { return constants::twoPi * std::sqrt(a * a * a / mu); }
double meanMotion(double a, double mu) noexcept { return std::sqrt(mu / (a * a * a)); }
double semiMajorAxisFromMeanMotion(double n, double mu) noexcept { return std::cbrt(mu / (n * n)); }

double specificEnergy(const StateVector& s, double mu) noexcept {
    return 0.5 * norm2(s.v) - mu / norm(s.r);
}

double flightPathAngle(const StateVector& s) noexcept {
    return std::asin(std::clamp(dot(s.r, s.v) / (norm(s.r) * norm(s.v)), -1.0, 1.0));
}

} // namespace astro
