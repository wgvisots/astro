/**********************************************************************
 * secular.cpp
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

#include "astro/orbit/secular.h"

#include <cmath>

#include "../core/enum_io.h"
#include "astro/core/angle.h"
#include "astro/core/error.h"
#include "astro/orbit/kepler.h"
#include "astro/orbit/numerical.h"

namespace astro {

SecularRates secularRates(double a, double e, double i, const ZonalParameters& zp, SecularOrder order) noexcept {
    const double n = std::sqrt(zp.mu / (a * a * a));
    const double eta2 = 1.0 - e * e;
    const double eta = std::sqrt(eta2);
    const double th = std::cos(i);
    const double th2 = th * th, th4 = th2 * th2;
    const double ra = zp.re / a;
    // Brouwer's small parameters.
    const double g2 = 0.5 * zp.J2 * ra * ra;           // gamma2
    const double g2p = g2 / (eta2 * eta2);             // gamma2'
    const double g4 = -0.375 * zp.J4 * ra * ra * ra * ra;  // gamma4
    const double g4p = g4 / (eta2 * eta2 * eta2 * eta2);   // gamma4'

    SecularRates r;
    r.meanAnomalyDot = n * (1.0 + 1.5 * g2p * eta * (-1.0 + 3.0 * th2));
    r.argpDot = n * 1.5 * g2p * (-1.0 + 5.0 * th2);
    r.raanDot = n * -3.0 * g2p * th;
    if (order == SecularOrder::J2SquaredJ4) {
        r.meanAnomalyDot += n * (3.0 / 32.0 * g2p * g2p * eta *
                                     (-15.0 + 16.0 * eta + 25.0 * eta2 + (30.0 - 96.0 * eta - 90.0 * eta2) * th2 +
                                      (105.0 + 144.0 * eta + 25.0 * eta2) * th4) +
                                 15.0 / 16.0 * g4p * eta * e * e * (3.0 - 30.0 * th2 + 35.0 * th4));
        r.argpDot += n * (3.0 / 32.0 * g2p * g2p *
                              (-35.0 + 24.0 * eta + 25.0 * eta2 + (90.0 - 192.0 * eta - 126.0 * eta2) * th2 +
                               (385.0 + 360.0 * eta + 45.0 * eta2) * th4) +
                          5.0 / 16.0 * g4p *
                              (21.0 - 9.0 * eta2 + (-270.0 + 126.0 * eta2) * th2 + (385.0 - 189.0 * eta2) * th4));
        r.raanDot += n * (3.0 / 8.0 * g2p * g2p * ((-5.0 + 12.0 * eta + 9.0 * eta2) * th + (-35.0 - 36.0 * eta - 5.0 * eta2) * th2 * th) +
                          5.0 / 4.0 * g4p * (5.0 - 3.0 * eta2) * th * (3.0 - 7.0 * th2));
    }
    return r;
}

KeplerianElements propagateSecular(const KeplerianElements& mean, double dt, const ZonalParameters& zp,
                                   SecularOrder order) {
    if (mean.e >= 1.0) throw Error("propagateSecular: requires an elliptic orbit");
    const SecularRates r = secularRates(mean.a, mean.e, mean.i, zp, order);
    KeplerianElements out = mean;
    out.raan = wrapTwoPi(mean.raan + r.raanDot * dt);
    out.argp = wrapTwoPi(mean.argp + r.argpDot * dt);
    out.nu = wrapTwoPi(meanToTrueAnomaly(mean.meanAnomaly() + r.meanAnomalyDot * dt, mean.e));
    return out;
}

StateVector propagateSecular(const StateVector& s0, double dt, const ZonalParameters& zp, SecularOrder order) {
    return elementsToState(propagateSecular(stateToElements(s0, zp.mu), dt, zp, order), zp.mu);
}

Propagation12 propagateSecular(const StateVector12& s0, double dt, const ZonalParameters& zp, SecularOrder order,
                               HigherOrderMode mode) {
    const GravityField field(order == SecularOrder::J2 ? std::vector<double>{0.0, 0.0, zp.J2}
                                                       : std::vector<double>{0.0, 0.0, zp.J2, 0.0, zp.J4},
                             zp.mu, zp.re);
    const bool unmodeled = mode == HigherOrderMode::UnmodeledForce;
    Propagation12 out;
    const Vec3 a0 = field.acceleration(0.0, s0.r, s0.v);
    out.accelerationResidual = s0.a - a0;
    out.jerkResidual = s0.j - field.jerk(0.0, s0.r, s0.v, unmodeled ? s0.a : a0);

    StateVector s = propagateSecular(s0.rv(), dt, zp, order);
    Vec3 da, dj;
    if (unmodeled) {
        const Vec3& a = out.accelerationResidual;
        const Vec3& j = out.jerkResidual;
        s.r += a * (0.5 * dt * dt) + j * (dt * dt * dt / 6.0);
        s.v += a * dt + j * (0.5 * dt * dt);
        da = a + j * dt;
        dj = j;
    }
    const Vec3 a = field.acceleration(dt, s.r, s.v) + da;
    out.state = {s.r, s.v, a, field.jerk(dt, s.r, s.v, a) + dj};
    return out;
}

double sunSynchronousInclination(double a, double e, const ZonalParameters& zp) {
    constexpr double tropicalYear = 365.242190 * 86400.0;
    const double target = constants::twoPi / tropicalYear;
    const double n = std::sqrt(zp.mu / (a * a * a));
    const double p = a * (1.0 - e * e);
    const double c = -target / (1.5 * n * zp.J2 * (zp.re / p) * (zp.re / p));
    if (std::abs(c) > 1.0) throw Error("sunSynchronousInclination: no sun-synchronous solution for this orbit");
    return std::acos(c);
}

double repeatGroundTrackSemiMajorAxis(int revs, int days, double e, double i, const ZonalParameters& zp) {
    if (revs <= 0 || days <= 0) throw Error("repeatGroundTrackSemiMajorAxis: revs and days must be positive");
    // Nodal period condition: revs * (2 pi) / (Mdot + argpDot) = days * 2 pi / (omegaE - raanDot)
    const double we = constants::earthRotationRate;
    double a = std::cbrt(zp.mu / std::pow(we * revs / days, 2.0));
    for (int it = 0; it < 50; ++it) {
        const SecularRates r = secularRates(a, e, i, zp, SecularOrder::J2);
        const double nodalRate = (we - r.raanDot) * revs / days;  // required (Mdot + argpDot)
        const double nNeeded = nodalRate - (r.meanAnomalyDot - std::sqrt(zp.mu / (a * a * a))) - r.argpDot;
        const double aNew = std::cbrt(zp.mu / (nNeeded * nNeeded));
        if (std::abs(aNew - a) < 1e-6) return aNew;
        a = aNew;
    }
    return a;
}

std::ostream& operator<<(std::ostream& os, SecularOrder value) {
    return detail::writeEnum(os, "SecularOrder", value);
}

} // namespace astro
