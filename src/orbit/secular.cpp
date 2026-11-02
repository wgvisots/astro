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
#include "astro/orbit/gravity.h"
#include "astro/orbit/kepler.h"
#include "gravity_detail.h"

namespace astro {

SecularRates secularRates(double a, double e, double i, const ZonalParameters& zp, SecularOrder order) noexcept {
    const double n    = std::sqrt(zp.mu / (a * a * a));
    const double eta2 = 1.0 - e * e;
    const double eta  = std::sqrt(eta2);
    const double th   = std::cos(i);
    const double th2 = th * th, th4 = th2 * th2;
    const double ra = zp.re / a;
    // Brouwer's small parameters.
    const double g2  = 0.5 * zp.J2 * ra * ra;              // gamma2
    const double g2p = g2 / (eta2 * eta2);                 // gamma2'
    const double g4  = -0.375 * zp.J4 * ra * ra * ra * ra; // gamma4
    const double g4p = g4 / (eta2 * eta2 * eta2 * eta2);   // gamma4'

    SecularRates r;
    r.meanAnomalyDot = n * (1.0 + 1.5 * g2p * eta * (-1.0 + 3.0 * th2));
    r.argpDot        = n * 1.5 * g2p * (-1.0 + 5.0 * th2);
    r.raanDot        = n * -3.0 * g2p * th;
    if (order == SecularOrder::J2SquaredJ4) {
        r.meanAnomalyDot += n
                          * (3.0 / 32.0 * g2p * g2p * eta
                                 * (-15.0 + 16.0 * eta + 25.0 * eta2 + (30.0 - 96.0 * eta - 90.0 * eta2) * th2
                                    + (105.0 + 144.0 * eta + 25.0 * eta2) * th4)
                             + 15.0 / 16.0 * g4p * eta * e * e * (3.0 - 30.0 * th2 + 35.0 * th4));
        r.argpDot +=
            n
            * (3.0 / 32.0 * g2p * g2p
                   * (-35.0 + 24.0 * eta + 25.0 * eta2 + (90.0 - 192.0 * eta - 126.0 * eta2) * th2
                      + (385.0 + 360.0 * eta + 45.0 * eta2) * th4)
               + 5.0 / 16.0 * g4p * (21.0 - 9.0 * eta2 + (-270.0 + 126.0 * eta2) * th2 + (385.0 - 189.0 * eta2) * th4));
        r.raanDot += n
                   * (3.0 / 8.0 * g2p * g2p
                          * ((-5.0 + 12.0 * eta + 9.0 * eta2) * th + (-35.0 - 36.0 * eta - 5.0 * eta2) * th2 * th)
                      + 5.0 / 4.0 * g4p * (5.0 - 3.0 * eta2) * th * (3.0 - 7.0 * th2));
    }
    return r;
}

KeplerianElements
propagateSecular(const KeplerianElements& mean, double dt, const ZonalParameters& zp, SecularOrder order) {
    if (mean.e >= 1.0) {
        throw Error("propagateSecular: requires an elliptic orbit");
    }
    const SecularRates r   = secularRates(mean.a, mean.e, mean.i, zp, order);
    KeplerianElements  out = mean;
    out.raan               = wrapTwoPi(mean.raan + r.raanDot * dt);
    out.argp               = wrapTwoPi(mean.argp + r.argpDot * dt);
    out.nu                 = wrapTwoPi(meanToTrueAnomaly(mean.meanAnomaly() + r.meanAnomalyDot * dt, mean.e));
    return out;
}

StateVector propagateSecular(const StateVector& s0, double dt, const ZonalParameters& zp, SecularOrder order) {
    return SecularPropagator(s0, zp, order).propagate(dt);
}

Propagation12 propagateSecular(const StateVector12&   s0,
                               double                 dt,
                               const ZonalParameters& zp,
                               SecularOrder           order,
                               HigherOrderMode        mode) {
    const SecularPropagator prop(s0, zp, order, mode);
    return {prop.propagate12(dt), prop.accelerationResidual(), prop.jerkResidual()};
}

// ------------------------------------------------------------ SecularPropagator

namespace {

/// sin and cos of x + d from those of x, for |d| <= 0.1, using Taylor series for sin d and cos d
/// (truncation error below 1e-19).
inline void rotateSmall(double& s, double& c, double d) noexcept {
    const double d2 = d * d;
    const double sd =
        d * (1.0 - d2 / 6.0 * (1.0 - d2 / 20.0 * (1.0 - d2 / 42.0 * (1.0 - d2 / 72.0 * (1.0 - d2 / 110.0)))));
    const double cd = 1.0 - d2 / 2.0 * (1.0 - d2 / 12.0 * (1.0 - d2 / 30.0 * (1.0 - d2 / 56.0 * (1.0 - d2 / 90.0))));
    const double ns = s * cd + c * sd;
    c               = c * cd - s * sd;
    s               = ns;
}

/// sin and cos of x + d from those of x and of d.
inline void rotate(double& s, double& c, double sd, double cd) noexcept {
    const double ns = s * cd + c * sd;
    c               = c * cd - s * sd;
    s               = ns;
}

/// Refine E (with its sin and cos) to solve Kepler's equation E - e sin E = M, for M in [-pi, pi].
///
/// Chebyshev's method (third order, one division per correction): d = d1 - f'' d1^2 / (2 f') with
/// d1 = -f / f'. Small corrections update sin E and cos E with rotateSmall instead of new library
/// calls, so a good starting E costs only a few multiplies.
inline void refineKepler(double M, double e, double& E, double& sinE, double& cosE) noexcept {
    for (int it = 0; it < 40; ++it) {
        const double f = E - e * sinE - M;
        if (std::abs(f) <= 2e-15) {
            return;
        }
        const double g   = 1.0 / (1.0 - e * cosE);
        const double d1  = -f * g;
        const double d   = d1 * (1.0 - 0.5 * d1 * e * sinE * g);
        E               += d;
        if (std::abs(d) <= 0.1) {
            rotateSmall(sinE, cosE, d);
        } else {
            sinE = std::sin(E);
            cosE = std::cos(E);
        }
    }
}

/// Solve Kepler's equation from scratch for M in [-pi, pi].
inline void solveKepler(double M, double e, double& E, double& sinE, double& cosE) noexcept {
    E    = e < 0.8 ? M : M + 0.85 * e * (M >= 0.0 ? 1.0 : -1.0);
    sinE = std::sin(E);
    cosE = std::cos(E);
    refineKepler(M, e, E, sinE, cosE);
}

inline double reduceAngle(double x) noexcept {
    constexpr double inv2Pi = 1.0 / constants::twoPi;
    return x - constants::twoPi * std::nearbyint(x * inv2Pi);
}

} // namespace

SecularPropagator::SecularPropagator(const KeplerianElements& mean, const ZonalParameters& zp, SecularOrder order)
    : mean_(mean)
    , zp_(zp)
    , order_(order) {
    init();
}

SecularPropagator::SecularPropagator(const StateVector& s0, const ZonalParameters& zp, SecularOrder order)
    : SecularPropagator(stateToElements(s0, zp.mu), zp, order) {}

SecularPropagator::SecularPropagator(const StateVector12&   s0,
                                     const ZonalParameters& zp,
                                     SecularOrder           order,
                                     HigherOrderMode        mode)
    : SecularPropagator(s0.rv(), zp, order) {
    unmodeled_ = mode == HigherOrderMode::UnmodeledForce;
    Vec3 a, j;
    fieldDerivatives(s0.r, s0.v, a, j); // point mass + zonal jerk does not depend on the acceleration
    da0_ = s0.a - a;
    dj0_ = s0.j - j;
}

void SecularPropagator::init() {
    if (!(mean_.e >= 0.0 && mean_.e < 1.0) || !(mean_.a > 0.0)) {
        throw Error("SecularPropagator: requires an elliptic orbit (a > 0, 0 <= e < 1)");
    }
    rates_   = secularRates(mean_.a, mean_.e, mean_.i, zp_, order_);
    J_       = order_ == SecularOrder::J2 ? std::vector<double>{0.0, 0.0, zp_.J2}
                                          : std::vector<double>{0.0, 0.0, zp_.J2, 0.0, zp_.J4};
    M0_      = mean_.meanAnomaly();
    b_       = mean_.a * std::sqrt(1.0 - mean_.e * mean_.e);
    sqrtMuA_ = std::sqrt(zp_.mu * mean_.a);
    bOverA_  = b_ / mean_.a;
    cosI_    = std::cos(mean_.i);
    sinI_    = std::sin(mean_.i);
}

void SecularPropagator::fieldDerivatives(const Vec3& r, const Vec3& v, Vec3& a, Vec3& j) const noexcept {
    detail::gravity<true>(r, v, J_.data(), J_.size(), zp_.mu, zp_.re, true, &a, &j);
}

KeplerianElements SecularPropagator::elements(double dt) const {
    return propagateSecular(mean_, dt, zp_, order_);
}

StateVector SecularPropagator::state(const Angles& x) const noexcept {
    // Perifocal position and velocity.
    const double xp = mean_.a * (x.cosE - mean_.e);
    const double yp = b_ * x.sinE;
    const double k  = sqrtMuA_ / (mean_.a * (1.0 - mean_.e * x.cosE)); // sqrt(mu a) / r
    const double vx = -k * x.sinE;
    const double vy = k * bOverA_ * x.cosE;

    // Perifocal -> inertial: P and Q unit vectors from the drifting node and perigee.
    const double cO = x.cosRaan, sO = x.sinRaan, cw = x.cosArgp, sw = x.sinArgp;
    const Vec3   P{cO * cw - sO * sw * cosI_, sO * cw + cO * sw * cosI_, sw * sinI_};
    const Vec3   Q{-cO * sw - sO * cw * cosI_, -sO * sw + cO * cw * cosI_, cw * sinI_};
    return {P * xp + Q * yp, P * vx + Q * vy};
}

SecularPropagator::Angles SecularPropagator::angles(double dt) const noexcept {
    Angles x;
    x.M = reduceAngle(M0_ + rates_.meanAnomalyDot * dt);
    solveKepler(x.M, mean_.e, x.E, x.sinE, x.cosE);
    const double raan = mean_.raan + rates_.raanDot * dt;
    const double argp = mean_.argp + rates_.argpDot * dt;
    x.cosRaan         = std::cos(raan);
    x.sinRaan         = std::sin(raan);
    x.cosArgp         = std::cos(argp);
    x.sinArgp         = std::sin(argp);
    return x;
}

StateVector SecularPropagator::propagate(double dt) const noexcept {
    return state(angles(dt));
}

StateVector12 SecularPropagator::complete(double dt, StateVector s) const noexcept {
    Vec3 da, dj;
    if (unmodeled_) {
        // Superpose the unmodeled force on the secular trajectory (see propagateSecular).
        s.r += da0_ * (0.5 * dt * dt) + dj0_ * (dt * dt * dt / 6.0);
        s.v += da0_ * dt + dj0_ * (0.5 * dt * dt);
        da   = da0_ + dj0_ * dt;
        dj   = dj0_;
    }
    Vec3 a, j;
    fieldDerivatives(s.r, s.v, a, j);
    return {s.r, s.v, a + da, j + dj};
}

StateVector12 SecularPropagator::propagate12(double dt) const noexcept {
    return complete(dt, propagate(dt));
}

template <class Emit>
void SecularPropagator::grid(double start, double step, std::size_t count, Emit&& emit) const noexcept {
    // The angles advance by constant increments, so their sines and cosines are stepped with
    // angle-addition formulas instead of library calls. Kepler's equation starts from the previous
    // E plus the mean-anomaly increment. Every kResync outputs the angles are recomputed from
    // scratch, which keeps the accumulated rounding below ~1e-14 rad.
    constexpr std::size_t kResync = 64;
    const double          dM      = rates_.meanAnomalyDot * step;
    const double          sdM = std::sin(dM), cdM = std::cos(dM);
    const double          sdO = std::sin(rates_.raanDot * step), cdO = std::cos(rates_.raanDot * step);
    const double          sdw = std::sin(rates_.argpDot * step), cdw = std::cos(rates_.argpDot * step);
    const bool            incremental = std::abs(dM) < 1.0; // otherwise every output is computed directly
    Angles                x;
    for (std::size_t k = 0; k < count; ++k) {
        const double dt = start + step * static_cast<double>(k);
        if (!incremental || k % kResync == 0) {
            x = angles(dt);
        } else {
            // Advance M and E together, wrapping both at +-pi so that E - e sin E - M stays small.
            x.M += dM;
            x.E += dM;
            if (x.M > constants::pi) {
                x.M -= constants::twoPi;
                x.E -= constants::twoPi;
            } else if (x.M < -constants::pi) {
                x.M += constants::twoPi;
                x.E += constants::twoPi;
            }
            rotate(x.sinE, x.cosE, sdM, cdM);
            refineKepler(x.M, mean_.e, x.E, x.sinE, x.cosE);
            rotate(x.sinRaan, x.cosRaan, sdO, cdO);
            rotate(x.sinArgp, x.cosArgp, sdw, cdw);
        }
        emit(k, dt, state(x));
    }
}

void SecularPropagator::propagateGrid(double start, double step, std::span<StateVector> out) const noexcept {
    grid(start, step, out.size(), [&](std::size_t k, double, const StateVector& s) { out[k] = s; });
}

void SecularPropagator::propagateGrid12(double start, double step, std::span<StateVector12> out) const noexcept {
    grid(start, step, out.size(), [&](std::size_t k, double dt, const StateVector& s) { out[k] = complete(dt, s); });
}

std::vector<StateVector> SecularPropagator::propagateGrid(double start, double step, std::size_t count) const {
    std::vector<StateVector> out(count);
    propagateGrid(start, step, out);
    return out;
}

std::vector<StateVector12> SecularPropagator::propagateGrid12(double start, double step, std::size_t count) const {
    std::vector<StateVector12> out(count);
    propagateGrid12(start, step, out);
    return out;
}

double sunSynchronousInclination(double a, double e, const ZonalParameters& zp) {
    constexpr double tropicalYear = 365.242190 * 86400.0;
    const double     target       = constants::twoPi / tropicalYear;
    const double     n            = std::sqrt(zp.mu / (a * a * a));
    const double     p            = a * (1.0 - e * e);
    const double     c            = -target / (1.5 * n * zp.J2 * (zp.re / p) * (zp.re / p));
    if (std::abs(c) > 1.0) {
        throw Error("sunSynchronousInclination: no sun-synchronous solution for this orbit");
    }
    return std::acos(c);
}

double repeatGroundTrackSemiMajorAxis(int revs, int days, double e, double i, const ZonalParameters& zp) {
    if (revs <= 0 || days <= 0) {
        throw Error("repeatGroundTrackSemiMajorAxis: revs and days must be positive");
    }
    // Nodal period condition: revs * (2 pi) / (Mdot + argpDot) = days * 2 pi / (omegaE - raanDot)
    const double we = constants::earthRotationRate;
    double       a  = std::cbrt(zp.mu / std::pow(we * revs / days, 2.0));
    for (int it = 0; it < 50; ++it) {
        const SecularRates r         = secularRates(a, e, i, zp, SecularOrder::J2);
        const double       nodalRate = (we - r.raanDot) * revs / days; // required (Mdot + argpDot)
        const double       nNeeded   = nodalRate - (r.meanAnomalyDot - std::sqrt(zp.mu / (a * a * a))) - r.argpDot;
        const double       aNew      = std::cbrt(zp.mu / (nNeeded * nNeeded));
        if (std::abs(aNew - a) < 1e-6) {
            return aNew;
        }
        a = aNew;
    }
    return a;
}

std::ostream& operator<<(std::ostream& os, SecularOrder value) {
    return detail::writeEnum(os, "SecularOrder", value);
}

} // namespace astro
