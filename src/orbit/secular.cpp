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
    // The same domain as SecularPropagator (a < 0 with e < 1 used to give NaN here).
    if (!(mean.e >= 0.0 && mean.e < 1.0) || !(mean.a > 0.0)) {
        throw Error("propagateSecular: requires an elliptic orbit (a > 0, 0 <= e < 1)");
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
///
/// Safeguarded: the root lies in [M - e, M + e] (f is increasing there), and the bracket shrinks
/// with every iterate. A step that leaves the bracket is replaced by bisection, so the method
/// converges from any start, including a poor grid guess near perigee of a very eccentric orbit.
inline void refineKepler(double M, double e, double& E, double& sinE, double& cosE) noexcept {
    double lo = M - e, hi = M + e;
    for (int it = 0; it < 80; ++it) {
        const double f = E - e * sinE - M;
        if (std::abs(f) <= 2e-15) {
            return;
        }
        if (f < 0.0) {
            lo = std::max(lo, E);
        } else {
            hi = std::min(hi, E);
        }
        if (hi - lo <= 4e-16 * (1.0 + std::abs(E))) {
            return; // bracket at rounding level
        }
        const double g  = 1.0 / (1.0 - e * cosE);
        const double d1 = -f * g;
        double       d  = d1 * (1.0 - 0.5 * d1 * e * sinE * g);
        if (!(E + d > lo && E + d < hi)) {
            d = 0.5 * (lo + hi) - E;
        }
        E += d;
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
    : SecularPropagator(stateToElements(s0, zp.mu), zp, order) {
    input_     = StateVector12(s0.r, s0.v, Vec3{}, Vec3{});
    fromState_ = true;
}

SecularPropagator::SecularPropagator(const StateVector12&   s0,
                                     const ZonalParameters& zp,
                                     SecularOrder           order,
                                     HigherOrderMode        mode)
    : SecularPropagator(s0.rv(), zp, order) {
    input_       = s0;
    mode_        = mode;
    fromState_   = true;
    fromState12_ = true;
    initResiduals(s0, mode);
}

SecularPropagator::SecularPropagator(RateVelocity,
                                     const KeplerianElements& mean,
                                     const ZonalParameters&   zp,
                                     SecularOrder             order)
    : mean_(mean)
    , zp_(zp)
    , order_(order)
    , rateVelocity_(true) {
    init();
}

SecularPropagator::SecularPropagator(RateVelocity           tag,
                                     const StateVector&     s0,
                                     const ZonalParameters& zp,
                                     SecularOrder           order)
    : SecularPropagator(tag, stateToElements(s0, zp.mu), zp, order) {
    input_     = StateVector12(s0.r, s0.v, Vec3{}, Vec3{});
    fromState_ = true;
}

SecularPropagator::SecularPropagator(RateVelocity           tag,
                                     const StateVector12&   s0,
                                     const ZonalParameters& zp,
                                     SecularOrder           order,
                                     HigherOrderMode        mode)
    : SecularPropagator(tag, s0.rv(), zp, order) {
    input_       = s0;
    mode_        = mode;
    fromState_   = true;
    fromState12_ = true;
    initResiduals(s0, mode);
}

void SecularPropagator::initResiduals(const StateVector12& s0, HigherOrderMode mode) {
    unmodeled_ = mode == HigherOrderMode::UnmodeledForce;
    // The model's velocity at the epoch: the input one, or its corrected form.
    const Vec3 v0 = rateVelocity_ ? propagate(0.0).v : s0.v;
    Vec3       a, j;
    fieldDerivatives(s0.r, v0, a, j); // point mass + zonal jerk does not depend on the acceleration
    da0_ = s0.a - a;
    dj0_ = s0.j - j;
}

SecularPropagatorCorrected::SecularPropagatorCorrected(const KeplerianElements& mean,
                                                       const ZonalParameters&   zp,
                                                       SecularOrder             order)
    : SecularPropagator(RateVelocity{}, mean, zp, order) {}

SecularPropagatorCorrected::SecularPropagatorCorrected(const StateVector&     s0,
                                                       const ZonalParameters& zp,
                                                       SecularOrder           order)
    : SecularPropagator(RateVelocity{}, s0, zp, order) {}

SecularPropagatorCorrected::SecularPropagatorCorrected(const StateVector12&   s0,
                                                       const ZonalParameters& zp,
                                                       SecularOrder           order,
                                                       HigherOrderMode        mode)
    : SecularPropagator(RateVelocity{}, s0, zp, order, mode) {}

void SecularPropagator::init() {
    if (!(mean_.e >= 0.0 && mean_.e < 1.0) || !(mean_.a > 0.0)) {
        throw Error("SecularPropagator: requires an elliptic orbit (a > 0, 0 <= e < 1)");
    }
    rates_     = secularRates(mean_.a, mean_.e, mean_.i, zp_, order_);
    J_         = order_ == SecularOrder::J2 ? std::vector<double>{0.0, 0.0, zp_.J2}
                                            : std::vector<double>{0.0, 0.0, zp_.J2, 0.0, zp_.J4};
    M0_        = mean_.meanAnomaly();
    b_         = mean_.a * std::sqrt(1.0 - mean_.e * mean_.e);
    sqrtMuA_   = std::sqrt(zp_.mu * mean_.a);
    bOverA_    = b_ / mean_.a;
    cosI_      = std::cos(mean_.i);
    sinI_      = std::sin(mean_.i);
    mDotOverN_ = rates_.meanAnomalyDot / std::sqrt(zp_.mu / (mean_.a * mean_.a * mean_.a));
}

void SecularPropagator::fieldDerivatives(const Vec3& r, const Vec3& v, Vec3& a, Vec3& j) const noexcept {
    detail::GravityPartials p;
    detail::gravity<true>(r, J_.data(), J_.size(), zp_.mu, zp_.re, true, nullptr, &p);
    a = p.acceleration;
    j = p.jerk(v);
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
    const Vec3   r = P * xp + Q * yp;
    const Vec3   v = P * vx + Q * vy;
    if (!rateVelocity_) {
        return {r, v};
    }
    // d r / dt: each angle's partial times its rate. dr/dM = v / n; the node turns the orbit about
    // z and the perigee about the orbit normal W = P x Q.
    const Vec3 W{sO * sinI_, -cO * sinI_, cosI_};
    const Vec3 zxr{-r.y, r.x, 0.0};
    return {r, mDotOverN_ * v + rates_.raanDot * zxr + rates_.argpDot * cross(W, r)};
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

namespace {

// Relative step of the state transition matrix's differences, |r0| and |v0| times this. The
// fourth-order stencil makes the truncation error ~(n dt * step)^4: 1e-5 keeps it below ~1e-5 over
// a week, and the rounding error of the differences (machine epsilon / step) near 1e-11. A step of
// 1e-6 gives 2e-10 everywhere in the study instead, and 1.6e-9 for e = 0.74 over a week.
constexpr double kStmStep = 1e-5;

// Offsets of the stencil, in steps, and their weights times 12 (f' = sum w f(x + o h) / (12 h)).
constexpr double kStencilOffset[4] = {1.0, -1.0, 2.0, -2.0};
constexpr double kStencilWeight[4] = {8.0, -8.0, -1.0, 1.0};

} // namespace

StateVector SecularPropagator::epochState() const {
    // Mean elements are converted on demand, so that constructions which never ask (every
    // propagateSecular call) do not pay for it.
    return fromState_ ? input_.rv() : elementsToState(mean_, zp_.mu);
}

SecularPropagator SecularPropagator::rebuilt(const StateVector& s) const {
    StateVector12 x = input_;
    x.r             = s.r;
    x.v             = s.v;
    if (fromState12_) {
        return rateVelocity_ ? SecularPropagator(RateVelocity{}, x, zp_, order_, mode_)
                             : SecularPropagator(x, zp_, order_, mode_);
    }
    return rateVelocity_ ? SecularPropagator(RateVelocity{}, x.rv(), zp_, order_)
                         : SecularPropagator(x.rv(), zp_, order_);
}

void SecularPropagator::trajectoryGrid(double start, double step, std::span<StateVector> out) const {
    if (!unmodeled_) {
        propagateGrid(start, step, out);
        return;
    }
    grid(start, step, out.size(), [&](std::size_t k, double dt, const StateVector& s) {
        out[k] = complete(dt, s).rv();
    });
}

void SecularPropagator::stmGrid(double start, double step, std::span<StateVector> states, std::span<Mat6> stm) const {
    trajectoryGrid(start, step, states);
    for (Mat6& m : stm) {
        m = Mat6{};
    }
    std::vector<StateVector> buf(states.size());
    const StateVector        epoch = epochState();
    for (std::size_t j = 0; j < 6; ++j) {
        const double h = kStmStep * (j < 3 ? norm(epoch.r) : norm(epoch.v));
        for (std::size_t p = 0; p < 4; ++p) {
            StateVector s0 = epoch;
            // Component j of the state. (Not by indexing from the address of r.x: that reaches
            // outside the object pointed to, and an optimizer may then drop the write.)
            (j < 3 ? s0.r : s0.v)[j % 3] += kStencilOffset[p] * h;
            rebuilt(s0).trajectoryGrid(start, step, buf);
            const double w = kStencilWeight[p] / (12.0 * h);
            for (std::size_t k = 0; k < buf.size(); ++k) {
                const double c[6] = {buf[k].r.x, buf[k].r.y, buf[k].r.z, buf[k].v.x, buf[k].v.y, buf[k].v.z};
                for (std::size_t i = 0; i < 6; ++i) {
                    stm[k][i][j] += w * c[i];
                }
            }
        }
    }
}

StateTransition SecularPropagator::propagateStm(double dt) const {
    StateTransition out;
    stmGrid(dt, 1.0, std::span<StateVector>(&out.state, 1), std::span<Mat6>(&out.stm, 1));
    return out;
}

StateCovariance SecularPropagator::propagateWithCovariance(double dt, const Mat6& P0) const {
    const StateTransition t = propagateStm(dt);
    return {t.state, congruence(t.stm, P0)};
}

void SecularPropagator::propagateGridWithCovariance(double                     start,
                                                    double                     step,
                                                    const Mat6&                P0,
                                                    std::span<StateCovariance> out) const {
    std::vector<StateVector> states(out.size());
    std::vector<Mat6>        stm(out.size());
    stmGrid(start, step, states, stm);
    for (std::size_t k = 0; k < out.size(); ++k) {
        out[k] = {states[k], congruence(stm[k], P0)};
    }
}

std::vector<StateCovariance>
SecularPropagator::propagateGridWithCovariance(double start, double step, std::size_t count, const Mat6& P0) const {
    std::vector<StateCovariance> out(count);
    propagateGridWithCovariance(start, step, P0, out);
    return out;
}

double sunSynchronousInclination(double a, double e, const ZonalParameters& zp) {
    // a = 0 used to give pi/2 (the formula below then divides by infinity), and e > 1 an angle.
    if (!(e >= 0.0 && e < 1.0) || !(a > 0.0)) {
        throw Error("sunSynchronousInclination: requires an elliptic orbit (a > 0, 0 <= e < 1)");
    }
    constexpr double tropicalYear = 365.242190 * 86400.0;
    const double     target       = constants::twoPi / tropicalYear;
    const double     n            = std::sqrt(zp.mu / (a * a * a));
    const double     p            = a * (1.0 - e * e);
    const double     c            = -target / (1.5 * n * zp.J2 * (zp.re / p) * (zp.re / p));
    if (!(std::abs(c) <= 1.0)) {
        throw Error("sunSynchronousInclination: no sun-synchronous solution for this orbit");
    }
    return std::acos(c);
}

double repeatGroundTrackSemiMajorAxis(int revs, int days, double e, double i, const ZonalParameters& zp) {
    if (revs <= 0 || days <= 0) {
        throw Error("repeatGroundTrackSemiMajorAxis: revs and days must be positive");
    }
    if (!(e >= 0.0 && e < 1.0)) {
        throw Error("repeatGroundTrackSemiMajorAxis: requires an elliptic orbit (0 <= e < 1)");
    }
    // Nodal period condition: revs * (2 pi) / (Mdot + argpDot) = days * 2 pi / (omegaE - raanDot)
    const double we = constants::earthRotationRate;
    double       a  = std::cbrt(zp.mu / std::pow(we * revs / days, 2.0));
    for (int it = 0; it < 50; ++it) {
        const SecularRates r         = secularRates(a, e, i, zp, SecularOrder::J2);
        const double       nodalRate = (we - r.raanDot) * revs / days; // required (Mdot + argpDot)
        const double       nNeeded   = nodalRate - (r.meanAnomalyDot - std::sqrt(zp.mu / (a * a * a))) - r.argpDot;
        if (!(nNeeded > 0.0)) {
            break; // a mean motion that is zero, negative or NaN: the iteration has left its domain
        }
        const double aNew = std::cbrt(zp.mu / (nNeeded * nNeeded));
        // Relative to a (1e-6 m at 10 000 km). An absolute tolerance is also met when the iteration
        // diverges towards a = 0, which it does when there is no solution: 2 revolutions a day at
        // e = 0.99 "converged" to 0.000 km.
        if (std::abs(aNew - a) < 1e-13 * aNew) {
            return aNew;
        }
        a = aNew;
    }
    throw Error("repeatGroundTrackSemiMajorAxis: did not converge (no such orbit in the J2 model)");
}

std::ostream& operator<<(std::ostream& os, SecularOrder value) {
    return detail::writeEnum(os, "SecularOrder", value);
}

} // namespace astro
