/**********************************************************************
 * secular.h
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

#pragma once

#include <iosfwd>
#include <span>
#include <vector>

#include <boost/describe/enum.hpp>

#include "astro/core/constants.h"
#include "astro/orbit/elements.h"
#include "astro/orbit/state12.h"
#include "astro/orbit/state_covariance.h"

/// \file
/// Analytic secular (mean-element) drift due to the Earth's zonal harmonics,
/// after Brouwer (1959): first order in J2, plus J2^2 and J4.
namespace astro {

/// Central-body constants for the secular theory. Defaults: the Earth, with the radius and zonal
/// coefficients of EGM2008 and the gravitational parameter of constants::gmEarth, which is the
/// TCG-compatible IERS value and 7.5e-10 larger than EGM2008's own (see its description).
struct ZonalParameters {
    double mu{constants::gmEarth};            ///< gravitational parameter (m^3/s^2)
    double re{constants::earthRadiusEgm2008}; ///< reference radius of the harmonics (m)
    double J2{constants::J2};                 ///< unnormalized zonal coefficient J2
    double J4{constants::J4};                 ///< unnormalized zonal coefficient J4
};

/// Which zonal terms enter the secular rates.
enum class SecularOrder {
    J2,         ///< first-order J2 only
    J2SquaredJ4 ///< J2 + J2^2 + J4 (Brouwer)
};

BOOST_DESCRIBE_ENUM(SecularOrder, J2, J2SquaredJ4)

/// Writes "SecularOrder::<enumerator>", or "SecularOrder::INVALID_VALUE" for any other value.
std::ostream& operator<<(std::ostream& os, SecularOrder value);

/// Secular rates (rad/s) of the mean elements.
struct SecularRates {
    double raanDot{0.0};        ///< rate of the right ascension of the ascending node
    double argpDot{0.0};        ///< rate of the argument of periapsis
    double meanAnomalyDot{0.0}; ///< total rate including the mean motion
};

/// Rates for mean semi-major axis a, eccentricity e and inclination i.
[[nodiscard]] SecularRates secularRates(double                 a,
                                        double                 e,
                                        double                 i,
                                        const ZonalParameters& zp    = {},
                                        SecularOrder           order = SecularOrder::J2SquaredJ4) noexcept;

/// Propagate mean elements by dt seconds with the secular rates (a, e, i constant).
/// The true anomaly is advanced through the mean anomaly. Throws astro::Error unless the orbit is
/// elliptic (a > 0, 0 <= e < 1), as SecularPropagator does.
[[nodiscard]] KeplerianElements propagateSecular(const KeplerianElements& mean,
                                                 double                   dt,
                                                 const ZonalParameters&   zp    = {},
                                                 SecularOrder             order = SecularOrder::J2SquaredJ4);

/// propagateSecular for an inertial state whose z axis is the Earth's pole. The osculating
/// elements of s0 are treated as mean elements (no short-period correction), so expect errors at the
/// short-period level (kilometres in LEO); use NumericalPropagator::j2j4 when that matters.
/// Throws astro::Error unless s0 is on an elliptic orbit, as SecularPropagator does.
[[nodiscard]] StateVector propagateSecular(const StateVector&     s0,
                                           double                 dt,
                                           const ZonalParameters& zp    = {},
                                           SecularOrder           order = SecularOrder::J2SquaredJ4);

/// propagateSecular for a 12-component state (position, velocity, acceleration, jerk), with the
/// same mean-element caveat. The force model for acceleration, jerk and residuals is the point mass
/// plus J2 (SecularOrder::J2) or plus J2 and J4 (SecularOrder::J2SquaredJ4).
///  - ModelOnly: secular propagation from position and velocity.
///  - UnmodeledForce: the unmodeled force da0 + dj0 t is superposed on the secular trajectory,
///    dr = da0 dt^2/2 + dj0 dt^3/6, dv = da0 dt + dj0 dt^2/2, ignoring how the displacement feeds
///    back into gravity. That is fine while dr stays small compared with the orbit; otherwise use
///    NumericalPropagator::j2j4 with HigherOrderMode::UnmodeledForce.
///
/// Throws astro::Error unless s0 is on an elliptic orbit, as SecularPropagator does.
[[nodiscard]] Propagation12 propagateSecular(const StateVector12&   s0,
                                             double                 dt,
                                             const ZonalParameters& zp    = {},
                                             SecularOrder           order = SecularOrder::J2SquaredJ4,
                                             HigherOrderMode        mode  = HigherOrderMode::ModelOnly);

/// Secular J2/J4 propagation of one orbit to many times.
///
/// Everything that does not change with time (elements, secular rates, sin/cos of the inclination,
/// residuals) is set up once in the constructor. Each later state then costs a Kepler solve and a
/// few sines and cosines, so use this class instead of the free propagateSecular functions when
/// propagating the same orbit more than once. Results are the same as the free functions (which use
/// this class). Instances are immutable and safe to share between threads.
///
/// The z axis of the state's frame must be the Earth's pole (TOD, TEME or CIRS). For a GCRS or
/// J2000 state, rotate it with poleOfDateMatrix (frames/transform.h) first and rotate the results
/// back.
///
/// The velocity is the two-body velocity of the mean elements at each time (as in Vallado's
/// PKEPLER), so a state converts back to exactly the propagated mean elements and the input state
/// is reproduced at the epoch. It is not the time derivative of the position, which also carries
/// the node and perigee drift and the J2/J4 part of the mean-anomaly rate (1e-3 relative: 5-20 m/s
/// in LEO, depending on the inclination). SecularPropagatorCorrected returns that derivative instead; see the propagation guide
/// for when each fits.
class SecularPropagator {
public:
    /// From mean elements (elliptic: e < 1). Throws astro::Error unless a > 0 and 0 <= e < 1.
    explicit SecularPropagator(const KeplerianElements& mean,
                               const ZonalParameters&   zp    = {},
                               SecularOrder             order = SecularOrder::J2SquaredJ4);
    /// From an inertial state whose z axis is the Earth's pole; its osculating elements are taken
    /// as mean elements (see propagateSecular). Throws astro::Error unless s0 is on an elliptic
    /// orbit: for a hyperbolic, parabolic or rectilinear state (r = 0, or r parallel to v).
    explicit SecularPropagator(const StateVector&     s0,
                               const ZonalParameters& zp    = {},
                               SecularOrder           order = SecularOrder::J2SquaredJ4);
    /// From a 12-component state; `mode` is as for propagateSecular(const StateVector12&, ...).
    /// Throws astro::Error as the constructor from a StateVector does.
    SecularPropagator(const StateVector12&   s0,
                      const ZonalParameters& zp    = {},
                      SecularOrder           order = SecularOrder::J2SquaredJ4,
                      HigherOrderMode        mode  = HigherOrderMode::ModelOnly);

    /// Mean elements dt seconds after the epoch.
    [[nodiscard]] KeplerianElements elements(double dt) const;
    /// State dt seconds after the epoch.
    [[nodiscard]] StateVector propagate(double dt) const noexcept;
    /// 12-component state dt seconds after the epoch. Acceleration and jerk are the point mass plus
    /// zonal field (exact jerk), plus the unmodeled force in HigherOrderMode::UnmodeledForce.
    [[nodiscard]] StateVector12 propagate12(double dt) const noexcept;

    /// States at dt = start + k * step for k = 0 .. out.size() - 1, written into `out`.
    void propagateGrid(double start, double step, std::span<StateVector> out) const noexcept;
    /// As propagateGrid, for 12-component states.
    void propagateGrid12(double start, double step, std::span<StateVector12> out) const noexcept;
    /// As propagateGrid, returning `count` states.
    [[nodiscard]] std::vector<StateVector> propagateGrid(double start, double step, std::size_t count) const;
    /// As propagateGrid12, returning `count` states.
    [[nodiscard]] std::vector<StateVector12> propagateGrid12(double start, double step, std::size_t count) const;

    /// State dt seconds after the epoch and the state transition matrix Phi(dt) = d(state at dt) /
    /// d(epochState()), for the covariance of the input state.
    ///
    /// Phi is the derivative of the whole map from the input state to the propagated one (input
    /// state -> mean elements -> secular drift -> state), by fourth-order central differences of
    /// propagations from perturbed input states, with steps of 1e-5 |r0| and 1e-5 |v0|. Relative
    /// accuracy 2e-5 or better: against the closed-form two-body Phi (propagateKeplerStm) with
    /// J2 = J4 = 0, 1e-11 to 5e-9 for the LEO test orbits over up to a week, and 1.6e-5 for e = 0.74
    /// over a week (astro_study_propagation).
    /// For SecularPropagatorCorrected Phi(0) is not the identity: its velocity at the epoch is not
    /// the input velocity. Costs 24 extra propagations (and constructions); for many times use
    /// propagateGridWithCovariance, which builds the perturbed propagators once. The state is
    /// propagate(dt), or propagate12(dt) for HigherOrderMode::UnmodeledForce, whose input
    /// acceleration and jerk are held fixed. Throws astro::Error if a perturbed input state is not
    /// elliptic, which happens within ~1e-5 of e = 1.
    [[nodiscard]] StateTransition propagateStm(double dt) const;
    /// State and covariance dt seconds after the epoch: P = Phi P0 Phi^T, Phi from propagateStm. P0 is
    /// the covariance of epochState() (m^2, m^2/s, m^2/s^2) in its frame. There is no process noise:
    /// the covariance shows only how the model maps the initial uncertainty (see
    /// NumericalPropagator::propagateGridWithCovariance for process noise). Throws astro::Error as
    /// propagateStm does.
    [[nodiscard]] StateCovariance propagateWithCovariance(double dt, const Mat6& P0) const;
    /// States and covariances at dt = start + k * step for k = 0 .. out.size() - 1, written into
    /// `out`: the grid form of propagateWithCovariance, 25 grid propagations in all (1.2 us per
    /// output; NumericalPropagator::propagateGridWithCovariance is faster on long grids). Throws
    /// astro::Error as propagateStm does.
    void propagateGridWithCovariance(double start, double step, const Mat6& P0, std::span<StateCovariance> out) const;
    /// As above, returning `count` states and covariances. Throws astro::Error as propagateStm does.
    [[nodiscard]] std::vector<StateCovariance>
    propagateGridWithCovariance(double start, double step, std::size_t count, const Mat6& P0) const;
    /// The input state the covariance methods refer to: the state given to the constructor, or for
    /// mean elements elementsToState(mean, zp.mu).
    [[nodiscard]] StateVector epochState() const;

    /// Mean elements at the epoch.
    [[nodiscard]] const KeplerianElements& meanElements() const noexcept { return mean_; }
    /// Secular rates of the mean elements.
    [[nodiscard]] const SecularRates& rates() const noexcept { return rates_; }
    /// Input acceleration minus the model's at the epoch (zero unless built from a StateVector12).
    [[nodiscard]] const Vec3& accelerationResidual() const noexcept { return da0_; }
    /// Input jerk minus the model's at the epoch (zero unless built from a StateVector12).
    [[nodiscard]] const Vec3& jerkResidual() const noexcept { return dj0_; }
    /// True if velocities are the time derivative of the positions (SecularPropagatorCorrected),
    /// false if they are the two-body velocity of the mean elements (SecularPropagator).
    [[nodiscard]] bool velocityIsPositionRate() const noexcept { return rateVelocity_; }

protected:
    /// Selects the velocity of SecularPropagatorCorrected in the constructors below.
    struct RateVelocity {};
    /// As the public constructors, with velocities that are the time derivative of the positions.
    /// Throws astro::Error as they do.
    SecularPropagator(RateVelocity, const KeplerianElements& mean, const ZonalParameters& zp, SecularOrder order);
    /// As above, from an inertial state whose osculating elements are taken as mean elements.
    /// Throws astro::Error as the public constructor does.
    SecularPropagator(RateVelocity, const StateVector& s0, const ZonalParameters& zp, SecularOrder order);
    /// As above, from a 12-component state. Throws astro::Error as the public constructor does.
    SecularPropagator(RateVelocity,
                      const StateVector12&   s0,
                      const ZonalParameters& zp,
                      SecularOrder           order,
                      HigherOrderMode        mode);

private:
    /// Time-dependent angles of one state: reduced mean anomaly, eccentric anomaly and the sines
    /// and cosines needed to build the state.
    struct Angles {
        double M{0.0}, E{0.0}, sinE{0.0}, cosE{1.0};
        double sinRaan{0.0}, cosRaan{1.0}, sinArgp{0.0}, cosArgp{1.0};
    };

    void                        init();
    void                        initResiduals(const StateVector12& s0, HigherOrderMode mode);
    void                        fieldDerivatives(const Vec3& r, const Vec3& v, Vec3& a, Vec3& j) const noexcept;
    [[nodiscard]] Angles        angles(double dt) const noexcept;
    [[nodiscard]] StateVector   state(const Angles& x) const noexcept;
    [[nodiscard]] StateVector12 complete(double dt, StateVector s) const noexcept;
    template <class Emit>
    void grid(double start, double step, std::size_t count, Emit&& emit) const noexcept;
    /// The same kind of propagator from another input state (keeping acceleration, jerk and mode).
    /// Throws astro::Error unless `s` is on an elliptic orbit.
    [[nodiscard]] SecularPropagator rebuilt(const StateVector& s) const;
    /// States of the trajectory (including an unmodeled force) on a grid.
    void trajectoryGrid(double start, double step, std::span<StateVector> out) const;
    /// States and state transition matrices on a grid. Throws astro::Error as propagateStm does.
    void stmGrid(double start, double step, std::span<StateVector> states, std::span<Mat6> stm) const;

    KeplerianElements   mean_;
    ZonalParameters     zp_;
    SecularOrder        order_;
    SecularRates        rates_;
    std::vector<double> J_;       // zonal coefficients of the model, by degree
    double              M0_{0.0}; // mean anomaly at the epoch
    double              b_{0.0};  // semi-minor axis a sqrt(1 - e^2)
    double              sqrtMuA_{0.0};
    double              bOverA_{1.0};
    double              cosI_{1.0};
    double              sinI_{0.0};
    bool                unmodeled_{false};
    bool                rateVelocity_{false}; // velocity = d(position)/dt (SecularPropagatorCorrected)
    double              mDotOverN_{1.0};      // mean-anomaly rate over the two-body mean motion
    Vec3                da0_;
    Vec3                dj0_;
    StateVector12       input_;                            // the input state, when built from one
    HigherOrderMode     mode_{HigherOrderMode::ModelOnly}; // mode of a StateVector12 start
    bool                fromState_{false};                 // built from a state (not mean elements)
    bool                fromState12_{false};               // built from a StateVector12
};

/// SecularPropagator whose velocity is the time derivative of its position.
///
/// Positions are identical to SecularPropagator's. The velocity adds what the two-body velocity of
/// the mean elements leaves out, the drift of the angles:
///
///     v = (Mdot / n) v_mean + raanDot (z x r) + argpDot (w x r)
///
/// with n the two-body mean motion and w the orbit normal. Use it wherever velocity and position
/// must agree: range rate and Doppler, closest-approach searches (Trajectory uses it), Hermite
/// interpolation or differencing of outputs. The price: the velocity at the epoch is not the input
/// velocity (it differs by the same ~1e-3 relative), and a propagated state no longer converts back
/// to the mean elements. Acceleration and jerk in propagate12 are those of the field at the
/// propagated position and (corrected) velocity; the jerk residual of a 12-component start is taken
/// against the corrected velocity, so propagate12(0) reproduces the input jerk.
///
/// It is a SecularPropagator (pass it wherever one is accepted) and costs ~2 % more per output.
class SecularPropagatorCorrected : public SecularPropagator {
public:
    /// From mean elements (elliptic: e < 1). Throws astro::Error unless a > 0 and 0 <= e < 1.
    explicit SecularPropagatorCorrected(const KeplerianElements& mean,
                                        const ZonalParameters&   zp    = {},
                                        SecularOrder             order = SecularOrder::J2SquaredJ4);
    /// From an inertial state whose z axis is the Earth's pole; its osculating elements are taken
    /// as mean elements, as for SecularPropagator. Throws astro::Error as SecularPropagator does:
    /// unless s0 is on an elliptic orbit.
    explicit SecularPropagatorCorrected(const StateVector&     s0,
                                        const ZonalParameters& zp    = {},
                                        SecularOrder           order = SecularOrder::J2SquaredJ4);
    /// From a 12-component state; `mode` is as for SecularPropagator. Throws astro::Error as
    /// SecularPropagator does: unless s0 is on an elliptic orbit.
    SecularPropagatorCorrected(const StateVector12&   s0,
                               const ZonalParameters& zp    = {},
                               SecularOrder           order = SecularOrder::J2SquaredJ4,
                               HigherOrderMode        mode  = HigherOrderMode::ModelOnly);
};

/// Inclination (rad) of a sun-synchronous orbit with mean a and e (first-order J2;
/// nodal rate 360 deg per tropical year). Throws astro::Error unless the orbit is elliptic
/// (a > 0, 0 <= e < 1), or if no solution exists.
[[nodiscard]] double sunSynchronousInclination(double a, double e = 0.0, const ZonalParameters& zp = {});

/// Mean semi-major axis of a repeat-ground-track orbit that completes `revs`
/// revolutions in `days` nodal days, for given e and i (J2, iterative). Throws astro::Error if
/// revs or days is not positive, unless 0 <= e < 1, and if the iteration does not converge (to
/// 1e-13 of a): there is then no such orbit in the J2 model, e.g. 2 revolutions a day at e = 0.99.
[[nodiscard]] double
repeatGroundTrackSemiMajorAxis(int revs, int days, double e, double i, const ZonalParameters& zp = {});

/// Critical inclination (rad) where the J2 apsidal rotation vanishes: acos(1/sqrt(5)).
inline constexpr double kCriticalInclination = 1.1071487177940904;

} // namespace astro
