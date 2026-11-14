/**********************************************************************
 * numerical.h
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

#include <functional>
#include <iosfwd>
#include <memory>
#include <optional>
#include <span>
#include <vector>

#include <boost/describe/enum.hpp>

#include "astro/core/constants.h"
#include "astro/core/mat6.h"
#include "astro/core/vec.h"
#include "astro/ephemeris/ephemeris.h"
#include "astro/frames/transform.h"
#include "astro/orbit/integrators.h"
#include "astro/orbit/state12.h"
#include "astro/orbit/state_covariance.h"
#include "astro/time/time.h"

/// \file
/// Numerical orbit propagation with pluggable force models.
namespace astro {

/// A perturbing (or central) acceleration. `t` is seconds since the propagator epoch;
/// r, v are in the propagator's inertial frame.
///
/// A force that depends on `t` alone, or has a part that does (a thrust profile, for example), is
/// invisible to the integrator's step-size control: see integrateRkf78Steps. Bound the step with
/// AdaptiveOptions::maxStep if it varies smoothly, and start a new propagation where it jumps.
class ForceModel {
public:
    virtual ~ForceModel() = default;
    /// Acceleration (m/s^2) at `t` for position r (m) and velocity v (m/s). Must be safe to call
    /// concurrently if the propagator using it is shared between threads.
    [[nodiscard]] virtual Vec3 acceleration(double t, const Vec3& r, const Vec3& v) const = 0;

    /// Time derivative of acceleration() along the trajectory, d/dt a(t, r(t), v(t)), given the
    /// total acceleration acting on the body. The default is a central difference along
    /// (1, v, totalAcceleration) with a step of 1e-5 |r| / |v| seconds (about 0.01 s in LEO), which
    /// is accurate to roughly 1e-9 relative for orbital forces. Override it for an exact form.
    [[nodiscard]] virtual Vec3 jerk(double t, const Vec3& r, const Vec3& v, const Vec3& totalAcceleration) const;

    /// acceleration(t, r, v), and its partial derivatives with respect to position (dadr) and
    /// velocity (dadv), for the state transition matrix. The default uses central differences
    /// (12 extra acceleration evaluations); override it with an exact form where one exists.
    virtual Vec3 partials(double t, const Vec3& r, const Vec3& v, Mat3& dadr, Mat3& dadv) const;
};

/// Central body gravity: point mass plus zonal harmonics J2..Jn.
///
/// `toBodyFrame` rotates the propagation frame into one whose z axis is the
/// rotation axis (e.g. the GCRS -> CIRS matrix at epoch); identity treats the
/// propagation frame's z axis as the pole.
class GravityField final : public ForceModel {
public:
    /// J indexed by degree (J[2] = J2, ...). Empty -> point mass only.
    explicit GravityField(std::vector<double> J           = {0.0, 0.0, constants::J2, constants::J3, constants::J4},
                          double              mu          = constants::gmEarth,
                          double              re          = constants::earthRadiusEgm2008,
                          const Mat3&         toBodyFrame = Mat3::identity());
    /// Point-mass plus zonal acceleration (m/s^2) in the propagation frame; t and v are unused.
    [[nodiscard]] Vec3 acceleration(double t, const Vec3& r, const Vec3& v) const override;
    /// Exact (point mass and zonal terms).
    [[nodiscard]] Vec3 jerk(double t, const Vec3& r, const Vec3& v, const Vec3& totalAcceleration) const override;
    /// Exact gravity gradient; no velocity dependence.
    Vec3 partials(double t, const Vec3& r, const Vec3& v, Mat3& dadr, Mat3& dadv) const override;
    /// Gravitational parameter of the field (m^3/s^2).
    [[nodiscard]] double mu() const noexcept { return mu_; }

private:
    std::vector<double> J_;
    double              mu_;
    double              re_;
    Mat3                R_;
    bool                identity_;
};

/// Point-mass perturbation from a third body (e.g. Sun, Moon) using an ephemeris.
/// The propagation frame is assumed to be the GCRS.
class ThirdBodyForce final : public ForceModel {
public:
    /// Perturbation by `body` with gravitational parameter `muBody` (m^3/s^2); its geocentric
    /// position comes from `ephemeris` at `epoch` + t. Throws astro::Error if `ephemeris` is null.
    ThirdBodyForce(Body body, double muBody, std::shared_ptr<const Ephemeris> ephemeris, const Time& epoch);
    /// Tidal acceleration (m/s^2) on a satellite at geocentric r: the body's pull on the satellite
    /// minus its pull on the Earth. v is unused.
    [[nodiscard]] Vec3 acceleration(double t, const Vec3& r, const Vec3& v) const override;
    /// Exact tidal gradient; no velocity dependence.
    Vec3 partials(double t, const Vec3& r, const Vec3& v, Mat3& dadr, Mat3& dadv) const override;

private:
    Body                             body_;
    double                           mu_;
    std::shared_ptr<const Ephemeris> eph_;
    Time                             epoch_;
};

/// An unmodeled inertial acceleration that varies linearly in time: a(t) = a0 + j0 t, with t in
/// seconds since the propagator epoch. HigherOrderMode::UnmodeledForce adds one of these.
class UnmodeledAcceleration final : public ForceModel {
public:
    /// `a0` is the acceleration at the epoch (m/s^2), `j0` its constant rate (m/s^3).
    UnmodeledAcceleration(const Vec3& a0, const Vec3& j0) noexcept
        : a0_(a0)
        , j0_(j0) {}
    /// a0 + j0 t (m/s^2); r and v are unused.
    [[nodiscard]] Vec3 acceleration(double t, const Vec3& r, const Vec3& v) const override;
    /// j0 (m/s^3).
    [[nodiscard]] Vec3 jerk(double t, const Vec3& r, const Vec3& v, const Vec3& totalAcceleration) const override;
    /// Zero: the unmodeled acceleration depends on time only.
    Vec3 partials(double t, const Vec3& r, const Vec3& v, Mat3& dadr, Mat3& dadv) const override;

private:
    Vec3 a0_;
    Vec3 j0_;
};

/// Which zero crossings of an event function count as events (see propagateToEvent).
enum class Crossing {
    Any,     ///< either direction
    Rising,  ///< the function goes from negative to positive as time increases
    Falling, ///< the function goes from positive to negative as time increases (e.g. height at impact)
};

BOOST_DESCRIBE_ENUM(Crossing, Any, Rising, Falling)

/// Writes "Crossing::<enumerator>", or "Crossing::INVALID_VALUE" for any other value.
std::ostream& operator<<(std::ostream& os, Crossing value);

/// Event function for propagateToEvent: a value at dt seconds after the epoch for the state there.
/// Events are its zero crossings (e.g. height above the ellipsoid, or the z coordinate for nodes).
using EventFunction = std::function<double(double dt, const StateVector& state)>;

/// An event found by NumericalPropagator::propagateToEvent.
struct PropagationEvent {
    double      dt{0.0}; ///< seconds after the epoch
    StateVector state;   ///< position (m) and velocity (m/s) at dt
};

/// Integrates the equations of motion with RKF7(8).
class NumericalPropagator {
public:
    /// Start from position (m) and velocity (m/s) at `epoch`, in the inertial frame the force
    /// models assume (all outputs are in the same frame). Throws astro::Error if `forces` is empty
    /// or holds a null pointer.
    NumericalPropagator(const Time&                                    epoch,
                        const StateVector&                             initial,
                        std::vector<std::shared_ptr<const ForceModel>> forces,
                        const AdaptiveOptions&                         options = {});

    /// Start from a 12-component state (position, velocity, acceleration, jerk). The trajectory
    /// starts from position and velocity; `mode` says what the input acceleration and jerk do.
    /// With HigherOrderMode::UnmodeledForce an UnmodeledAcceleration holding the residuals is
    /// added to `forces`.
    NumericalPropagator(const Time&                                    epoch,
                        const StateVector12&                           initial,
                        std::vector<std::shared_ptr<const ForceModel>> forces,
                        HigherOrderMode                                mode    = HigherOrderMode::ModelOnly,
                        const AdaptiveOptions&                         options = {});

    /// Point mass + J2 + J4 for a state in the celestial frame `frame` (GCRS, J2000, MOD, TOD, TEME
    /// or CIRS) at `epoch`: the zonal field is aligned with the true pole at the epoch (in GCRS the
    /// pole is ~0.15 deg off z in 2026, which would cost ~270 m over 4 h in LEO). The frame itself
    /// is treated as inertial, frozen at the epoch; for of-date frames (MOD, TOD, TEME, CIRS) that
    /// drifts from the frame of each output time by ~0.2 m/h in LEO. EstimateTracker converts
    /// outputs to the frame of their own time.
    [[nodiscard]] static NumericalPropagator
    j2j4(const Time& epoch, const StateVector& initial, Frame frame, const AdaptiveOptions& options = {});
    /// As above, from a 12-component state.
    [[nodiscard]] static NumericalPropagator j2j4(const Time&            epoch,
                                                  const StateVector12&   initial,
                                                  Frame                  frame,
                                                  HigherOrderMode        mode    = HigherOrderMode::ModelOnly,
                                                  const AdaptiveOptions& options = {});

    /// Point mass + J2 + J4 (the zonal J3 is omitted) in a frame whose z axis is the true pole
    /// (TOD, TEME, CIRS). For GCRS or J2000 states use the overload taking a Frame.
    [[nodiscard]] static NumericalPropagator j2j4(const Time&        epoch,
                                                  const StateVector& initial,
                                                  double             mu = constants::gmEarth,
                                                  double             re = constants::earthRadiusEgm2008,
                                                  double             J2 = constants::J2,
                                                  double             J4 = constants::J4);
    /// As above, from a 12-component state.
    [[nodiscard]] static NumericalPropagator j2j4(const Time&          epoch,
                                                  const StateVector12& initial,
                                                  HigherOrderMode      mode = HigherOrderMode::ModelOnly,
                                                  double               mu   = constants::gmEarth,
                                                  double               re   = constants::earthRadiusEgm2008,
                                                  double               J2   = constants::J2,
                                                  double               J4   = constants::J4);

    // Every propagation call takes an optional `stats` out-parameter for the integrator's step
    // counts. The propagator itself is immutable, so all member functions may be called
    // concurrently on a shared instance (provided its force models are thread-safe; the library's
    // are).

    /// State dt seconds after the epoch (negative dt propagates backward).
    [[nodiscard]] StateVector propagate(double dt, AdaptiveStats* stats = nullptr) const;
    /// State at time t (any scale).
    [[nodiscard]] StateVector propagate(const Time& t, AdaptiveStats* stats = nullptr) const {
        return propagate(t - epoch_, stats);
    }
    /// States at several offsets (any order), integrating once in each direction.
    [[nodiscard]] std::vector<StateVector> propagate(const std::vector<double>& dts,
                                                     AdaptiveStats*             stats = nullptr) const;

    /// 12-component state dt seconds after the epoch; acceleration and jerk come from the force
    /// models (including any unmodeled acceleration).
    [[nodiscard]] StateVector12 propagate12(double dt, AdaptiveStats* stats = nullptr) const;
    /// 12-component state at time t (any scale).
    [[nodiscard]] StateVector12 propagate12(const Time& t, AdaptiveStats* stats = nullptr) const {
        return propagate12(t - epoch_, stats);
    }
    /// 12-component states at several offsets (any order), integrating once in each direction.
    [[nodiscard]] std::vector<StateVector12> propagate12(const std::vector<double>& dts,
                                                         AdaptiveStats*             stats = nullptr) const;

    /// States at dt = start + k * step for k = 0 .. out.size() - 1 (step != 0), written into `out`.
    ///
    /// Use this for many equally spaced outputs. The integrator takes its own (error-controlled)
    /// steps instead of stopping on every output time, and outputs between steps come from a
    /// degree-7 Hermite polynomial that matches position, velocity, acceleration and jerk at both
    /// ends of the step. The cost therefore hardly depends on `step`. Positions agree with
    /// propagate(dts) to about the integrator tolerance (micrometres in LEO with the defaults).
    ///
    /// The grid may lie anywhere relative to the epoch (dt = 0) and run in either direction. The
    /// integration only ever moves away from the epoch: to the grid's nearest point and on through
    /// the rest or, if the grid spans the epoch, outwards from it once in each direction. So each
    /// output carries the error of the direct path to its own time, and an output at dt = 0 is the
    /// initial state exactly.
    void propagateGrid(double start, double step, std::span<StateVector> out, AdaptiveStats* stats = nullptr) const;
    /// As propagateGrid, with acceleration and jerk from the force models at each output.
    void propagateGrid12(double start, double step, std::span<StateVector12> out, AdaptiveStats* stats = nullptr) const;
    /// As propagateGrid, returning `count` states.
    [[nodiscard]] std::vector<StateVector>
    propagateGrid(double start, double step, std::size_t count, AdaptiveStats* stats = nullptr) const;
    /// As propagateGrid12, returning `count` states.
    [[nodiscard]] std::vector<StateVector12>
    propagateGrid12(double start, double step, std::size_t count, AdaptiveStats* stats = nullptr) const;

    /// States and covariances at dt = start + k * step for k = 0 .. out.size() - 1 (step != 0).
    ///
    /// P0 is the covariance of the initial state: position and velocity in the propagation frame
    /// (m^2, m^2/s, m^2/s^2). The state transition matrix Phi(t, epoch) is integrated with the
    /// state from each force model's partials(), and P(t) = Phi P0 Phi^T + Q(t).
    ///
    /// The integrator takes exactly the steps of propagateGrid, so the states are identical to it.
    /// Between steps Phi and Q are interpolated with quintic Hermite polynomials (values, first and
    /// second derivatives at both ends of the step): relative error about 2e-9 at the default
    /// tolerance and 4e-7 at relTol = 1e-9.
    ///
    /// Q(t) is the process noise from white, isotropic acceleration noise with power spectral
    /// density accelerationNoisePsd (m^2/s^3), standing in for forces the model leaves out (drag,
    /// solar pressure, higher harmonics). 0 gives the purely deterministic propagation. Q is the
    /// noise accumulated between the epoch and t, so it adds uncertainty in both directions:
    /// propagating backward does not shrink the covariance.
    ///
    /// As in propagateGrid the integration only moves away from the epoch, so Q only accumulates:
    /// an output at dt = 0 has Phi = I and Q = 0 exactly (its covariance is P0), and on a grid that
    /// spans the epoch every Q is positive definite.
    void propagateGridWithCovariance(double                     start,
                                     double                     step,
                                     const Mat6&                P0,
                                     std::span<StateCovariance> out,
                                     double                     accelerationNoisePsd = 0.0,
                                     AdaptiveStats*             stats                = nullptr) const;
    /// As above, returning `count` states and covariances.
    [[nodiscard]] std::vector<StateCovariance> propagateGridWithCovariance(double         start,
                                                                           double         step,
                                                                           std::size_t    count,
                                                                           const Mat6&    P0,
                                                                           double         accelerationNoisePsd = 0.0,
                                                                           AdaptiveStats* stats = nullptr) const;

    /// States and state transition matrices Phi(t, epoch) at dt = start + k * step, as for
    /// propagateGridWithCovariance. states and stm must have the same size.
    void propagateGridStm(double                 start,
                          double                 step,
                          std::span<StateVector> states,
                          std::span<Mat6>        stm,
                          AdaptiveStats*         stats = nullptr) const;

    /// First zero crossing of the event function g, in the given direction, between the epoch and
    /// maxDt (negative: search backward in time). nullopt if there is none.
    ///
    /// The integrator takes its usual error-controlled steps. g is sampled at the end and at three
    /// interior points of each step (from the degree-7 Hermite interpolant used by propagateGrid),
    /// and a sign change is refined with Brent's method on the interpolant to `timeTolerance`
    /// seconds; the integration then stops. A zero of g at the epoch itself is not an event, so a
    /// trajectory that starts on the surface finds the impact, not the launch. Crossings that
    /// enter and leave within a quarter step can be missed (steps are ~100 s in LEO).
    [[nodiscard]] std::optional<PropagationEvent> propagateToEvent(const EventFunction& g,
                                                                   double               maxDt,
                                                                   Crossing             crossing      = Crossing::Any,
                                                                   double               timeTolerance = 1e-6,
                                                                   AdaptiveStats*       stats         = nullptr) const;

    /// Acceleration and jerk of the force models at state s, t seconds after the epoch.
    [[nodiscard]] StateVector12 withDerivatives(double t, const StateVector& s) const;

    /// Input acceleration minus the force models' at the epoch (zero for a 6-component start).
    [[nodiscard]] const Vec3& accelerationResidual() const noexcept { return da0_; }
    /// Input jerk minus the force models' at the epoch (zero for a 6-component start).
    [[nodiscard]] const Vec3& jerkResidual() const noexcept { return dj0_; }

    /// Epoch of the initial state; dt arguments are seconds since it.
    [[nodiscard]] const Time& epoch() const noexcept { return epoch_; }

private:
    [[nodiscard]] OdeState<6> derivative(double t, const OdeState<6>& y) const;
    template <class Emit>
    void grid(double start, double step, std::size_t count, AdaptiveStats* stats, Emit&& emit) const;
    template <std::size_t N, class Emit>
    void
    variationalGrid(double start, double step, std::size_t count, double q, AdaptiveStats* stats, Emit&& emit) const;

    Time                                           epoch_;
    StateVector                                    s0_;
    std::vector<std::shared_ptr<const ForceModel>> forces_;
    AdaptiveOptions                                opt_;
    Vec3                                           da0_;
    Vec3                                           dj0_;
};

} // namespace astro
