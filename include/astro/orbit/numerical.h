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

#include "astro/core/mat6.h"
#include "astro/core/vec.h"
#include "astro/ephemeris/ephemeris.h"
#include "astro/frames/transform.h"
#include "astro/orbit/central_body.h"
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
    /// concurrently if the propagator using it is shared between threads. An implementation may
    /// throw (ThirdBodyForce does, where its ephemeris has no data); the exception passes through
    /// the propagation that called it.
    [[nodiscard]] virtual Vec3 acceleration(double t, const Vec3& r, const Vec3& v) const = 0;

    /// Time derivative of acceleration() along the trajectory, d/dt a(t, r(t), v(t)), given the
    /// total acceleration acting on the body. The default is a central difference along
    /// (1, v, totalAcceleration) with a step of 1e-5 |r| / |v| seconds (about 0.01 s in LEO), which
    /// is accurate to roughly 1e-9 relative for orbital forces. Override it for an exact form.
    /// The default throws what acceleration() throws.
    [[nodiscard]] virtual Vec3 jerk(double t, const Vec3& r, const Vec3& v, const Vec3& totalAcceleration) const;

    /// acceleration(t, r, v), and its partial derivatives with respect to position (dadr) and
    /// velocity (dadv), for the state transition matrix. The default uses central differences
    /// (12 extra acceleration evaluations); override it with an exact form where one exists.
    /// The default throws what acceleration() throws.
    virtual Vec3 partials(double t, const Vec3& r, const Vec3& v, Mat3& dadr, Mat3& dadv) const;
};

/// Central body gravity: point mass plus zonal harmonics J2..Jn.
///
/// `toBodyFrame` rotates the propagation frame into one whose z axis is the
/// rotation axis (e.g. the GCRS -> CIRS matrix at epoch); identity treats the
/// propagation frame's z axis as the pole.
class GravityField final : public ForceModel {
public:
    /// The point mass and every zonal term of `body` (an empty body.J: point mass only).
    explicit GravityField(const CentralBody& body = CentralBody::earth(), const Mat3& toBodyFrame = Mat3::identity());
    /// Point-mass plus zonal acceleration (m/s^2) in the propagation frame; t and v are unused.
    [[nodiscard]] Vec3 acceleration(double t, const Vec3& r, const Vec3& v) const override;
    /// Exact (point mass and zonal terms).
    [[nodiscard]] Vec3 jerk(double t, const Vec3& r, const Vec3& v, const Vec3& totalAcceleration) const override;
    /// Exact gravity gradient; no velocity dependence.
    Vec3 partials(double t, const Vec3& r, const Vec3& v, Mat3& dadr, Mat3& dadv) const override;
    /// The central body of the field.
    [[nodiscard]] const CentralBody& body() const noexcept { return body_; }

private:
    CentralBody body_;
    Mat3        R_;
    bool        identity_;
};

/// Point-mass perturbation from a third body (e.g. Sun, Moon) using an ephemeris.
/// The propagation frame is assumed to be the GCRS.
class ThirdBodyForce final : public ForceModel {
public:
    /// Perturbation by `body` with gravitational parameter `muBody` (m^3/s^2); its geocentric
    /// position comes from `ephemeris` at `epoch` + t. Throws astro::Error if `ephemeris` is null.
    ThirdBodyForce(Body body, double muBody, std::shared_ptr<const Ephemeris> ephemeris, const Time& epoch);
    /// Tidal acceleration (m/s^2) on a satellite at geocentric r: the body's pull on the satellite
    /// minus its pull on the Earth. v is unused. Throws astro::Error if the ephemeris has no data
    /// for the body at `epoch` + t (see Ephemeris::state).
    [[nodiscard]] Vec3 acceleration(double t, const Vec3& r, const Vec3& v) const override;
    /// Exact tidal gradient; no velocity dependence. Throws astro::Error as acceleration() does.
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
    /// added to `forces`. Throws astro::Error as the constructor above does. An exception thrown
    /// by a force model at the epoch passes through.
    NumericalPropagator(const Time&                                    epoch,
                        const StateVector12&                           initial,
                        std::vector<std::shared_ptr<const ForceModel>> forces,
                        HigherOrderMode                                mode    = HigherOrderMode::ModelOnly,
                        const AdaptiveOptions&                         options = {});

    /// Point mass + J2 + J4 (CentralBody::earthJ2J4()) for a state in the celestial frame `frame`
    /// (GCRS, J2000, MOD, TOD, TEME or CIRS) at `epoch`: the zonal field is aligned with the true
    /// pole at the epoch (in GCRS the pole is ~0.15 deg off z in 2026, which would cost ~270 m over 4 h in LEO). The frame itself
    /// is treated as inertial, frozen at the epoch; for of-date frames (MOD, TOD, TEME, CIRS) that
    /// drifts from the frame of each output time by ~0.2 m/h in LEO. EstimateTracker converts
    /// outputs to the frame of their own time. Throws astro::Error as poleOfDateMatrix does: for an
    /// Earth-fixed `frame` (TIRS, PEF, ITRS), and if `epoch` cannot be converted to TT.
    [[nodiscard]] static NumericalPropagator
    j2j4(const Time& epoch, const StateVector& initial, Frame frame, const AdaptiveOptions& options = {});
    /// As above, from a 12-component state. Throws astro::Error as the overload above does.
    [[nodiscard]] static NumericalPropagator j2j4(const Time&            epoch,
                                                  const StateVector12&   initial,
                                                  Frame                  frame,
                                                  HigherOrderMode        mode    = HigherOrderMode::ModelOnly,
                                                  const AdaptiveOptions& options = {});

    /// Point mass and every zonal term of `body`, for a state in the celestial frame `frame` (GCRS,
    /// J2000, MOD, TOD, TEME or CIRS) at `epoch`: the zonal field is aligned with the true pole at the
    /// epoch, frozen there, exactly as for j2j4 (which is zonal with CentralBody::earthJ2J4()).
    /// Throws astro::Error as j2j4 does: for an Earth-fixed `frame` (TIRS, PEF, ITRS), and if
    /// `epoch` cannot be converted to TT.
    [[nodiscard]] static NumericalPropagator zonal(const Time&            epoch,
                                                   const StateVector&     initial,
                                                   const CentralBody&     body,
                                                   Frame                  frame,
                                                   const AdaptiveOptions& options = {});
    /// As above, from a 12-component state. Throws astro::Error as the overload above does.
    [[nodiscard]] static NumericalPropagator zonal(const Time&            epoch,
                                                   const StateVector12&   initial,
                                                   const CentralBody&     body,
                                                   Frame                  frame,
                                                   HigherOrderMode        mode    = HigherOrderMode::ModelOnly,
                                                   const AdaptiveOptions& options = {});

    // Every propagation call takes an optional `stats` out-parameter for the integrator's step
    // counts. The propagator itself is immutable, so all member functions may be called
    // concurrently on a shared instance (provided its force models are thread-safe; the library's
    // are).

    /// State dt seconds after the epoch (negative dt propagates backward). Throws astro::Error
    /// when the integration fails, as integrateRkf78Steps does: for a dt that is not finite,
    /// unusable AdaptiveOptions, a state or an acceleration that stops being finite, an error that
    /// stays too large at AdaptiveOptions::minStep, or more than AdaptiveOptions::maxSteps steps.
    /// An exception thrown by a force model passes through.
    [[nodiscard]] StateVector propagate(double dt, AdaptiveStats* stats = nullptr) const;
    /// State at time t (any scale). Throws astro::Error as propagate(dt) does, and if t cannot be
    /// subtracted from the epoch (one of the two in UT1, or a UTC time before 1960; see Time).
    [[nodiscard]] StateVector propagate(const Time& t, AdaptiveStats* stats = nullptr) const {
        return propagate(t - epoch_, stats);
    }
    /// States at several offsets (any order), integrating once in each direction. Throws
    /// astro::Error as propagate(dt) does, for any of the offsets.
    [[nodiscard]] std::vector<StateVector> propagate(const std::vector<double>& dts,
                                                     AdaptiveStats*             stats = nullptr) const;

    /// 12-component state dt seconds after the epoch; acceleration and jerk come from the force
    /// models (including any unmodeled acceleration). Throws astro::Error as propagate(dt) does.
    [[nodiscard]] StateVector12 propagate12(double dt, AdaptiveStats* stats = nullptr) const;
    /// 12-component state at time t (any scale). Throws astro::Error as propagate(t) does.
    [[nodiscard]] StateVector12 propagate12(const Time& t, AdaptiveStats* stats = nullptr) const {
        return propagate12(t - epoch_, stats);
    }
    /// 12-component states at several offsets (any order), integrating once in each direction.
    /// Throws astro::Error as propagate(dt) does, for any of the offsets.
    [[nodiscard]] std::vector<StateVector12> propagate12(const std::vector<double>& dts,
                                                         AdaptiveStats*             stats = nullptr) const;

    /// States at dt = start + k * step for k = 0 .. out.size() - 1 (step != 0), written into `out`.
    ///
    /// Use this for many equally spaced outputs. The integrator takes its own (error-controlled)
    /// steps instead of stopping on every output time, and outputs between steps come from a
    /// degree-7 Hermite polynomial that matches position, velocity, acceleration and jerk at both
    /// ends of the step. The cost therefore hardly depends on `step`. Positions agree with
    /// propagate(dts), which takes different steps, to about the accuracy of the integration
    /// itself: 0.15 mm after 4 h in LEO with the defaults, 0.2 m at relTol = 1e-9.
    ///
    /// The grid may lie anywhere relative to the epoch (dt = 0) and run in either direction. The
    /// integration only ever moves away from the epoch: to the grid's nearest point and on through
    /// the rest or, if the grid spans the epoch, outwards from it once in each direction. So each
    /// output carries the error of the direct path to its own time, and an output at dt = 0 is the
    /// initial state exactly.
    ///
    /// Throws astro::Error if `start` is not finite or `step` is zero or not finite (an empty `out`
    /// is accepted whatever they are), and as propagate(dt) does.
    void propagateGrid(double start, double step, std::span<StateVector> out, AdaptiveStats* stats = nullptr) const;
    /// As propagateGrid, with acceleration and jerk from the force models at each output. Throws
    /// astro::Error as propagateGrid does.
    void propagateGrid12(double start, double step, std::span<StateVector12> out, AdaptiveStats* stats = nullptr) const;
    /// As propagateGrid, returning `count` states. Throws astro::Error as propagateGrid does.
    [[nodiscard]] std::vector<StateVector>
    propagateGrid(double start, double step, std::size_t count, AdaptiveStats* stats = nullptr) const;
    /// As propagateGrid12, returning `count` states. Throws astro::Error as propagateGrid does.
    [[nodiscard]] std::vector<StateVector12>
    propagateGrid12(double start, double step, std::size_t count, AdaptiveStats* stats = nullptr) const;

    /// State dt seconds after the epoch and the state transition matrix Phi(dt, epoch) = d(state at
    /// dt) / d(initial state), integrated with the state as in propagateGridStm: the same steps and
    /// the same result as a grid of one output at dt. Phi(0) is the identity. Throws astro::Error as
    /// propagate(dt) does.
    [[nodiscard]] StateTransition propagateStm(double dt, AdaptiveStats* stats = nullptr) const;
    /// State and covariance dt seconds after the epoch: P = Phi P0 Phi^T + Q(dt), as in
    /// propagateGridWithCovariance (the same steps and the same result as a grid of one output at dt).
    /// P0 is the covariance of the initial state in the propagation frame (m^2, m^2/s, m^2/s^2);
    /// accelerationNoisePsd is the power spectral density of the white acceleration noise (m^2/s^3),
    /// 0 for none. Throws astro::Error as propagate(dt) does, and if accelerationNoisePsd is negative
    /// or not finite.
    [[nodiscard]] StateCovariance propagateWithCovariance(double         dt,
                                                          const Mat6&    P0,
                                                          double         accelerationNoisePsd = 0.0,
                                                          AdaptiveStats* stats                = nullptr) const;

    /// States and covariances at dt = start + k * step for k = 0 .. out.size() - 1 (step != 0).
    ///
    /// P0 is the covariance of the initial state: position and velocity in the propagation frame
    /// (m^2, m^2/s, m^2/s^2). The state transition matrix Phi(t, epoch) is integrated with the
    /// state from each force model's partials(), and P(t) = Phi P0 Phi^T + Q(t).
    ///
    /// The integrator takes exactly the steps of propagateGrid, so the states are identical to it.
    /// Between steps Phi and Q are interpolated with quintic Hermite polynomials (values, first and
    /// second derivatives at both ends of the step). Against stopping exactly at each output,
    /// over 4 h at 10 s: about 1e-8 relative at the default tolerance and 1e-6 at relTol = 1e-9
    /// for near-circular orbits, ten times that for a 1000 x 20000 km orbit.
    ///
    /// Q(t) is the process noise from white, isotropic acceleration noise (the continuous
    /// white-noise acceleration model; Bar-Shalom, Li & Kirubarajan 2001) with power spectral
    /// density accelerationNoisePsd (m^2/s^3), standing in for forces the model leaves out (drag,
    /// solar pressure, higher harmonics). 0 gives the purely deterministic propagation. Q is the
    /// noise accumulated between the epoch and t, so it adds uncertainty in both directions:
    /// propagating backward does not shrink the covariance.
    ///
    /// As in propagateGrid the integration only moves away from the epoch, so Q only accumulates:
    /// an output at dt = 0 has Phi = I and Q = 0 exactly (its covariance is P0), and on a grid that
    /// spans the epoch every Q is positive definite.
    ///
    /// Throws astro::Error as propagateGrid does, and if accelerationNoisePsd is negative or not
    /// finite.
    void propagateGridWithCovariance(double                     start,
                                     double                     step,
                                     const Mat6&                P0,
                                     std::span<StateCovariance> out,
                                     double                     accelerationNoisePsd = 0.0,
                                     AdaptiveStats*             stats                = nullptr) const;
    /// As above, returning `count` states and covariances. Throws astro::Error as the overload
    /// above does.
    [[nodiscard]] std::vector<StateCovariance> propagateGridWithCovariance(double         start,
                                                                           double         step,
                                                                           std::size_t    count,
                                                                           const Mat6&    P0,
                                                                           double         accelerationNoisePsd = 0.0,
                                                                           AdaptiveStats* stats = nullptr) const;

    /// States and state transition matrices Phi(t, epoch) at dt = start + k * step, as for
    /// propagateGridWithCovariance. states and stm must have the same size. Throws astro::Error if
    /// they differ in size, and as propagateGrid does.
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
    ///
    /// Throws astro::Error if g is empty, if maxDt is not finite, if timeTolerance is not positive
    /// and finite, if g returns a value that is not finite, and as propagate(dt) does. An exception
    /// thrown by g passes through.
    [[nodiscard]] std::optional<PropagationEvent> propagateToEvent(const EventFunction& g,
                                                                   double               maxDt,
                                                                   Crossing             crossing      = Crossing::Any,
                                                                   double               timeTolerance = 1e-6,
                                                                   AdaptiveStats*       stats         = nullptr) const;

    /// Acceleration and jerk of the force models at state s, t seconds after the epoch. An
    /// exception thrown by a force model passes through.
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
