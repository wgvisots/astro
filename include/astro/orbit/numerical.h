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

#include <memory>
#include <span>
#include <vector>

#include "astro/core/constants.h"
#include "astro/core/mat6.h"
#include "astro/core/vec.h"
#include "astro/ephemeris/ephemeris.h"
#include "astro/frames/transform.h"
#include "astro/orbit/integrators.h"
#include "astro/orbit/state12.h"
#include "astro/time/time.h"

/// Numerical orbit propagation with pluggable force models.
namespace astro {

/// A perturbing (or central) acceleration. `t` is seconds since the propagator epoch;
/// r, v are in the propagator's inertial frame.
class ForceModel {
public:
    virtual ~ForceModel()                                                                 = default;
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
    [[nodiscard]] Vec3 acceleration(double t, const Vec3& r, const Vec3& v) const override;
    /// Exact (point mass and zonal terms).
    [[nodiscard]] Vec3 jerk(double t, const Vec3& r, const Vec3& v, const Vec3& totalAcceleration) const override;
    /// Exact gravity gradient; no velocity dependence.
    Vec3                 partials(double t, const Vec3& r, const Vec3& v, Mat3& dadr, Mat3& dadv) const override;
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
    ThirdBodyForce(Body body, double muBody, std::shared_ptr<const Ephemeris> ephemeris, const Time& epoch);
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
    UnmodeledAcceleration(const Vec3& a0, const Vec3& j0) noexcept
        : a0_(a0)
        , j0_(j0) {}
    [[nodiscard]] Vec3 acceleration(double t, const Vec3& r, const Vec3& v) const override;
    [[nodiscard]] Vec3 jerk(double t, const Vec3& r, const Vec3& v, const Vec3& totalAcceleration) const override;
    /// Zero: the unmodeled acceleration depends on time only.
    Vec3 partials(double t, const Vec3& r, const Vec3& v, Mat3& dadr, Mat3& dadv) const override;

private:
    Vec3 a0_;
    Vec3 j0_;
};

/// A state and its 6x6 position/velocity covariance at one time.
struct StateCovariance {
    StateVector state;
    Mat6        covariance;
};

/// Integrates the equations of motion with RKF7(8).
class NumericalPropagator {
public:
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
    [[nodiscard]] StateVector propagate(const Time& t, AdaptiveStats* stats = nullptr) const {
        return propagate(t - epoch_, stats);
    }
    /// States at several offsets (any order), integrating once in each direction.
    [[nodiscard]] std::vector<StateVector> propagate(const std::vector<double>& dts,
                                                     AdaptiveStats*             stats = nullptr) const;

    /// 12-component state dt seconds after the epoch; acceleration and jerk come from the force
    /// models (including any unmodeled acceleration).
    [[nodiscard]] StateVector12 propagate12(double dt, AdaptiveStats* stats = nullptr) const;
    [[nodiscard]] StateVector12 propagate12(const Time& t, AdaptiveStats* stats = nullptr) const {
        return propagate12(t - epoch_, stats);
    }
    [[nodiscard]] std::vector<StateVector12> propagate12(const std::vector<double>& dts,
                                                         AdaptiveStats*             stats = nullptr) const;

    /// States at dt = start + k * step for k = 0 .. out.size() - 1 (step != 0), written into `out`.
    ///
    /// Use this for many equally spaced outputs. The integrator takes its own (error-controlled)
    /// steps instead of stopping on every output time, and outputs between steps come from a
    /// degree-7 Hermite polynomial that matches position, velocity, acceleration and jerk at both
    /// ends of the step. The cost therefore hardly depends on `step`. Positions agree with
    /// propagate(dts) to about the integrator tolerance (micrometres in LEO with the defaults).
    void propagateGrid(double start, double step, std::span<StateVector> out, AdaptiveStats* stats = nullptr) const;
    /// As propagateGrid, with acceleration and jerk from the force models at each output.
    void propagateGrid12(double start, double step, std::span<StateVector12> out, AdaptiveStats* stats = nullptr) const;
    [[nodiscard]] std::vector<StateVector>
    propagateGrid(double start, double step, std::size_t count, AdaptiveStats* stats = nullptr) const;
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
    /// solar pressure, higher harmonics). 0 gives the purely deterministic propagation.
    void                                       propagateGridWithCovariance(double                     start,
                                                                           double                     step,
                                                                           const Mat6&                P0,
                                                                           std::span<StateCovariance> out,
                                                                           double                     accelerationNoisePsd = 0.0,
                                                                           AdaptiveStats*             stats = nullptr) const;
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

    /// Acceleration and jerk of the force models at state s, t seconds after the epoch.
    [[nodiscard]] StateVector12 withDerivatives(double t, const StateVector& s) const;

    /// Input acceleration minus the force models' at the epoch (zero for a 6-component start).
    [[nodiscard]] const Vec3& accelerationResidual() const noexcept { return da0_; }
    /// Input jerk minus the force models' at the epoch (zero for a 6-component start).
    [[nodiscard]] const Vec3& jerkResidual() const noexcept { return dj0_; }

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
