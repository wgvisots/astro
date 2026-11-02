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
#include "astro/core/vec.h"
#include "astro/ephemeris/ephemeris.h"
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
    [[nodiscard]] Vec3   jerk(double t, const Vec3& r, const Vec3& v, const Vec3& totalAcceleration) const override;
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

private:
    Vec3 a0_;
    Vec3 j0_;
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

    /// Convenience: point mass + J2 + J4 (the zonal J3 is omitted), in a frame with z along the pole.
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

    /// State dt seconds after the epoch (negative dt propagates backward).
    [[nodiscard]] StateVector propagate(double dt) const;
    [[nodiscard]] StateVector propagate(const Time& t) const { return propagate(t - epoch_); }
    /// States at several offsets (any order), integrating once in each direction.
    [[nodiscard]] std::vector<StateVector> propagate(const std::vector<double>& dts) const;

    /// 12-component state dt seconds after the epoch; acceleration and jerk come from the force
    /// models (including any unmodeled acceleration).
    [[nodiscard]] StateVector12              propagate12(double dt) const;
    [[nodiscard]] StateVector12              propagate12(const Time& t) const { return propagate12(t - epoch_); }
    [[nodiscard]] std::vector<StateVector12> propagate12(const std::vector<double>& dts) const;

    /// States at dt = start + k * step for k = 0 .. out.size() - 1 (step != 0), written into `out`.
    ///
    /// Use this for many equally spaced outputs. The integrator takes its own (error-controlled)
    /// steps instead of stopping on every output time, and outputs between steps come from a
    /// degree-7 Hermite polynomial that matches position, velocity, acceleration and jerk at both
    /// ends of the step. The cost therefore hardly depends on `step`. Positions agree with
    /// propagate(dts) to about the integrator tolerance (micrometres in LEO with the defaults).
    void propagateGrid(double start, double step, std::span<StateVector> out) const;
    /// As propagateGrid, with acceleration and jerk from the force models at each output.
    void propagateGrid12(double start, double step, std::span<StateVector12> out) const;
    [[nodiscard]] std::vector<StateVector>   propagateGrid(double start, double step, std::size_t count) const;
    [[nodiscard]] std::vector<StateVector12> propagateGrid12(double start, double step, std::size_t count) const;

    /// Acceleration and jerk of the force models at state s, t seconds after the epoch.
    [[nodiscard]] StateVector12 withDerivatives(double t, const StateVector& s) const;

    /// Input acceleration minus the force models' at the epoch (zero for a 6-component start).
    [[nodiscard]] const Vec3& accelerationResidual() const noexcept { return da0_; }
    /// Input jerk minus the force models' at the epoch (zero for a 6-component start).
    [[nodiscard]] const Vec3& jerkResidual() const noexcept { return dj0_; }

    [[nodiscard]] const Time&          epoch() const noexcept { return epoch_; }
    [[nodiscard]] const AdaptiveStats& lastStats() const noexcept { return stats_; }

private:
    [[nodiscard]] OdeState<6> derivative(double t, const OdeState<6>& y) const;
    template <class Emit>
    void grid(double start, double step, std::size_t count, Emit&& emit) const;

    Time                                           epoch_;
    StateVector                                    s0_;
    std::vector<std::shared_ptr<const ForceModel>> forces_;
    AdaptiveOptions                                opt_;
    Vec3                                           da0_;
    Vec3                                           dj0_;
    mutable AdaptiveStats                          stats_;
};

} // namespace astro
