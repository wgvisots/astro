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
#include <vector>

#include "astro/core/constants.h"
#include "astro/core/vec.h"
#include "astro/ephemeris/ephemeris.h"
#include "astro/orbit/integrators.h"
#include "astro/time/time.h"

/// Numerical orbit propagation with pluggable force models.
namespace astro {

/// A perturbing (or central) acceleration. `t` is seconds since the propagator epoch;
/// r, v are in the propagator's inertial frame.
class ForceModel {
public:
    virtual ~ForceModel() = default;
    [[nodiscard]] virtual Vec3 acceleration(double t, const Vec3& r, const Vec3& v) const = 0;
};

/// Central body gravity: point mass plus zonal harmonics J2..Jn.
///
/// `toBodyFrame` rotates the propagation frame into one whose z axis is the
/// rotation axis (e.g. the GCRS -> CIRS matrix at epoch); identity treats the
/// propagation frame's z axis as the pole.
class GravityField final : public ForceModel {
public:
    /// J indexed by degree (J[2] = J2, ...). Empty -> point mass only.
    explicit GravityField(std::vector<double> J = {0.0, 0.0, constants::J2, constants::J3, constants::J4},
                          double mu = constants::gmEarth, double re = constants::earthRadiusEgm2008,
                          const Mat3& toBodyFrame = Mat3::identity());
    [[nodiscard]] Vec3 acceleration(double t, const Vec3& r, const Vec3& v) const override;
    [[nodiscard]] double mu() const noexcept { return mu_; }

private:
    std::vector<double> J_;
    double mu_;
    double re_;
    Mat3 R_;
    bool identity_;
};

/// Point-mass perturbation from a third body (e.g. Sun, Moon) using an ephemeris.
/// The propagation frame is assumed to be the GCRS.
class ThirdBodyForce final : public ForceModel {
public:
    ThirdBodyForce(Body body, double muBody, std::shared_ptr<const Ephemeris> ephemeris, const Time& epoch);
    [[nodiscard]] Vec3 acceleration(double t, const Vec3& r, const Vec3& v) const override;

private:
    Body body_;
    double mu_;
    std::shared_ptr<const Ephemeris> eph_;
    Time epoch_;
};

/// Integrates the equations of motion with RKF7(8).
class NumericalPropagator {
public:
    NumericalPropagator(const Time& epoch, const StateVector& initial,
                        std::vector<std::shared_ptr<const ForceModel>> forces, const AdaptiveOptions& options = {});

    /// Convenience: point mass + J2 + J4 (the zonal J3 is omitted), in a frame with z along the pole.
    [[nodiscard]] static NumericalPropagator j2j4(const Time& epoch, const StateVector& initial,
                                                  double mu = constants::gmEarth,
                                                  double re = constants::earthRadiusEgm2008, double J2 = constants::J2,
                                                  double J4 = constants::J4);

    /// State dt seconds after the epoch (negative dt propagates backward).
    [[nodiscard]] StateVector propagate(double dt) const;
    [[nodiscard]] StateVector propagate(const Time& t) const { return propagate(t - epoch_); }
    /// States at several offsets (any order), integrating once in each direction.
    [[nodiscard]] std::vector<StateVector> propagate(const std::vector<double>& dts) const;

    [[nodiscard]] const Time& epoch() const noexcept { return epoch_; }
    [[nodiscard]] const AdaptiveStats& lastStats() const noexcept { return stats_; }

private:
    [[nodiscard]] OdeState<6> derivative(double t, const OdeState<6>& y) const;

    Time epoch_;
    StateVector s0_;
    std::vector<std::shared_ptr<const ForceModel>> forces_;
    AdaptiveOptions opt_;
    mutable AdaptiveStats stats_;
};

} // namespace astro
