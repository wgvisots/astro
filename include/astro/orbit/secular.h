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

/// Analytic secular (mean-element) drift due to the Earth's zonal harmonics,
/// after Brouwer (1959): first order in J2, plus J2^2 and J4.
namespace astro {

struct ZonalParameters {
    double mu{constants::gmEarth};
    double re{constants::earthRadiusEgm2008};
    double J2{constants::J2};
    double J4{constants::J4};
};

enum class SecularOrder {
    J2,         ///< first-order J2 only
    J2SquaredJ4 ///< J2 + J2^2 + J4 (Brouwer)
};

BOOST_DESCRIBE_ENUM(SecularOrder, J2, J2SquaredJ4)

/// Writes "SecularOrder::<enumerator>", or "SecularOrder::INVALID_VALUE" for any other value.
std::ostream& operator<<(std::ostream& os, SecularOrder value);

/// Secular rates (rad/s) of the mean elements.
struct SecularRates {
    double raanDot{0.0};
    double argpDot{0.0};
    double meanAnomalyDot{0.0}; ///< total rate including the mean motion
};

/// Rates for mean semi-major axis a, eccentricity e and inclination i.
[[nodiscard]] SecularRates secularRates(double                 a,
                                        double                 e,
                                        double                 i,
                                        const ZonalParameters& zp    = {},
                                        SecularOrder           order = SecularOrder::J2SquaredJ4) noexcept;

/// Propagate mean elements by dt seconds with the secular rates (a, e, i constant).
/// The true anomaly is advanced through the mean anomaly.
[[nodiscard]] KeplerianElements propagateSecular(const KeplerianElements& mean,
                                                 double                   dt,
                                                 const ZonalParameters&   zp    = {},
                                                 SecularOrder             order = SecularOrder::J2SquaredJ4);

/// propagateSecular for an inertial state whose z axis is the Earth's pole. The osculating
/// elements of s0 are treated as mean elements (no short-period correction), so expect errors at the
/// short-period level (kilometres in LEO); use NumericalPropagator::j2j4 when that matters.
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
class SecularPropagator {
public:
    /// From mean elements (elliptic: e < 1).
    explicit SecularPropagator(const KeplerianElements& mean,
                               const ZonalParameters&   zp    = {},
                               SecularOrder             order = SecularOrder::J2SquaredJ4);
    /// From an inertial state whose z axis is the Earth's pole; its osculating elements are taken
    /// as mean elements (see propagateSecular).
    explicit SecularPropagator(const StateVector&     s0,
                               const ZonalParameters& zp    = {},
                               SecularOrder           order = SecularOrder::J2SquaredJ4);
    /// From a 12-component state; `mode` is as for propagateSecular(const StateVector12&, ...).
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
    void propagateGrid12(double start, double step, std::span<StateVector12> out) const noexcept;
    [[nodiscard]] std::vector<StateVector>   propagateGrid(double start, double step, std::size_t count) const;
    [[nodiscard]] std::vector<StateVector12> propagateGrid12(double start, double step, std::size_t count) const;

    [[nodiscard]] const KeplerianElements& meanElements() const noexcept { return mean_; }
    [[nodiscard]] const SecularRates&      rates() const noexcept { return rates_; }
    /// Input acceleration minus the model's at the epoch (zero unless built from a StateVector12).
    [[nodiscard]] const Vec3& accelerationResidual() const noexcept { return da0_; }
    /// Input jerk minus the model's at the epoch (zero unless built from a StateVector12).
    [[nodiscard]] const Vec3& jerkResidual() const noexcept { return dj0_; }

private:
    /// Time-dependent angles of one state: reduced mean anomaly, eccentric anomaly and the sines
    /// and cosines needed to build the state.
    struct Angles {
        double M{0.0}, E{0.0}, sinE{0.0}, cosE{1.0};
        double sinRaan{0.0}, cosRaan{1.0}, sinArgp{0.0}, cosArgp{1.0};
    };

    void                        init();
    void                        fieldDerivatives(const Vec3& r, const Vec3& v, Vec3& a, Vec3& j) const noexcept;
    [[nodiscard]] Angles        angles(double dt) const noexcept;
    [[nodiscard]] StateVector   state(const Angles& x) const noexcept;
    [[nodiscard]] StateVector12 complete(double dt, StateVector s) const noexcept;
    template <class Emit>
    void grid(double start, double step, std::size_t count, Emit&& emit) const noexcept;

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
    Vec3                da0_;
    Vec3                dj0_;
};

/// Inclination (rad) of a sun-synchronous orbit with mean a and e (first-order J2;
/// nodal rate 360 deg per tropical year). Throws astro::Error if no solution exists.
[[nodiscard]] double sunSynchronousInclination(double a, double e = 0.0, const ZonalParameters& zp = {});

/// Mean semi-major axis of a repeat-ground-track orbit that completes `revs`
/// revolutions in `days` nodal days, for given e and i (J2, iterative).
[[nodiscard]] double
repeatGroundTrackSemiMajorAxis(int revs, int days, double e, double i, const ZonalParameters& zp = {});

/// Critical inclination (rad) where the J2 apsidal rotation vanishes: acos(1/sqrt(5)).
inline constexpr double kCriticalInclination = 1.1071487177940904;

} // namespace astro
