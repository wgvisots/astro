/**********************************************************************
 * impact.h
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
#include <optional>
#include <vector>

#include "astro/core/mat6.h"
#include "astro/core/vec.h"
#include "astro/frames/eop.h"
#include "astro/frames/rotation_series.h"
#include "astro/frames/transform.h"
#include "astro/geodesy/ellipsoid.h"
#include "astro/geodesy/geodetic.h"
#include "astro/orbit/central_body.h"
#include "astro/orbit/integrators.h"
#include "astro/orbit/numerical.h"
#include "astro/time/time.h"

/// \file
/// Ballistic impact prediction: where and when a trajectory meets the ellipsoid, and the impact
/// error ellipse from the state covariance.
namespace astro {

/// An error ellipse in the local horizontal plane.
struct ErrorEllipse {
    double semiMajor{0.0}; ///< semi-major axis (m)
    double semiMinor{0.0}; ///< semi-minor axis (m)
    /// Rotation of the semi-major axis, clockwise from north, in [0, pi) (rad); 0 for a circle.
    double rotation{0.0};

    /// The same ellipse scaled to contain `probability` (in (0, 1)) of a 2-D Gaussian, e.g. 0.95.
    /// Throws astro::Error for a probability outside (0, 1).
    [[nodiscard]] ErrorEllipse scaled(double probability) const;
};

/// The 1-sigma ellipse of an east/north covariance (m^2). A 1-sigma ellipse contains 39.3 % of a
/// 2-D Gaussian; use ErrorEllipse::scaled for other probabilities.
[[nodiscard]] ErrorEllipse errorEllipse(double eastVariance, double northVariance, double eastNorthCovariance) noexcept;

/// Factor k by which a 1-sigma ellipse is scaled to contain `probability` of a 2-D Gaussian:
/// k = sqrt(-2 ln(1 - p)), e.g. 1.1774 for 50 %, 2.4477 for 95 %. Throws astro::Error for a
/// probability outside (0, 1).
[[nodiscard]] double ellipseScale(double probability);

/// Uncertainty of an impact, from the covariance of the initial state.
struct ImpactUncertainty {
    /// Covariance of position and velocity at the nominal impact time, in the input frame
    /// (m^2, m^2/s, m^2/s^2).
    Mat6         covariance{};
    double       eastVariance{0.0};        ///< variance of the impact point, east (m^2)
    double       northVariance{0.0};       ///< variance of the impact point, north (m^2)
    double       eastNorthCovariance{0.0}; ///< east/north covariance of the impact point (m^2)
    ErrorEllipse ellipse;                  ///< 1-sigma impact error ellipse on the surface
    double       timeSigma{0.0};           ///< standard deviation of the impact time (s)
};

/// Where and when a trajectory meets the surface.
struct Impact {
    Time                             time;            ///< impact time (in the time scale of the initial epoch)
    double                           flightTime{0.0}; ///< seconds from the initial epoch to the impact
    StateVector                      state;        ///< position and velocity at impact, in the input frame of that time
    StateVector                      itrs;         ///< position and ground-relative velocity at impact, ITRS ("ECEF")
    Geodetic                         point;        ///< geodetic latitude, longitude, and height (= surfaceHeight)
    double                           speed{0.0};   ///< speed relative to the ground (m/s)
    double                           angle{0.0};   ///< impact angle: flight path below the local horizontal (rad)
    double                           heading{0.0}; ///< direction of travel over the ground, clockwise from north (rad)
    std::optional<ImpactUncertainty> uncertainty;  ///< present when a covariance was given
};

/// Settings for predictImpact.
struct ImpactOptions {
    /// Longest flight to search (s); no impact by then gives std::nullopt.
    double maxFlightTime{7200.0};
    /// Height of the impact surface above the ellipsoid (m): 0 is the ellipsoid itself; use the
    /// local terrain or geoid height for the real ground.
    double    surfaceHeight{0.0};
    Ellipsoid ellipsoid{Ellipsoid::wgs84()}; ///< reference ellipsoid of the surface
    /// Central body of the gravity model: point mass and every zonal term of body.J (by default
    /// J2, J3 and J4), aligned with the true pole at the epoch.
    CentralBody body{CentralBody::earth()};
    /// Further forces (e.g. drag), evaluated with t in seconds after the epoch and GCRS r, v.
    std::vector<std::shared_ptr<const ForceModel>> additionalForces;
    /// Integrator settings. The default (relTol 1e-9) keeps the impact point of a 40-minute
    /// ballistic flight within a few millimetres of a tight-tolerance solution.
    AdaptiveOptions integrator{[] {
        AdaptiveOptions o;
        o.relTol = 1e-9;
        o.absTol = 1e-6;
        return o;
    }()};
    /// White acceleration process noise added to the covariance (m^2/s^3); 0 for none.
    double accelerationNoisePsd{0.0};
    /// Impact time tolerance (s). At 7 km/s, 1e-7 s is 0.7 mm along the trajectory.
    double timeTolerance{1e-7};
    /// Precession-nutation model for GCRS -> ITRS. IAU 2006/2000B is within 1 mas of 2000A
    /// (3 cm at the surface) and about 20x cheaper to set up.
    FrameOptions frameOptions{PrecessionNutationModel::IAU2006_2000B, true};
    /// Precomputed GCRS <-> ITRS rotations covering [epoch, epoch + maxFlightTime], shared between
    /// calls for speed (many predictions in the same time window). nullptr: built per call. When
    /// set, the series is used as it is: the `eop` argument of predictImpact and `frameOptions`
    /// above play no part.
    const ItrsRotationSeries* rotations{nullptr};
};

/// Ballistic impact of a trajectory on the surface `options.surfaceHeight` above the ellipsoid.
///
/// `state` is position (m) and velocity (m/s) at `epoch` in the celestial frame `frame` (GCRS,
/// J2000/EME2000, MOD, TOD, TEME or CIRS; of-date frames are taken at their own time). It is
/// propagated in GCRS with the gravity model of `options`, and its height above the ellipsoid is
/// followed in the ITRS (from `eop`) until it first falls through the surface. Starting on the
/// surface (a launch) is fine: only a downward crossing counts. Returns std::nullopt if there is
/// no impact within options.maxFlightTime (e.g. an orbit).
///
/// The Earth orientation data matter: without UT1 - UTC the impact point can be wrong by up to
/// ~420 m east-west, and without polar motion by up to ~15 m (pass ConstantEop{} only if that is
/// acceptable). Throws astro::Error for an Earth-fixed `frame` or invalid options; if the Earth
/// orientation data in use (`eop`, or options.rotations when given) do not cover the window from
/// `epoch` to `epoch` + options.maxFlightTime; and as NumericalPropagator::propagateToEvent does
/// (the integration can fail).
[[nodiscard]] std::optional<Impact> predictImpact(const Time&          epoch,
                                                  const StateVector&   state,
                                                  Frame                frame,
                                                  const EopProvider&   eop,
                                                  const ImpactOptions& options = {});

/// As above, with the 6x6 covariance of `state` (position then velocity, in `frame`; m^2, m^2/s,
/// m^2/s^2). The covariance is propagated to the nominal impact time (with process noise from
/// options.accelerationNoisePsd) and mapped onto the surface: a position error moves the impact
/// along the ground-relative velocity until it reaches the surface again. The result is linear,
/// which holds while the errors are small next to the trajectory's curvature (see the accuracy
/// page for a Monte Carlo comparison). Throws astro::Error as the overload above does.
[[nodiscard]] std::optional<Impact> predictImpact(const Time&          epoch,
                                                  const StateVector&   state,
                                                  const Mat6&          covariance,
                                                  Frame                frame,
                                                  const EopProvider&   eop,
                                                  const ImpactOptions& options = {});

} // namespace astro
