# Public API inventory for 1.0

| | |
|---|---|
| Status | Decided by the owner, and carried out in 0.4.0 |
| Date | 2026-11-15, from `main` at 015bdd1 |
| Scope | every header under `include/astro/`: 48 headers, 385 top-level declarations, 501 exported symbols |

This is the review of the public API before 1.0: each header, the recommendation and its reasons,
and the owner's decision in the last column (an empty cell agrees with the recommendation). The
[CHANGELOG](https://github.com/wgvisots/astro/blob/main/CHANGELOG.md) lists the resulting changes, with what to use instead.

## Outcome

- Every recommendation was accepted, except one: `tleChecksum` stays public, for services that
  check a two-line element set before passing it on.
- The two headers left to the owner, `precession_nutation.h` and `rotation.h`, stay public, as do
  `decodeAlpha5` and `encodeAlpha5`.
- Changes C1 to C4 were made.
- Afterwards: 47 public headers (`rootfind.h` and `integrators.h` became internal, and
  `orbit/integration.h` is new), and 492 exported symbols.

## The rule I used

From 1.0, the [compatibility policy](../compatibility.md) promises source compatibility. So the two
directions cost different amounts:

| Change after 1.0 | Cost |
|---|---|
| Make an internal helper public | free: an addition, allowed in any 1.x release |
| Make a public function internal, or change its signature | a major version (2.0) |

So **when in doubt, internal**. Public means: it is part of what astro is for, a user would
reasonably call it, and we are willing to keep its signature for all of 1.x.

"Internal" means the declaration moves to a header under `src/`, which is not installed. The
library and the tests can still use it (the tests get `src/` on their include path); users cannot.

## Summary

| Recommendation | Headers | Notes |
|---|---|---|
| Keep public, unchanged | 38 | the bulk of the library |
| Keep public, with a change | 3 | `transform.h` (C2), `numerical.h` (C1), `integrators.h` (C3: two structs stay public in a new header, the rest becomes internal) |
| Split: part public, part internal | 4 | `two_body.h`, `gravity.h`, `estimate_tracker.h`, `tle.h` |
| Internal | 1 | `rootfind.h` |
| Needs your call | 2 | `precession_nutation.h`, `rotation.h`. My recommendation is keep, for both. |

**Must stay public:** `propagateStm`, `propagateGridStm`, `propagateKeplerStm` and `StateTransition`,
because orbit determination builds on them ([design](orbit-determination.md)).

## Changes the review should make, beyond moving headers

These are the items that would be breaking after 1.0, so they need deciding now.

| # | Item | Recommendation | Opinion | Your decision |
|---|---|---|---|---|
| C1 | `ForceModel::partials(t, r, v, Mat3& dadr, Mat3& dadv)` returns two of its results through out-parameters | Return a struct, `ForcePartials {Vec3 acceleration; Mat3 dadr; Mat3 dadv;}` | The error policy (rule 6) says no out-parameters, and this is the one public case left. It is a virtual function that users override, so its signature is frozen at 1.0. Orbit determination will add force models that implement it, so fixing it now also saves churn later. |Return a struct | 
| C2 | `eciToEcef` and `ecefToEci` (4 inline aliases of `gcrsToItrs` and `itrsToGcrs`) | Remove | "ECI" is ambiguous (GCRS? J2000? TEME? true of date?), which is exactly the confusion the frame-named functions avoid, and the `xToY` naming rule names frames precisely. They are used by no example or guide page. Keep the words "ECI" and "ECEF" in the docs, pointing to the GCRS and ITRS functions. |Agreed, stick to gcrsToItrs, ... | 
| C3 | `AdaptiveOptions` and `AdaptiveStats` live in `integrators.h`, which public headers include only for these two | Move them to a new public header `astro/orbit/integration.h`, and make the rest of `integrators.h` internal | Users need the two structs: they set the tolerances of `NumericalPropagator` and read its statistics. They do not need the integrator templates. |Agreed, please move | 
| C4 | `j2Acceleration`, documented as "for reference/testing" | Move into the tests | Its own documentation says it is not for users. |Agreed, move to tests | 

## Inventory by header

Columns:
- **Used in examples or guides:** whether any example program or guide page uses it. If so, we
  have already shown it to users.
- **Recommendation:** Keep, Internal, Split, Change, or Decide (your call; my lean is in the opinion).

### Core

| Header | Items | Used in examples or guides | Recommendation | Opinion | Your decision |
|---|---|---|---|---|---|
| `core/vec.h` | `Vec3`, `StateVector`, `StateVector12`, `Mat3`, their operators, `dot`, `cross`, `norm`, `norm2`, `normalized`, `angleBetween`, `outer` | yes | Keep | The basic types of every interface. |keep |
| `core/mat6.h` | `Mat6`, operators, `congruence`, `rotateCovariance`, `mahalanobis2` | yes | Keep | Covariance work needs all of them, and the guide shows them. |keep |
| `core/angle.h` | unit conversions, `wrapPi`, `wrapTwoPi`, `wrap360`, `wrap180`, sexagesimal parse, format and decompose | yes | Keep | Small, stable, and wanted by everyone. |keep |
| `core/constants.h` | 52 constants: pi, unit factors, time constants, Earth and GM values | yes | Keep | Each one documents its source. |keep |
| `core/error.h` | `Error` | yes | Keep | The one exception type, from the error policy. |keep |
| `core/rotation.h` | `rotationAboutAxis`, `axisAngleToMatrix`, `Quat`, `slerp`, `EulerSequence`, `EulerAngles`, `eulerToMatrix`, `matrixToEuler`, `orthonormalize` | 2 uses in the guide | **Decide** (lean: keep) | General attitude mathematics, used very little inside astro. Keeping it is cheap: the conventions are written down and tested, and people working with frames often need exactly these conversions. The case for internal is that it is not astrodynamics. |keep |
| `core/rootfind.h` | `bisect`, `brentRoot`, `brentMinimize`, `MinimizeResult`, `newtonRoot` | no | **Internal** | Generic numerical tools that astro uses internally. Users have better general libraries for this, and keeping these public would freeze their template signatures for no gain. |internal |

### Time

| Header | Items | Used in examples or guides | Recommendation | Opinion | Your decision |
|---|---|---|---|---|---|
| `time/time.h` | `Time`, `TimeScale`, `toString`, `parseTimeScale`, `GnssSystem`, `GnssWeekTow`, `tdbMinusTt`, `resolveWeekRollover` | yes | Keep | The centre of the time module. `tdbMinusTt` is the documented Fairhead-Bretagnon series, which some users call directly. |keep |
| `time/julian_date.h` | `JulianDate`, `CalendarDate`, `YearMonthDay`, `isLeapYear`, `daysInMonth`, `dayOfYear`, `dayOfYearToDate`, `dayOfWeek`, `julianDayNumber`, `julianDayNumberToDate` | yes | Keep | Calendar arithmetic that users of a time library expect to find. |keep |
| `time/leap_seconds.h` | `LeapSecondEntry`, `LeapSecondTable`, `leapSecondTable`, `setLeapSecondTable` | yes | Keep | Needed for loading newer tables (the data-files guide). |keep |

### Reference frames

| Header | Items | Used in examples or guides | Recommendation | Opinion | Your decision |
|---|---|---|---|---|---|
| `frames/transform.h` | `Frame`, `FrameTransformer`, `FrameOptions`, `isEarthFixed`, `gcrsToCelestialMatrix`, `poleOfDateMatrix`, the one-off GCRS, ITRS and TEME conversions | yes | **Change** (C2) | Keep all of it except the four `eciToEcef` and `ecefToEci` aliases. |agreed |
| `frames/eop.h` | `EopValues`, `EopProvider`, `ConstantEop`, `EopRecord`, `EopInterpolation`, `EopTable` | yes | Keep | Users load and pass these everywhere. |keep |
| `frames/earth_rotation.h` | `earthRotationAngle`, `gmst06`, `gast06`, `gmst82`, `equationOfOrigins`, `polarMotionMatrix` | no | Keep | Sidereal time and the Earth rotation angle are standard quantities that people ask for by name. They implement fixed IAU definitions, so their signatures will not need to change. | |
| `frames/precession_nutation.h` | 20 functions and 4 types: `fundamentalArguments`, `nutationIau2000A`, `nutationIau2000B`, `nutation`, `meanObliquity06`, `fukushimaWilliams06`, `fwToMatrix`, `frameBiasMatrix06`, `biasPrecessionMatrix06`, `precessionMatrix06`, `nutationMatrix`, `npbMatrix06`, `cipXY06`, `matrixToCipXY`, `cipXY`, `cioLocatorS06`, `tioLocatorSp`, `celestialToIntermediateMatrix`, `poleOffsetsToNutation`, `PrecessionNutationModel` | `nutation` only | **Decide** (lean: keep) | These are the building blocks that SOFA and ERFA expose, each validated against ERFA. For keep: they implement fixed IAU standards, so their signatures will not change, and astronomy users value them. For internal: `FrameTransformer` covers the usual needs, and 20 functions is a large promise for a few users. A middle way: keep the matrices, `nutation`, `cipXY` and the CIO and TIO locators public, and make the intermediate pieces internal (`fundamentalArguments`, `fukushimaWilliams06`, `fwToMatrix`, `cipXY06`, `matrixToCipXY`). |keep |
| `frames/rotation_series.h` | `ItrsRotationSeries` | yes | Keep | The documented fast path for many epochs. |keep |

### Geodesy

| Header | Items | Used in examples or guides | Recommendation | Opinion | Your decision |
|---|---|---|---|---|---|
| `geodesy/ellipsoid.h` | `Ellipsoid` | yes | Keep | |keep |
| `geodesy/geodetic.h` | `Geodetic`, `geodeticToEcef`, `ecefToGeodetic`, the latitude conversions, `Spherical`, `cartesianToSpherical`, `sphericalToCartesian` | yes | Keep | |keep |
| `geodesy/geodesic.h` | `geodesicInverse`, `geodesicDirect` and their result types, `haversineDistance` | yes | Keep | |keep |
| `geodesy/helmert.h` | `HelmertParams`, `HelmertTimeDependent`, `HelmertConvention`, two `helmertTransform` | yes | Keep | |keep |
| `geodesy/local.h` | ENU and NED matrices and conversions, `LookAngles`, `enuToLookAngles`, `lookAnglesToEnu`, `ecefToLookAngles`, `lookAnglesToEcef` | yes | Keep | Settled in 0.3.0. |keep |

### Ephemerides and astrometry

| Header | Items | Used in examples or guides | Recommendation | Opinion | Your decision |
|---|---|---|---|---|---|
| `ephemeris/ephemeris.h` | `Body`, `toString`, `naifId`, `Ephemeris` | yes | Keep | |keep |
| `ephemeris/low_precision.h` | `LowPrecisionEphemeris`, `sunPositionLowPrecision`, `moonPositionLowPrecision` | yes (the class) | Keep | The two free functions repeat what the class gives. They are harmless and convenient for a quick Sun or Moon position, so I would keep them, but they could go if you prefer one way to do each thing. |keep |
| `ephemeris/spk.h` | `SpkEphemeris`, `SpkSegmentInfo`, `ephemerisTime` | yes | Keep | |keep |
| `astrometry/apparent.h` | `CatalogStar`, `ObserverEphemeris`, `observerEphemeris`, `lightDeflection`, `aberration`, star and body places, `gcrsToTopocentric`, sidereal times, `hourAngle` | yes | Keep | `lightDeflection` and `aberration` are the ERFA building blocks; orbit determination's line-of-sight model will use `aberration`. |keep |
| `astrometry/refraction.h` | `Atmosphere`, `RefractionConstants`, `refractionConstants`, `refractTrueToObserved`, `refractObservedToTrue`, `refractionBennett`, `refractionSaemundsson` | yes | Keep | |keep |
| `astrometry/rise_set.h` | `findElevationEvents`, rise, set and twilight, `moonPhase`, and their types | yes | Keep | |keep |
| `astrometry/spherical.h` | `RaDec`, conversions, separations, ecliptic and galactic rotations | yes | Keep | |keep |

### Orbits and propagation

| Header | Items | Used in examples or guides | Recommendation | Opinion | Your decision |
|---|---|---|---|---|---|
| `orbit/state_covariance.h` | `StateCovariance`, `StateTransition` | yes | **Keep** (required) | `StateTransition` must stay public for orbit determination. |keep |
| `orbit/state12.h` | `HigherOrderMode`, `Propagation12` | yes | Keep | The 12-component state is a deliberate feature, with public overloads in three propagators. |keep |
| `orbit/central_body.h` | `CentralBody` | yes | Keep | |keep |
| `orbit/elements.h` | Keplerian and equinoctial elements and conversions, `stateToEquinoctialJacobian`, period, mean motion, energy, flight-path angle, apsides | yes | Keep | The equinoctial Jacobian is what a user needs to convert a covariance into elements. |keep |
| `orbit/kepler.h` | 14 anomaly conversions (elliptic, hyperbolic, parabolic, any) | no | Keep | Standard, small, and stable. Users ask for them by name. |keep |
| `orbit/two_body.h` | `propagateKepler` (6 and 12), `propagateKeplerStm`, `propagateKeplerWithCovariance`, `LambertSolution`, `lambertIzzo`, `lambert` | yes | Keep | `propagateKeplerStm` must stay public for orbit determination. |keep |
| `orbit/two_body.h` | `Stumpff`, `stumpff` | no | **Split: internal** | An implementation detail of the universal-variable propagator. |internal |
| `orbit/secular.h` | `SecularOrder`, `SecularRates`, `secularRates`, `propagateSecular` (3), `SecularPropagator`, `SecularPropagatorCorrected`, `sunSynchronousInclination`, `repeatGroundTrackSemiMajorAxis` | yes | Keep | |keep |
| `orbit/numerical.h` | `ForceModel`, `GravityField`, `ThirdBodyForce`, `UnmodeledAcceleration`, `Crossing`, `EventFunction`, `PropagationEvent`, `NumericalPropagator` | yes | **Change** (C1) | Keep all of it. Change `ForceModel::partials` to return a struct. `propagateStm` and `propagateGridStm` must stay public. |Agreed, return a struct |
| `orbit/integrators.h` | `AdaptiveOptions`, `AdaptiveStats` | yes | **Change** (C3) | Move to a new public header, `orbit/integration.h`. |Agreed, move |
| `orbit/integrators.h` | `OdeState`, `axpy`, `rk4Step`, `integrateRk4`, `integrateRkf78`, `integrateRkf78Steps` | no | **Internal** | The integrator templates are machinery. Users get integration through `NumericalPropagator`. |internal |
| `orbit/gravity.h` | `pointMassAcceleration`, `zonalAcceleration`, `thirdBodyAcceleration` | no | **Split: keep** | Useful to anyone writing a `ForceModel` (the flyout tool, for one). Small and stable. |keep |
| `orbit/gravity.h` | `pointMassJerk`, `zonalJerk`, `zonalAccelerationAndJerk`, `pointMassGradient`, `zonalGradient` | no | **Split: internal** | Machinery for the 12-component states and the state transition matrix. `zonalAccelerationAndJerk` also uses out-parameters. They can be made public later if a user needs them. |internal |
| `orbit/gravity.h` | `j2Acceleration` | no | **Internal** (C4) | Into the tests. |internal |
| `orbit/sgp4.h` | `Sgp4`, `Sgp4Error`, `Sgp4Result`, `Sgp4Status`, `Sgp4Gravity`, `Sgp4OpsMode` | yes | Keep | Settled in 0.3.0. |keep |
| `orbit/tle.h` | `Tle` | yes | Keep | |keep |
| `orbit/tle.h` | `tleChecksum` | no | **Split: internal** | `Tle::parse` checks the checksum and `Tle::line1` and `line2` write it, so users never need it. |keep - useful for services that use it to validate data before giving it to us |
| `orbit/tle.h` | `decodeAlpha5`, `encodeAlpha5` | no | **Decide** (lean: keep) | Users who handle catalogue numbers above 99999 from other sources (conjunction messages, catalogue files) need exactly this conversion, and it is a fixed format. |keep |
| `orbit/trajectory.h` | `Trajectory`, `UncertainTrajectory`, `TrajectoryModel`, `sgp4Trajectory`, `bodyTrajectory` | yes | Keep | Settled in 0.3.0. |keep |
| `orbit/close_approach.h` | closest approach, conjunction assessment, collision probability, `probabilityInCircle` | yes | Keep | |keep |
| `orbit/impact.h` | `predictImpact`, `Impact`, `ImpactOptions`, `ImpactUncertainty`, `ErrorEllipse`, `errorEllipse`, `ellipseScale` | yes | Keep | |keep |

### Satellite geometry and tracking

| Header | Items | Used in examples or guides | Recommendation | Opinion | Your decision |
|---|---|---|---|---|---|
| `satellite/look.h` | `findPasses`, `Pass`, `dopplerShift`, `dopplerShiftRelativistic` | yes | Keep | |keep |
| `satellite/ground.h` | sub-satellite point, ground track, coverage, `lineOfSight` | yes | Keep | |keep |
| `satellite/eclipse.h` | `earthShadow`, `findEclipses`, `betaAngle`, and their types | yes | Keep | Orbit determination's radiation-pressure model will use `earthShadow`. |keep |
| `satellite/orbital_frames.h` | RIC, LVLH, VNC and NTW matrices, relative states | yes | Keep | |keep |
| `satellite/estimate_tracker.h` | `EstimateTracker`, `StateEstimate`, `TrackerOptions`, `ConsistencyTest`, `ConsistencyCheck`, `checkConsistency`, `trackerPropagator` | yes | Keep | |keep |
| `satellite/estimate_tracker.h` | `chiSquareSurvival3`, `chiSquareSurvival6` | no | **Split: internal** | `ConsistencyCheck` already reports the p-value, so users do not need these. A general chi-square function belongs in a statistics library. |internal |

### Umbrella

| Header | Items | Used in examples or guides | Recommendation | Opinion | Your decision |
|---|---|---|---|---|---|
| `astro.h` | includes every public header | yes | Keep | Updated to match whatever is decided. |keep |

## What happens after your decisions

One pull request, released as 0.4.0:

1. Move each internal declaration to a header under `src/`, and point the library and the tests at it.
2. Make changes C1 to C4.
3. Update `astro.h`, the docs and the export check. The number of exported symbols goes down; I'll
   report the new count.
4. A CHANGELOG section under "Breaking changes" listing every removal, with what to use instead.
