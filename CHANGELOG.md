# Changelog

All notable changes to astro. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and the version numbers follow [semantic versioning](https://semver.org/spec/v2.0.0.html): before 1.0, a
change that breaks existing code raises the minor number.

## [0.3.0] - unreleased

### Breaking changes

- **A trajectory carries its frame.** "A state as a function of time" is one type, `Trajectory`,
  now in its own header `astro/orbit/trajectory.h` (with `TrajectoryModel` and
  `UncertainTrajectory`; `close_approach.h` includes it). Every constructor takes a `Frame`, with no
  default, and `frame()` returns it. The two-body and secular models are inertial and throw
  `astro::Error` for an Earth-fixed frame (TIRS, PEF, ITRS); a function may be in any frame.
  `UncertainTrajectory(epoch, state, covariance, frame, ...)` takes the frame after the covariance.
- **The search functions take a `Trajectory` and check its frame**, throwing `astro::Error` that
  names the frame given and the one needed: `findPasses` and `groundTrack` need `Frame::ITRS`;
  `findEclipses(sat, sun, ...)` needs one frame for both (any frame; it uses positions); and
  `closestApproach`, `closeApproaches`, `assessConjunction` and `assessConjunctions` need one inertial
  frame for both objects. The results for valid input are unchanged.
- **Removed:** `ItrsStateFunction`, `sgp4ItrsFunction`, and the `std::function<Vec3(const Time&)>`
  parameters of `groundTrack` and `findEclipses`. New helpers make the usual trajectories:
  `sgp4Trajectory(sgp4, eop, frame)` (SGP4's TEME output converted to `frame`; in `Frame::ITRS` it
  gives the same numbers as `sgp4ItrsFunction` did) and `bodyTrajectory(eph, target, center)` (an
  ephemeris state, in `Frame::GCRS`).

Before / after:

```cpp
// 0.2.0
findPasses(sgp4ItrsFunction(sgp4, eop), site, start, end, mask);
groundTrack([&](const Time& t) { return itrsPosition(t); }, start, end, 60.0);
findEclipses(satGcrsPosition, sunGcrsPosition, start, end);
closestApproach({epoch, state1}, {epoch, elements2}, start, end);
UncertainTrajectory a{epoch, state1, P1};
Trajectory numerical([&prop](const Time& t) { return prop.propagate(t); });

// 0.3.0
findPasses(sgp4Trajectory(sgp4, eop, Frame::ITRS), site, start, end, mask);
groundTrack({Frame::ITRS, [&](const Time& t) { return itrsState(t); }}, start, end, 60.0);
findEclipses(sgp4Trajectory(sgp4, eop, Frame::GCRS), bodyTrajectory(eph, Body::Sun), start, end);
closestApproach({epoch, state1, Frame::GCRS}, {epoch, elements2, Frame::GCRS}, start, end);
UncertainTrajectory a{epoch, state1, P1, Frame::GCRS};
Trajectory numerical(Frame::GCRS, [&prop](const Time& t) { return prop.propagate(t); });
```

## [0.2.0] - unreleased

The first tagged release. It follows a top-to-bottom review of the library (November 2026), which
fixed its defects, made its interface consistent, and strengthened its tests and documentation.
The pull request of each item is given in brackets.

### Breaking changes

Every change below makes old code fail to compile, **except the first**, which compiles and changes
meaning.

- **`Sgp4::propagate(double)` takes seconds since the element set's epoch, not minutes.** The minute
  form of the reference code is `Sgp4::propagateMinutes`. The non-throwing overload changed the
  same way. [#13]
- **One type for a central body, `CentralBody`** (`astro/orbit/central_body.h`), with the named values
  `CentralBody::earth()`, `earthJ2J4()` and `pointMass(mu)`. It replaces `ZonalParameters`, the
  `(J, mu, re)` arguments of `GravityField` and of the functions in `gravity.h`, and the `J`, `mu`,
  `re` members of `TrackerOptions` and `ImpactOptions` (now `body`). `NumericalPropagator::zonal`
  replaces the `j2j4` overloads that took loose constants; `GravityField::mu()` is `body().mu`. [#12]
- **Propagator argument order:** `propagateKeplerWithCovariance(s0, dt, P0, mu)` (the covariance
  after the time). [#13]
- **Names** follow one set of rules, listed in the conventions guide [#14]:
  - conversions are `xToY`: `deg2rad` is `degToRad` (and likewise `rad2deg`, `arcsec2rad`,
    `mas2rad`, ...); `raDecFromVector` is `vectorToRaDec`, `vectorFromRaDec` is `raDecToVector`,
    `cipXYFromMatrix` is `matrixToCipXY`, `semiMajorAxisFromMeanMotion` is
    `meanMotionToSemiMajorAxis`, `deputyFromRelativeRic` is `relativeRicToDeputy`,
    `topocentricFromGcrs` is `gcrsToTopocentric`, `monthDayFromDayOfYear` is `dayOfYearToMonthDay`,
    `dateFromJulianDayNumber` is `julianDayNumberToDate`, `timeScaleFromString` is `parseTimeScale`;
  - public constants are in `astro::constants`: `constants::earthRotationRateEra`,
    `criticalInclination`, `sunriseElevation`; `earthRadiusIers2010` is removed (use `iers2010A`);
  - accessors are nouns: `JulianDate::calendar()`, `Time::sysTime()`;
  - data before time: `ItrsRotationSeries::gcrsToItrs(r, dt)` and `itrsToGcrs(s, dt)`;
  - no generic free names: `HelmertParams::inverse()`, `trackerPropagator`, `earthShadow`,
    `twilightEvents`, `coverageFootprint`;
  - `JulianDate(double)` and the `StateVector12` constructors of the secular propagators are
    `explicit`.
- **`Mat6` is a struct** used like `Mat3`: `P(i, j)`, `a * b`, `a + b`, `a.transpose()`,
  `Mat6::identity()`, `P.block(i, j)`. The free `identity6`, `transpose`, `multiply`, `add` and
  `block` are removed. [#15]
- **Behaviour changes that throw where a wrong answer was returned** [#8, #10]:
  - equinoctial elements are defined for ellipses only (hyperbolic input throws, or gives NaN from
    `keplerianToEquinoctial`);
  - the search functions (`findPasses`, `findEclipses`, `groundTrack`, `findElevationEvents`) throw
    for a reversed window or a bad step;
  - the TLE, leap-second, Earth-orientation and angle parsers are strict: a damaged field throws,
    naming the line and the field;
  - `Time::iso8601` covers the years 0000 to 9999 only;
  - the integrator throws for non-finite times, unusable options and non-finite derivatives;
  - `EstimateTracker` validates its options when it is built.
- **Numbers that change because they were wrong** [#8]:
  - `helmertTransform` with rates takes and returns velocity in m/s (it treated it as m/yr);
  - `hmsToRad` / `dmsToRad` honour a negative first argument;
  - `brentMinimize`'s `tol` is absolute, with Brent's sqrt(epsilon) floor;
  - `refractionBennett` and `refractionSaemundsson` are zero at the zenith and finite everywhere;
  - TLE epochs on a leap-second day move by up to 1 s.

### Added

- Closest approach and conjunction lists for any two trajectories, collision probability (Foster and
  Estes; checked against NASA CARA's test cases), and `assessConjunction`, which combines them with
  propagated covariances. [#5]
- Covariance through the two-body (closed-form state transition matrix) and secular propagators;
  `SecularPropagatorCorrected`. [#5]
- `NumericalPropagator::propagateStm` and `propagateWithCovariance` at a single time. [#13]
- A Lambert solver written from Izzo (2015), with an accuracy study. [#7]
- Ballistic impact prediction with an impact error ellipse; event detection in the numerical
  propagator. [#3]
- `astro/version.h`: `ASTRO_VERSION_MAJOR`, `_MINOR`, `_PATCH`, `ASTRO_VERSION_STRING`,
  `ASTRO_VERSION` and `astro::versionString`.
- Debian and RPM packages with CPack (`cpack -G DEB`, `cpack -G RPM`). [#11]
- Fuzz targets for every parser (`fuzz/`, option `ASTRO_BUILD_FUZZERS`), and a line-coverage build
  (option `ASTRO_COVERAGE`, report with `tools/coverage_report.py`). [#10]

### Fixed

About 60 defects, each with a regression test that failed first [#6, #8, #10]. Among them:

- `propagateKepler` threw "did not converge" on valid orbits near e = 1 and for short steps;
- numerical grids that span the epoch were 8 m off at the epoch and had negative process-noise
  variances;
- ISO-8601 UTC offsets were 1 s wrong across a leap second; UT1 was 3 ns off on the two shortened
  UTC days of 1961 and 1968;
- TLE blank-padded fields were read 10 times too large, and small B* values were written wrongly;
- `findPasses` reported a pass cut by the window at the mask elevation; transits in the first and
  last sample interval were lost;
- Euler angles lost 10 digits near gimbal lock; `slerp` was off by 1e-6 rad; `stateToElements` had
  a 1.5e-8 rad floor; the hyperbolic Kepler solver threw for near-parabolic orbits;
- the SPK reader returned data from segments in other frames, and NaN from damaged records;
- undefined behaviour in the secular state transition matrix that one compiler setting turned into a
  wrong matrix.

### Documentation

- Citations checked against their sources; wrong figures corrected (for example the angle between
  true-of-date frames and GCRS, 0.27-0.42 degrees); every public function says what it throws;
  references added. [#9]
- Third-party notices for the ERFA routines and the SGP4 reference code that parts of the library
  follow. [#7]

### Infrastructure

- CI jobs for AddressSanitizer with UBSan, fuzzing of the parsers, line coverage (98 %), and the
  packages; 538 tests. [#10, #11]

## [0.1.0]

The first milestone [#1]: time scales, geodesy, reference frames, ephemerides, astrometry, orbit
propagation and satellite geometry. Not tagged.

[0.2.0]: https://github.com/wgvisots/astro/compare/db6a4b2...v0.2.0
[0.1.0]: https://github.com/wgvisots/astro/pull/1
