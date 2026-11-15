# Changelog

All notable changes to astro. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and the version numbers follow [semantic versioning](https://semver.org/spec/v2.0.0.html): before 1.0, a
change that breaks existing code raises the minor number.

## [Unreleased]

### Breaking changes

The review of the public API before 1.0 ([inventory and decisions](docs/design/api-review-1.0.md)).
From 1.0, making something public is an addition, but making it internal needs 2.0. So what users
have no need for became internal now, and the interfaces that 1.0 will freeze were fixed first.

- **`ForceModel::partials` returns a struct.** `partials(t, r, v)` returns
  `ForcePartials {acceleration, dadr, dadv}`, instead of returning the acceleration and filling two
  `Mat3&` out-parameters. A force model that overrides it changes its signature:

```cpp
// 0.3
Vec3 partials(double t, const Vec3& r, const Vec3& v, Mat3& dadr, Mat3& dadv) const override;
// 0.4
ForcePartials partials(double t, const Vec3& r, const Vec3& v) const override;
```

- **Removed `eciToEcef` and `ecefToEci`.** Use `gcrsToItrs` and `itrsToGcrs`, which they called.
  "ECI" does not say which frame is meant.
- **`AdaptiveOptions` and `AdaptiveStats` moved to the new header `astro/orbit/integration.h`.**
  `numerical.h`, `impact.h` and `estimate_tracker.h` include it, so code that reaches them through
  those headers is unaffected.
- **Now internal**, in headers under `src/`, which are not installed:

| No longer public | Use instead |
|---|---|
| `astro/core/rootfind.h`: `bisect`, `brentRoot`, `brentMinimize`, `MinimizeResult`, `newtonRoot` | Boost.Math (header-only): `boost::math::tools::bisect`, `toms748_solve`, `brent_find_minima` |
| `astro/orbit/integrators.h`: `OdeState`, `rk4Step`, `integrateRk4`, `integrateRkf78`, `integrateRkf78Steps` | `NumericalPropagator`; for other equations, Boost.Odeint |
| `stumpff`, `Stumpff` | nothing: they were a detail of `propagateKepler` |
| `pointMassJerk`, `zonalJerk`, `pointMassGradient`, `zonalGradient`; `zonalAccelerationAndJerk` (removed) | `GravityField::jerk` and `GravityField::partials` (point mass and zonal terms together) |
| `j2Acceleration`, now a reference in the tests | `zonalAcceleration` with a `CentralBody` whose `J` holds J2 only |
| `chiSquareSurvival3`, `chiSquareSurvival6` | the p-values in `ConsistencyCheck` |

- The shared library exports 492 symbols, 9 fewer.

### Platforms

- **Windows with MSVC is supported** (docs/compatibility.md): Visual Studio 2022 and newer, static
  library and DLL, and the Conan package. CI builds and tests the static library with Visual Studio
  2022 and 2026 and the DLL with 2026, warnings as errors, and runs `conan create` for both library
  types. The getting-started guide explains the Windows build.
- **Fixed:** 25 functions were declared `ASTRO_EXPORT` twice (0.3.0's export change), which MSVC
  rejects; GCC and Clang accepted it.
- The generated `astro/export.h` also defines `ASTRO_DLL_WARNINGS_PUSH` and `ASTRO_DLL_WARNINGS_POP`,
  which silence MSVC's warnings about exported classes with standard-library members (C4251, C4275)
  in astro's own headers only.
- Tests: references computed in `long double` are skipped where it is no wider than `double`
  (MSVC); the Lambert closure test draws its cases in a way that is the same on every standard
  library.

### Build and packaging

- **A release tag publishes the GitHub release.** When CI passes on a tag `vX.Y.Z`, its new last
  job attaches the `.deb` and `.rpm` that run built, with `SHA256SUMS`, and takes the notes from
  the version's CHANGELOG section (`tools/release_notes.py`). It refuses a tag that does not match
  `project(VERSION)` and an undated section. On every other push it shows the notes it would use,
  without publishing. The steps are in CONTRIBUTING.

### Documentation

- **[Public API review](docs/design/api-review-1.0.md):** every public header, whether it stays
  public, and why.
- **[Orbit determination design](docs/design/orbit-determination.md):** the agreed plan for fitting
  orbits to position and velocity reports and to lines of sight, after 1.0: inputs, force models
  needed, methods, validation, order of work and open questions.

## [0.3.0] - 2026-11-14

### Breaking changes

- **A trajectory carries its frame.** "A state as a function of time" is one type, `Trajectory`,
  now in its own header `astro/orbit/trajectory.h` (with `TrajectoryModel` and
  `UncertainTrajectory`; `close_approach.h` includes it). Every constructor takes a `Frame`, with no
  default, and `frame()` returns it. The two-body and secular models are inertial and throw
  `astro::Error` for an Earth-fixed frame (TIRS, PEF, ITRS); a function may be in any frame.
  `UncertainTrajectory(epoch, state, covariance, frame, ...)` takes the frame after the covariance.
  [#17]
- **The search functions take a `Trajectory` and check its frame**, throwing `astro::Error` that
  names the frame given and the one needed: `findPasses` and `groundTrack` need `Frame::ITRS`;
  `findEclipses(sat, sun, ...)` needs one frame for both (any frame; it uses positions); and
  `closestApproach`, `closeApproaches`, `assessConjunction` and `assessConjunctions` need one inertial
  frame for both objects. The results for valid input are unchanged. [#17]
- **Removed:** `ItrsStateFunction`, `sgp4ItrsFunction`, and the `std::function<Vec3(const Time&)>`
  parameters of `groundTrack` and `findEclipses`. New helpers make the usual trajectories:
  `sgp4Trajectory(sgp4, eop, frame)` (SGP4's TEME output converted to `frame`; in `Frame::ITRS` it
  gives the same numbers as `sgp4ItrsFunction` did) and `bodyTrajectory(eph, target)` (an
  ephemeris state, in `Frame::GCRS`). [#17]
- **One look-angle type, `LookAngles`** (`astro/geodesy/local.h`), with `az`, `el`, `range` and
  `rangeRate`. It replaces `Aer` and moves from `astro/satellite/look.h` (which includes it).
  `rangeRate` is NaN unless the function was given a velocity, so a missing value never reads as
  "not moving". The functions follow the `xToY` rule: `enuToAer` is `enuToLookAngles`, `aerToEnu`
  is `lookAnglesToEnu`, `ecefToAer` is `ecefToLookAngles`, `aerToEcef` is `lookAnglesToEcef`, and
  `lookAngles(satItrs, site, ell)` is the `StateVector` overload `ecefToLookAngles(satItrs, site,
  ell)`. `starObserved`, `bodyObserved` and `gcrsToTopocentric` return `LookAngles` with
  `rangeRate` NaN. Every computed number is unchanged. [#18]
- **One error style**, written down at the top of `astro/core/error.h`: a failure throws
  `astro::Error`; a message starts with the function's name (`"Function: what is wrong (got
  <value>)"`; the messages of the time scales, the leap-second table, `gcrsToCelestialMatrix`,
  `poleOfDateMatrix`, `LowPrecisionEphemeris`, `parseHms`, `parseDms` and SGP4 that did not are
  reworded, so code that matched their text must be checked); a `noexcept`
  function is a formula that does not check its input; no out-parameters. So: `stumpff(psi)` returns
  a `Stumpff {c2, c3}`; `julianDayNumberToDate(jdn)` returns a `YearMonthDay {year, month, day}`, and
  `dayOfYearToMonthDay(year, doy, month, day)` is `dayOfYearToDate(year, doy)`, returning one too.
  SGP4's non-throwing calls are `Sgp4::tryPropagate(dt)` and `tryPropagateMinutes(minutes)`, which
  return an `Sgp4Result {state, status}`; the status is the new `enum class Sgp4Status`, whose
  values are the reference code's error numbers (`Decayed` is 6, `TimeNotFinite` 7), and
  `Sgp4Error::code()` returns it too. `keplerianToEquinoctial` throws `astro::Error` for an orbit
  that is not an ellipse instead of returning NaN, and is no longer `noexcept`; `mahalanobis2`
  throws for `n` outside 1 to 6 instead of returning `nullopt`. `toString(Frame)`,
  `toString(TimeScale)` and `toString(Body)` give `"INVALID_VALUE"` (as `operator<<` does) instead
  of `"?"` for a value that is not an enumerator. Every computed number is unchanged. [#19]

Before / after:

```cpp
// 0.2.0
findPasses(sgp4ItrsFunction(sgp4, eop), site, start, end, mask);
groundTrack([&](const Time& t) { return itrsPosition(t); }, start, end, 60.0);
findEclipses(satGcrsPosition, sunGcrsPosition, start, end);
closestApproach({epoch, state1}, {epoch, elements2}, start, end);
UncertainTrajectory a{epoch, state1, P1};
Trajectory numerical([&prop](const Time& t) { return prop.propagate(t); });
Aer        aer  = ecefToAer(target, site);
LookAngles look = lookAngles(satItrs, site);
double c2, c3;
stumpff(psi, c2, c3);
int y, m, d;
julianDayNumberToDate(jdn, y, m, d);
dayOfYearToMonthDay(2024, 60, m, d);
StateVector s;
int         code = sgp4.propagate(dt, s);    // 0 on success
if (code == 6) { /* decayed */ }
EquinoctialElements q = keplerianToEquinoctial(hyperbola); // every element NaN

// 0.3.0
findPasses(sgp4Trajectory(sgp4, eop, Frame::ITRS), site, start, end, mask);
groundTrack({Frame::ITRS, [&](const Time& t) { return itrsState(t); }}, start, end, 60.0);
findEclipses(sgp4Trajectory(sgp4, eop, Frame::GCRS), bodyTrajectory(eph, Body::Sun), start, end);
closestApproach({epoch, state1, Frame::GCRS}, {epoch, elements2, Frame::GCRS}, start, end);
UncertainTrajectory a{epoch, state1, P1, Frame::GCRS};
Trajectory numerical(Frame::GCRS, [&prop](const Time& t) { return prop.propagate(t); });
LookAngles aer  = ecefToLookAngles(target, site); // aer.rangeRate is NaN
LookAngles look = ecefToLookAngles(satItrs, site);
const auto [c2, c3] = stumpff(psi);
const YearMonthDay date = julianDayNumberToDate(jdn);
const YearMonthDay leap = dayOfYearToDate(2024, 60); // {2024, 2, 29}
const Sgp4Result r = sgp4.tryPropagate(dt);          // r.state, r.status
if (r.status == Sgp4Status::Decayed) { /* decayed */ }
keplerianToEquinoctial(hyperbola);                   // throws astro::Error
```

### Added

- **Conan 2 recipe** (`conanfile.py`, with `test_package/`), for static and shared builds. It gets
  Boost's headers from ConanCenter. CI builds and tests both kinds of package. [#21]
- Example `data_files`: loads whichever data files are present, falls back to what is built in, and
  reports how long the loaded data remains valid. [#21]

### Documentation

- **[Compatibility](docs/compatibility.md)**, the official policy [#21]:
  - what a version number promises: 0.x is in flux; 1.x keeps source compatibility, not binary
    compatibility
  - what counts as public API, and what not to rely on
  - deprecation in 1.x, removal only in 2.0
  - the supported platforms
- **[Accuracy and validation](docs/accuracy.md)** completed [#21]:
  - a summary by module
  - every number re-measured against today's tests (several tolerances had been tightened since)
  - rows for every component checked against an independent reference
  - separate lists of what is checked only for consistency and what is not yet validated
- Study `astro_study_low_precision`: the low-precision Sun and Moon against the full DE440s kernel,
  1950-2050. It confirms the documented 30" and 6' (measured 28.6" and 5.76'). [#21]
- **[Data files](docs/guide/data-files.md)**: what is built in, where each file comes from, how often
  it changes, what is lost without it, how to load it, and how each loader behaves when its data runs
  out (EOP throws; a leap-second table goes quietly stale, so check its expiry at start-up). [#21]

### Build and packaging

- **The shared library exports only the public API.** It is built with hidden visibility; public
  functions and classes are marked `ASTRO_EXPORT`, from the generated header `astro/export.h`. [#20]
  - No source changes for users.
  - Internal symbols (`detail` namespaces, implementation classes) are no longer exported: 501
    symbols instead of 625.
  - `tools/check_exports.py` (a test in shared builds) fails if an internal symbol leaks.
- **The soname carries the minor version**, `libastro.so.MAJOR.MINOR` (`libastro.so.0.3` for this
  release), because a minor release may change the binary interface (see the compatibility policy).
  Programs linked against a shared 0.2 build must be relinked. [#20]
- **`astro/astro.h`** includes every public header, for convenience. A test checks that it is
  complete, and the build checks that every public header compiles on its own. [#20]
- CI builds and tests a shared library too, and builds a program against the installed result.
  [#20]

## [0.2.0] - 2026-11-14

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

[0.3.0]: https://github.com/wgvisots/astro/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/wgvisots/astro/compare/db6a4b2...v0.2.0
[0.1.0]: https://github.com/wgvisots/astro/pull/1
