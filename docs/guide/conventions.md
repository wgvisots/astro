# Conventions

## Units and angles

SI units throughout: metres, seconds, m/s, m/s^2, m^3/s^2 for GM. Angles are radians; convert at the
edges with `degToRad`, `radToDeg`, `arcsecToRad`, and parse or format with `parseHms`, `parseDms`,
`formatHms`, `formatDms`.

<!-- snippet: examples/conventions.cpp#units -->
```cpp
// SI units and radians everywhere; convert at the edges.
const double lat = degToRad(38.8895);
std::cout << formatDms(lat) << '\n';
```

A few inputs keep the unit of their source instead (a reference code, a file format or a catalogue);
each says so in its name or its documentation:

| Where                                                                       | Unit                                                       |
|-----------------------------------------------------------------------------|------------------------------------------------------------|
| `Sgp4::propagateMinutes(minutes)`                                           | minutes since the TLE epoch, as in the SGP4 reference code |
| `Tle::meanMotion`, `meanMotionDot`, `meanMotionDdot`                        | rev/day, rev/day^2, rev/day^3, as printed in the TLE       |
| `Tle::bstar`                                                                | 1 / Earth radii                                            |
| `tdbMinusTt(tt, ut1Fraction, elong, uKm, vKm)`                              | observer distances `uKm`, `vKm` in km (as in ERFA)         |
| `lightDeflection(p, q, e, emAu, massSolar)`                                 | distance in au, mass in solar masses (as in ERFA)          |
| `aberration(p, vOverC, sunDistanceAu)`                                      | velocity in units of c, distance in au (as in ERFA)        |
| `Atmosphere{pressureHpa, temperatureC, relativeHumidity, wavelengthMicron}` | hPa, degrees Celsius, 0 to 1, micrometres                  |
| `refractionBennett`, `refractionSaemundsson(el, pressureHpa, temperatureC)` | hPa, degrees Celsius                                       |
| `CatalogStar::pmRaCosDec`, `pmDec`                                          | rad per Julian year                                        |
| `HelmertTimeDependent` rates, `epoch`; `helmertTransform(..., decimalYear)` | per year; decimal years                                    |
| `EopProvider::at(mjdUtc)`, `LeapSecondTable::taiMinusUtc(mjdUtc)`           | Modified Julian Date (MJD) in UTC, days                    |
| `CatalogStar::epochJdTdb`                                                   | Julian Date in TDB, days                                   |
| `SpkEphemeris::stateById(target, center, et)`                               | ephemeris time: TDB seconds since J2000.0                  |
| precession-nutation functions (`frames/precession_nutation.h`), `t`         | Julian centuries of TT since J2000.0                       |

## Names

Names follow one pattern, so that a reader who has learned one part of the library can guess the rest:

- **Conversions are named `xToY`**, source first: `degToRad`, `geodeticToEcef`, `vectorToRaDec`,
  `matrixToCipXY`. A static factory is named `from...` (`Time::fromMjd`, `Quat::fromMatrix`), and
  text is read by a `parse...` function (`parseHms`, `parseTimeScale`).
- **Constants live in `astro::constants` and are `camelCase`**: `constants::gmEarth`,
  `constants::earthRotationRateEra`, `constants::criticalInclination`, `constants::sunriseElevation`.
  The angle factors (`constants::degToRad`, a number) share their names with the functions that
  apply them (`astro::degToRad(x)`). Standard symbols keep their usual form (`c`, `J2`, `LG`, `TDB0`).
- **Accessors are nouns; `to` means a conversion** to another type or scale: `Time::calendar()`,
  `Time::sysTime()`, `Time::unixSeconds()`, `JulianDate::calendar()`, but `Time::to(TimeScale)`,
  `Quat::toMatrix()` and `toString(TimeScale)`.
- **Data, then time.** A transformation takes the data first and the time after it:
  `gcrsToItrs(r, t, eop)` and `ItrsRotationSeries::gcrsToItrs(r, dt)`, as `propagateKepler(s0, dt)`.
  (Force models, the right-hand side of an ordinary differential equation, keep the usual `f(t, y)`
  order: `ForceModel::acceleration(t, r, v)`.)
- **No generic free names in `astro::`.** A free function says what it computes (`earthShadow`,
  `twilightEvents`, `coverageFootprint`, `trackerPropagator`); an operation that belongs to one type
  is a member (`HelmertParams::inverse()`).
- **A constructor from one number or state is `explicit`** (`JulianDate(jd1, jd2 = 0)`,
  `SecularPropagator(s0, ...)`), so a bare value never converts silently. (A `Trajectory` is made
  implicitly from a function of time, so that one can be passed where a trajectory is expected.)
- **Every struct member has a default value**, so a default-constructed struct holds no garbage.

## Rotations

Rotation matrices are passive (frame) rotations, as in the Explanatory Supplement
[Urban & Seidelmann (2013)](https://aa.usno.navy.mil/publications/exp_supp) and SOFA: `v_new = R3(theta) * v_old` re-expresses a fixed vector in
axes turned by theta. `Quat` converts to and from matrices (Shepperd's method
[Shepperd (1978)](https://doi.org/10.2514/3.55767b)); matrices convert to and from all twelve Euler sequences
(`eulerToMatrix`, `matrixToEuler`), so an Euler sequence reaches a quaternion through its matrix.

## Times and frames are explicit

A `Time` carries its time scale and converts on demand. States carry no frame tag, so the API
asks for the frame where it matters (`NumericalPropagator::j2j4(..., Frame::GCRS)`,
`TrackerOptions::frame`).

## Propagators

Every propagator (`propagateKepler`, `SecularPropagator`, `NumericalPropagator`, `Sgp4`) follows one
pattern:

- **Times are seconds from the epoch.** A function named `propagate...` with a `double` time
  argument always takes seconds. The one exception is named for it: `Sgp4::propagateMinutes`, for
  the minutes of the SGP4 reference code.
- **Arguments come in one order:** the initial state (for the free functions), the times, the
  covariance, other inputs, the output buffer, then the options (arguments with defaults, such as
  `accelerationNoisePsd` and `AdaptiveStats*`).
- **The same names for the same results.** At one time, `propagateStm(dt)` gives the state and its
  state transition matrix and `propagateWithCovariance(dt, P0)` the state and its covariance
  (`propagateKeplerStm` and `propagateKeplerWithCovariance` for the two-body functions). The
  `propagateGrid...(start, step, ...)` forms give results at `start + k * step`.

## Enums

Every public enum is reflected with Boost.Describe and prints as `Enum::Value`:

<!-- snippet: examples/conventions.cpp#enums -->
```cpp
// Every public enum prints as Enum::Value and is reflected with Boost.Describe.
std::cout << TimeScale::TT << '\n';          // TimeScale::TT
std::cout << static_cast<Frame>(42) << '\n'; // Frame::INVALID_VALUE
boost::mp11::mp_for_each<boost::describe::describe_enumerators<Frame>>(
    [](auto d) { std::cout << d.name << (isEarthFixed(d.value) ? " (Earth-fixed)\n" : "\n"); });
```

## Errors

Bad input and missing data (no EOP for an epoch, a TLE checksum mismatch, a non-converging solution)
throw `astro::Error`, derived from `std::runtime_error`. Functions marked `noexcept` cannot fail.

<!-- snippet: examples/conventions.cpp#errors -->
```cpp
// Bad input and missing data throw astro::Error (a std::runtime_error).
try {
    (void)Time::fromIso8601("not a time");
} catch (const Error& e) {
    std::cout << "astro::Error: " << e.what() << '\n';
}
```

## 6x6 matrices

`Mat6` (`core/mat6.h`) is a `std::array<std::array<double, 6>, 6>` for covariances and state
transition matrices, with named helpers (`multiply`, `transpose`, `congruence` for A P A^T,
`block`, `mahalanobis2`): operators declared for a `std::array` would not be found outside the
`astro` namespace.

## Threads

The propagators, ephemerides and other objects are immutable after construction (the SPK reader
locks internally), so their const member functions may be called concurrently on a shared instance.
Propagation calls report integrator statistics through an optional `AdaptiveStats*` argument rather
than storing them. `setLeapSecondTable` swaps the table atomically. `EstimateTracker::update` needs
exclusive access, like a standard container. The test suite runs clean under ThreadSanitizer.
