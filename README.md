# astro
Space and coordinate conversion tools.

This project mostly had to do with translations between Celestial Reference Systems and Terrestrial Reference Systems.  In particular, I try to use GCRS (Geocentric Celestial Reference System), that is the inertial coordinate system with its origin at the center of mass of the Earth.

In short, ECI <--> ECEF.

This project is based on the "Explanatory Supplement to the Astronomical Almanac 3rd Edition" {Urban and Seidelmann}.  This is an excellent book, and I highly recommend it.  All errors in the project are my own.

## What's in it

A C++20 library whose only dependency is header-only Boost (≥ 1.77, for Boost.Describe). SI units and
radians everywhere.

| Module | Highlights |
|---|---|
| `core` | `Vec3`/`Mat3`/`Quat`, passive `R1/R2/R3`, all 12 Euler sequences, angle wrapping, HMS/DMS parse/format, root finders, constants (WGS84, GRS80, IERS 2010, EGM2008 zonals, DE440 GMs) |
| `time` | Two-part Julian dates; `Time` over UTC, TAI, TT, TDB, TCG, TCB, UT1, GPS, Galileo, BeiDou; leap seconds (built in, or loaded from `Leap_Second.dat` / `leap-seconds.list` / `tai-utc.dat`); Unix and `std::chrono`; GNSS week/TOW with rollover; ISO-8601 |
| `geodesy` | Geodetic <-> ECEF (closed form, exact everywhere), ENU/NED/AER, geocentric/reduced latitude, Vincenty geodesics, 7/14-parameter Helmert |
| `frames` | IAU 2006/2000A (and 2000B) precession-nutation, ERA/GMST/GAST, polar motion; GCRS, J2000, MOD, TOD, TEME, CIRS, TIRS, PEF, ITRS with position+velocity; IERS EOP from `finals2000A` or C04 |
| `ephemeris` | Low-precision Sun/Moon (no data files); JPL SPK (`.bsp`, e.g. DE440) reader |
| `astrometry` | RA/Dec, ecliptic, galactic; apparent/observed places (space motion, parallax, light deflection, aberration, refraction); sidereal time; rise/set/transit, twilight, Moon phase |
| `orbit` | Keplerian and equinoctial elements; Kepler's equation; universal-variable propagation; Lambert (Izzo, multi-rev); J2/J2²/J4 secular rates, sun-synchronous and repeat-track design; zonal gravity, RK4/RKF7(8), numerical propagator; 12-component (position, velocity, acceleration, jerk) propagation with optional unmodeled force; TLE parse/format and SGP4/SDP4 |
| `satellite` | Look angles, range rate, Doppler; pass prediction; estimate tracking with consistency checks; eclipses (cylindrical/conical) and beta angle; sub-satellite point, ground track, footprint, line of sight; RIC/LVLH/VNC/NTW frames, relative motion, covariance rotation |

## Building

```sh
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
cmake --install build --prefix /your/prefix
```

Consumers use `find_package(astro)` and link `astro::astro` (Boost is found as a transitive dependency). Tests use an installed GoogleTest if one
is found, otherwise CMake fetches it. Set `-DASTRO_BUILD_TESTS=OFF` to skip them.

Code is formatted with clang-format 18 or newer using `.clang-format`: `clang-format -i $(git ls-files '*.h' '*.cpp')`.
Generated tables listed in `.clang-format-ignore` are skipped.

## Examples

```cpp
#include "astro/core/angle.h"
#include "astro/frames/transform.h"
#include "astro/geodesy/geodetic.h"
#include "astro/orbit/sgp4.h"
#include "astro/satellite/look.h"
using namespace astro;

// WGS84 lat/lon/alt (rad, rad, m) -> ECEF (m) and back.
Vec3 ecef = geodeticToEcef({deg2rad(38.8895), deg2rad(-77.0353), 100.0});
Geodetic g = ecefToGeodetic(ecef);

// Time scales.
Time utc = Time::fromIso8601("2024-03-20T12:00:00Z");
Time tt  = utc.to(TimeScale::TT);
double unix = utc.unixSeconds();
GnssWeekTow gps = utc.gnssWeekTow(GnssSystem::GPS);

// ECEF <-> ECI (ITRS <-> GCRS). Use real EOP for full accuracy.
EopTable eop = EopTable::loadFinals2000A("finals2000A.all");
StateVector eci = itrsToGcrs(StateVector{ecef, {}}, utc, eop);
FrameTransformer ft(utc, eop);                        // reuse for many vectors at one epoch
Vec3 teme = ft.transform(eci.r, Frame::GCRS, Frame::TEME);

// TLE -> SGP4 -> passes over a ground station.
Tle tle = Tle::parse(line1, line2);
Sgp4 sgp4(tle);
auto passes = findPasses(sgp4ItrsFunction(sgp4, eop), g, utc, utc + 86400.0, deg2rad(10.0));

// 12-component state (m, m/s, m/s^2, m/s^3, inertial). ModelOnly propagates from r and v and
// reports how far the input a and j are from the model; UnmodeledForce keeps that difference as
// an extra force da0 + dj0 t (e.g. thrust). Output a and j are at the new time.
StateVector12 s12{r, v, a, j};
Propagation12 p = propagateKepler(s12, 600.0, constants::gmEarth, HigherOrderMode::UnmodeledForce);
Propagation12 q = propagateSecular(s12, 86400.0);                       // J2 + J2² + J4 mean drift
auto num = NumericalPropagator::j2j4(utc, s12, Frame::GCRS, HigherOrderMode::UnmodeledForce);
StateVector12 later = num.propagate12(3600.0);
```

## Performance

To propagate one orbit to many equally spaced times, use the grid APIs: `SecularPropagator`
(analytic J2/J4) and `NumericalPropagator::propagateGrid`. LEO, one day of output, one core of an
i7-1260P:

| Method                                           | Y = 60 s     | Y = 10 s     |
|--------------------------------------------------|--------------|--------------|
| `SecularPropagator::propagateGrid`               | 42 ns/output | 43 ns/output |
| `SecularPropagator::propagateGrid12`             | 70 ns/output | 71 ns/output |
| `NumericalPropagator::j2j4(...).propagateGrid`   | 0.51 ms/day  | 0.56 ms/day  |
| `NumericalPropagator::j2j4(...).propagateGrid12` | 0.56 ms/day  | 0.89 ms/day  |

The numerical grid lets the integrator take its own error-controlled steps (about 96 s in LEO at the
default tolerance) and interpolates the outputs with a degree-7 Hermite polynomial, so its cost
barely depends on Y. Its accuracy is set by `AdaptiveOptions`: the defaults give about 2 mm after a
day in LEO, and each 10x looser `relTol` cuts the run time by roughly a quarter and costs about
10x in accuracy. Run `bench/propagation_bench.cpp` (`-DASTRO_BUILD_BENCHMARKS=ON`) with your own X
and Y:

```sh
cmake -B build -DCMAKE_BUILD_TYPE=Release -DASTRO_BUILD_BENCHMARKS=ON && cmake --build build -j
build/bench/astro_propagation_bench 86400 10 100      # X seconds, Y seconds, number of states
```

### Covariance

`propagateGridWithCovariance` carries a 6x6 position/velocity covariance along the grid. The state
transition matrix is integrated with the state using exact gravity partials (same steps, so the
states are identical to `propagateGrid`), and P(t) = Phi P0 Phi^T + Q(t). Q is optional process
noise: white acceleration noise with a power spectral density in m^2/s^3. `propagateGridStm`
returns Phi itself. Expect about 6x the cost of the states alone (8x with process noise).

```cpp
AdaptiveOptions opt;
opt.relTol = 1e-9;   // ~0.2 m worst case over 4 h
opt.absTol = 1e-6;
auto prop = NumericalPropagator::j2j4(epoch, state, Frame::GCRS, opt);   // or Frame::TOD, ...
std::vector<StateCovariance> out = prop.propagateGridWithCovariance(0.0, 10.0, 1441, P0, 1e-9);
```

### Following an upstream filter

`EstimateTracker` (in `satellite/estimate_tracker.h`) keeps the latest delivered state and
covariance, and checks each new delivery against the previous one propagated to its epoch. Successive
outputs of one filter are correlated, so they are not fused; the newest replaces the previous.

```cpp
TrackerOptions options;
options.frame = Frame::TOD;                           // frame of the deliveries (default GCRS)
EstimateTracker tracker(options);                     // J2 + J4, relTol 1e-9
if (auto check = tracker.update({epoch, state, P}); check && !check->consistent) {
    // check->positionResidualRic vs check->positionSigmaRic says where and by how many sigma
}
auto ahead = tracker.predictGrid(0.0, 10.0, 1441);    // states + covariances, 4 h at 10 s
```

The check reports residuals and predicted 1-sigma along radial / in-track / cross-track, the
semi-major axis residual, and Mahalanobis distances with chi-square p-values for position (3-D),
position + velocity (6-D, Cartesian) and position + velocity in equinoctial elements (6-D).
`TrackerOptions::test` picks which one sets `consistent`. The default, `Elements`, kept its nominal
false-alarm rate in Monte Carlo runs out to 4 h for LEO, SSO, HEO, near-equatorial, polar and
retrograde orbits, and catches a burn even just before the delivery. The Cartesian 6-D test is as
sensitive but over-alarms over long gaps, and the position test misses fresh velocity changes.

## External data

The library never touches the network. `scripts/fetch_iers_data.sh` downloads the files below
(`--c04` and `--de440s` add the optional ones; `--max-age DAYS` makes it cron-friendly). Each file is
checked before it replaces the previous copy:

```sh
scripts/fetch_iers_data.sh -o data            # finals2000A.all, Leap_Second.dat, leap-seconds.list
scripts/fetch_iers_data.sh -o data --all      # plus eopc04.1962-now and de440s.bsp
```

You supply:
- **Leap seconds**: a table through 2017 (TAI−UTC = 37 s) is built in. Install a newer one with
  `setLeapSecondTable(LeapSecondTable::load("Leap_Second.dat"))`.
- **Earth orientation (UT1−UTC, polar motion, dX/dY)**: `EopTable::loadFinals2000A(...)` or `loadC04(...)`,
  from the IERS. You can also pass values directly with `ConstantEop`. `ConstantEop{}` (all zero) means
  UT1 = UTC with no polar motion, which is good to about 1 s of time (roughly 400 m at the equator).
- **High-precision Sun/Moon/planets**: a JPL SPK kernel such as `de440s.bsp` from NAIF, read by
  `SpkEphemeris`. `LowPrecisionEphemeris` needs no files (Sun < 30″, Moon < 6′).

## Conventions

- Rotation matrices are passive (frame) rotations, as in the Explanatory Supplement and SOFA:
  `v_new = R3(θ) v_old`.
- `Time` carries its scale. UTC uses the SOFA "quasi-JD" convention, so leap seconds such as
  `23:59:60.5` can be represented and arithmetic counts them.
- Every public enum is reflected with `BOOST_DESCRIBE_ENUM` and has an `operator<<` that writes
  `TimeScale::UTC` (or `TimeScale::INVALID_VALUE` for a value outside the enumeration).
- Zonal gravity (J2, J4) acts about the Earth's true spin axis. In GCRS/J2000 that axis is about
  0.15 deg from z (2026), which over 4 h in LEO would cost ~270 m, so say which frame a state is in:
  `NumericalPropagator::j2j4(epoch, state, Frame::GCRS)`, `TrackerOptions::frame`, or
  `poleOfDateMatrix` for a hand-built `GravityField`. TOD, TEME and CIRS need no rotation.
  `EstimateTracker` also expresses predictions in the of-date frame of their own time.
- Thread safety: propagators, ephemerides and the other library objects are immutable after
  construction (the SPK reader locks internally), so const member functions may be called from
  many threads on a shared instance. `setLeapSecondTable` swaps the table atomically.
  `EstimateTracker::update` needs exclusive access. ThreadSanitizer runs clean on the test suite.
- Celestial frames use the IERS 2010 CIO-based chain. The equinox-based frames are also available
  for interoperability.

## Validation

Tests compare against independent references:
- ERFA/SOFA for time scales, precession-nutation, GCRS→ITRS, apparent places, refraction and
  geodesy (generated by `tools/gen_reference.py`)
- JPL DE440 via jplephem
- Vallado's SGP4 verification set (all cases, within 0.12 mm)
- Worked examples from Vallado and Curtis

To regenerate the reference data, run the dev-only scripts in `tools/` (they need `numpy`, `pyerfa`,
`sgp4` and `jplephem`).

Third-party coefficient tables and test data are credited in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
