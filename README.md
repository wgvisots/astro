# astro
Space and coordinate conversion tools.

This project mostly had to do with translations between Celestial Reference Systems and Terrestrial Reference Systems.  In particular, I try to use GCRS (Geocentric Celestial Reference System), that is the inertial coordinate system with its origin at the center of mass of the Earth.

In short, ECI <--> ECEF.

This project is based on the "Explanatory Supplement to the Astronomical Almanac 3rd Edition" {Urban and Seidelmann}.  This is an excellent book, and I highly recommend it.  All errors in the project are my own.

## What's in it

A C++20 library whose only dependency is header-only Boost (≥ 1.77, for Boost.Describe). SI units and
radians everywhere. Validated against ERFA/SOFA, JPL DE440 and Vallado's SGP4 verification set; see
[Accuracy](docs/accuracy.md).

| Module | Highlights |
|---|---|
| `core` | `Vec3`/`Mat3`/`Quat`, passive `R1/R2/R3`, all 12 Euler sequences, angle wrapping, HMS/DMS parse/format, root finders, constants (WGS84, GRS80, IERS 2010, EGM2008 zonals, DE440 GMs) |
| `time` | Two-part Julian dates; `Time` over UTC, TAI, TT, TDB, TCG, TCB, UT1, GPS, Galileo, BeiDou; leap seconds (built in, or loaded from `Leap_Second.dat` / `leap-seconds.list` / `tai-utc.dat`); Unix and `std::chrono`; GNSS week/TOW with rollover; ISO-8601 |
| `geodesy` | Geodetic <-> ECEF (closed form, exact everywhere), ENU/NED/AER, geocentric/reduced latitude, Vincenty geodesics, 7/14-parameter Helmert |
| `frames` | IAU 2006/2000A (and 2000B) precession-nutation, ERA/GMST/GAST, polar motion; GCRS, J2000, MOD, TOD, TEME, CIRS, TIRS, PEF, ITRS with position+velocity; fast GCRS <-> ITRS over a time window; IERS EOP from `finals2000A` or C04 |
| `ephemeris` | Low-precision Sun/Moon (no data files); JPL SPK (`.bsp`, e.g. DE440) reader |
| `astrometry` | RA/Dec, ecliptic, galactic; apparent/observed places (space motion, parallax, light deflection, aberration, refraction); sidereal time; rise/set/transit, twilight, Moon phase |
| `orbit` | Keplerian and equinoctial elements; Kepler's equation; universal-variable propagation; Lambert (Izzo, multi-rev); J2/J2²/J4 secular rates, sun-synchronous and repeat-track design; zonal gravity, RK4/RKF7(8), numerical propagator with event detection; 12-component (position, velocity, acceleration, jerk) propagation with optional unmodeled force; closest approach, conjunction lists and collision probability for any two of elements, states or trajectories; ballistic impact on the WGS84 ellipsoid with impact error ellipse; TLE parse/format and SGP4/SDP4 |
| `satellite` | Look angles, range rate, Doppler; pass prediction; estimate tracking with consistency checks; eclipses (cylindrical/conical) and beta angle; sub-satellite point, ground track, footprint, line of sight; RIC/LVLH/VNC/NTW frames, relative motion, covariance rotation |

## Documentation

- **[User guide](docs/index.md)**: getting started, time scales, geodesy, reference frames,
  ephemerides and astrometry, propagation, covariance and tracking, impact prediction, satellite
  geometry, conventions.
- **[Accuracy](docs/accuracy.md)** and **[Performance](docs/performance.md)**: what is validated
  against what, measured propagation and consistency-check accuracy, timings, and the programs
  that reproduce them.
- **[References](docs/references.md)**: the standards, papers and books the code follows.
- **[Glossary](docs/glossary.md)**: every acronym and term used in the docs, spelled out.
- **API reference**: `cmake --build build --target docs` builds it (with the pages above) into
  `build/docs/html`; it is also published at <https://wgvisots.github.io/astro/>.
- **[Examples](examples)**: small complete programs, built and run as tests.

## Building

```sh
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
cmake --install build --prefix /your/prefix
```

Consumers use `find_package(astro)` and link `astro::astro` (Boost is found as a transitive
dependency). See [Getting started](docs/guide/getting-started.md) for the CMake options.

## A taste

```cpp
#include "astro/core/angle.h"
#include "astro/frames/eop.h"
#include "astro/frames/transform.h"
#include "astro/geodesy/geodetic.h"
#include "astro/orbit/numerical.h"
using namespace astro;

// WGS84 latitude, longitude (rad) and height (m) -> ECEF (m) and back.
const Vec3     ecef = geodeticToEcef({deg2rad(38.8895), deg2rad(-77.0353), 100.0});
const Geodetic g    = ecefToGeodetic(ecef);

// Time scales.
const Time utc = Time::fromIso8601("2024-03-20T12:00:00Z");
const Time tt  = utc.to(TimeScale::TT);

// ECEF <-> ECI (ITRS <-> GCRS); load IERS Earth orientation data for full accuracy.
const EopTable    eop = EopTable::loadFinals2000A("finals2000A.all");
const StateVector eci = itrsToGcrs(StateVector{ecef, {}}, utc, eop);

// Numerical J2/J4 for 4 h, output every 10 s, from a GCRS StateVector satelliteState (m, m/s);
// the zonal field is aligned with the true pole.
AdaptiveOptions opt;
opt.relTol = 1e-9;                          // ~0.2 m over 4 h in LEO
const auto prop  = NumericalPropagator::j2j4(utc, satelliteState, Frame::GCRS, opt);
const auto track = prop.propagateGrid(0.0, 10.0, 1441);
```

## External data

The library never touches the network. `scripts/fetch_iers_data.sh -o data` downloads the IERS Earth
orientation and leap-second files (`--de440s` adds the JPL ephemeris); each file is checked before
it replaces the previous copy, so it is cron-friendly. See [Getting started](docs/guide/getting-started.md#external-data).

## Contributing

Code is formatted with clang-format 18 (`.clang-format`; `clang-format -i $(git ls-files '*.h' '*.cpp')`).
CI builds with GCC and Clang with warnings as errors, runs the tests under ThreadSanitizer, checks
formatting, and builds the documentation with Doxygen warnings as errors. The code in `docs/` must
match the examples: `tools/check_doc_snippets.py --fix` refreshes it.

## License

Apache License 2.0 ([LICENSE](LICENSE)). Third-party coefficient tables and test data are credited
in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
