# astro

A C++20 library for time scales, geodesy, reference frames, ephemerides and orbits: IAU 2006/2000A
frames, SGP4, numerical and analytic propagation with covariance, Lambert's problem, conjunction
assessment and collision probability, impact prediction, and satellite geometry.

It began with conversions between celestial and terrestrial reference systems: GCRS (the Geocentric
Celestial Reference System, the inertial frame with its origin at the Earth's centre of mass) and
ITRS, or in short, ECI <--> ECEF. It has grown from there.

This project is based on the "Explanatory Supplement to the Astronomical Almanac 3rd Edition" {Urban and Seidelmann}.  This is an excellent book, and I highly recommend it.  All errors in the project are my own.

## What's in it

A C++20 library whose only dependency is header-only Boost (≥ 1.77, for Boost.Describe). SI units and
radians everywhere. Validated against ERFA/SOFA, JPL DE440 and Vallado's SGP4 verification set; see
[Accuracy](docs/accuracy.md).

| Module | Highlights |
|---|---|
| `core` | `Vec3`/`Mat3`/`Mat6`/`Quat`, passive `R1/R2/R3`, all 12 Euler sequences, angle wrapping, HMS/DMS parse/format, root finders, constants (WGS84, GRS80, IERS 2010, EGM2008 zonals, DE440 GMs) |
| `time` | Two-part Julian dates; `Time` over UTC, TAI, TT, TDB, TCG, TCB, UT1, GPS, Galileo, BeiDou; leap seconds (built in, or loaded from `Leap_Second.dat` / `leap-seconds.list` / `tai-utc.dat`); Unix and `std::chrono`; GNSS week/TOW with rollover; ISO-8601 |
| `geodesy` | Geodetic <-> ECEF (closed form, exact everywhere), ENU/NED, look angles (azimuth, elevation, range, range rate), geocentric/reduced latitude, Vincenty geodesics, 7/14-parameter Helmert |
| `frames` | IAU 2006/2000A (and 2000B) precession-nutation, ERA/GMST/GAST, polar motion; GCRS, J2000, MOD, TOD, TEME, CIRS, TIRS, PEF, ITRS with position+velocity; fast GCRS <-> ITRS over a time window; IERS EOP from `finals2000A` or C04 |
| `ephemeris` | Low-precision Sun/Moon (no data files); JPL SPK (`.bsp`, e.g. DE440) reader |
| `astrometry` | RA/Dec, ecliptic, galactic; apparent/observed places (space motion, parallax, light deflection, aberration, refraction); sidereal time; rise/set/transit, twilight, Moon phase |
| `orbit` | Keplerian and equinoctial elements; Kepler's equation; universal-variable propagation; covariance through two-body (closed-form state transition matrix), secular and numerical propagators; Lambert (Izzo, multi-rev); J2/J2²/J4 secular rates, sun-synchronous and repeat-track design; zonal gravity, RK4/RKF7(8), numerical propagator with event detection; 12-component (position, velocity, acceleration, jerk) propagation with optional unmodeled force; closest approach, conjunction lists and collision probability for any two of elements, states or trajectories, or in one call from states and covariances; ballistic impact on the WGS84 ellipsoid with impact error ellipse; TLE parse/format and SGP4/SDP4 |
| `satellite` | Doppler; pass prediction; estimate tracking with consistency checks; eclipses (cylindrical/conical) and beta angle; sub-satellite point, ground track, footprint, line of sight; RIC/LVLH/VNC/NTW frames, relative motion |

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

To make a Debian or RPM package of the library, headers and CMake files instead (it installs under
`/usr` and depends on the Boost headers):

```sh
cmake -B build -DCMAKE_BUILD_TYPE=Release -DASTRO_BUILD_TESTS=OFF -DASTRO_BUILD_EXAMPLES=OFF
cmake --build build -j
cd build && cpack -G DEB     # astro_<version>_<arch>.deb; or -G RPM (needs rpmbuild)
```

## A taste

```cpp
#include "astro/core/angle.h"
#include "astro/frames/eop.h"
#include "astro/frames/transform.h"
#include "astro/geodesy/geodetic.h"
#include "astro/orbit/numerical.h"
using namespace astro;

// WGS84 latitude, longitude (rad) and height (m) -> ECEF (m) and back.
const Vec3     ecef = geodeticToEcef({degToRad(38.8895), degToRad(-77.0353), 100.0});
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

See [CONTRIBUTING.md](CONTRIBUTING.md): how to build and test, what a change needs, the style, and
what CI checks. Changes between releases are listed in [CHANGELOG.md](CHANGELOG.md).

## License

Apache License 2.0 ([LICENSE](LICENSE)). Third-party coefficient tables, test data and code that
parts of the library follow are credited in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
