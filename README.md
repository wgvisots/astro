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
| `orbit` | Keplerian and equinoctial elements; Kepler's equation; universal-variable propagation; Lambert (Izzo, multi-rev); J2/J2²/J4 secular rates, sun-synchronous and repeat-track design; zonal gravity, RK4/RKF7(8), numerical propagator; TLE parse/format and SGP4/SDP4 |
| `satellite` | Look angles, range rate, Doppler; pass prediction; eclipses (cylindrical/conical) and beta angle; sub-satellite point, ground track, footprint, line of sight; RIC/LVLH/VNC/NTW frames, relative motion, covariance rotation |

## Building

```sh
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
cmake --install build --prefix /your/prefix
```

Consumers use `find_package(astro)` and link `astro::astro` (Boost is found as a transitive dependency). Tests use an installed GoogleTest if one
is found, otherwise CMake fetches it. Set `-DASTRO_BUILD_TESTS=OFF` to skip them.

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
```

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
