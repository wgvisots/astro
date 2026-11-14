# Getting started

## Requirements

| | |
|---|---|
| Compiler | C++20 (tested with GCC 13 and Clang 18) |
| Build | CMake 3.20 or newer |
| Dependency | Boost 1.77 or newer, headers only (for Boost.Describe) |
| Tests | GoogleTest (an installed copy is used if found, otherwise CMake downloads it) |
| Docs (optional) | Doxygen 1.9 or newer; Graphviz for diagrams |

## Build, test, install

```sh
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
cmake --install build --prefix /your/prefix
```

| CMake option | Default | |
|---|---|---|
| `ASTRO_BUILD_TESTS` | on (top-level build) | unit tests |
| `ASTRO_BUILD_EXAMPLES` | on (top-level build) | the programs in `examples/`, run as tests |
| `ASTRO_BUILD_BENCHMARKS` | off | timing programs in `bench/` ([Performance](../performance.md)) |
| `ASTRO_BUILD_STUDIES` | off | accuracy studies in `studies/` ([Accuracy](../accuracy.md)) |
| `ASTRO_WARNINGS_AS_ERRORS` | off | `-Werror` (or `/WX`) |
| `ASTRO_BUILD_FUZZERS` | off | libFuzzer targets for the parsers in `fuzz/` (Clang only) |
| `ASTRO_COVERAGE` | off | line-coverage instrumentation (GCC); report with `tools/coverage_report.py` |

`cmake --build build --target docs` builds the HTML documentation into `build/docs/html` when
Doxygen is installed.

### Packages

`cpack -G DEB` or `cpack -G RPM`, run in a build directory, makes a Debian or RPM package of what
`cmake --install` installs: the library, the headers, the CMake package files and the licence files,
under `/usr`. It depends on the Boost headers (`libboost-dev`, `boost-devel`). The RPM generator
needs `rpmbuild`. The vendor, maintainer and home page come from the cache variables
`ASTRO_PACKAGE_VENDOR`, `ASTRO_PACKAGE_CONTACT` and `ASTRO_PACKAGE_HOMEPAGE`.

## Use it from your project

```cmake
find_package(astro REQUIRED)     # finds Boost as a dependency
target_link_libraries(myapp PRIVATE astro::astro)
```

or add the source tree with `add_subdirectory(astro)`. Headers are under `astro/<module>/`.

## A first program

This converts a UTC time to the other time scales:

<!-- snippet: examples/time_scales.cpp#convert -->
```cpp
// A UTC instant inside the leap second at the end of 2016.
const Time utc = Time::fromIso8601("2016-12-31T23:59:60.5Z");
std::cout << utc.iso8601(3) << '\n'                     // 2016-12-31T23:59:60.500Z
          << utc.to(TimeScale::TAI).iso8601(3) << '\n'  // 2017-01-01T00:00:36.500 TAI
          << utc.to(TimeScale::TT).iso8601(3) << '\n'   // TAI + 32.184 s
          << utc.to(TimeScale::TDB).iso8601(6) << '\n'  // TT + periodic terms (< 2 ms)
          << utc.to(TimeScale::GPS).iso8601(3) << '\n'; // TAI - 19 s
```

The examples in [`examples/`](https://github.com/wgvisots/astro/blob/main/examples) are complete programs; each guide page quotes one.

## External data

The library never touches the network. Built in: leap seconds through 2017 (TAI - UTC = 37 s),
precession-nutation series, and a low-precision Sun and Moon. For full accuracy you supply:

| Data | Needed for | Load with |
|---|---|---|
| Earth orientation (`finals2000A.all` or `eopc04.1962-now`) | UT1 and ECEF <-> ECI at the metre level | `EopTable::loadFinals2000A` / `loadC04` |
| Leap seconds (`Leap_Second.dat`, `leap-seconds.list`) | UTC after a future leap second | `setLeapSecondTable(LeapSecondTable::load(...))` |
| JPL ephemeris (`de440s.bsp`) | high-precision Sun, Moon, planets | `SpkEphemeris` |

`scripts/fetch_iers_data.sh -o data` downloads the first two (`--de440s` adds the ephemeris);
each file is checked before it replaces the previous copy, so it is safe to run from cron.
