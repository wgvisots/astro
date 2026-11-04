# astro documentation

astro is a C++20 library for time scales, geodesy, celestial and terrestrial reference frames,
ephemerides, astrometry, orbit propagation (two-body, analytic and numerical J2/J4, SGP4), state
covariance, ballistic impact prediction and satellite geometry. Units are SI and angles are radians throughout. The only
dependency is header-only Boost (Describe).

## Where to start

| If you want to... | Read |
|---|---|
| build it, link it, and run a first program | [Getting started](guide/getting-started.md) |
| convert between UTC, TAI, TT, TDB, UT1, GPS, Unix time | [Time scales](guide/time.md) |
| go between latitude/longitude/height, ECEF and local ENU / az-el | [Geodesy](guide/geodesy.md) |
| go between ECEF and ECI (ITRS, GCRS, TEME, TOD, ...) | [Reference frames](guide/frames.md) |
| get Sun and Moon positions, rise/set times, apparent places | [Ephemerides and astrometry](guide/ephemerides-astrometry.md) |
| propagate an orbit (Kepler, J2/J4, SGP4) quickly and accurately | [Propagation](guide/propagation.md) |
| carry a covariance along, or follow an upstream filter's estimates | [Covariance and tracking](guide/covariance-tracking.md) |
| find where and when a ballistic trajectory hits the ground, with its error ellipse | [Impact prediction](guide/impact.md) |
| compute passes, look angles, Doppler, eclipses, orbital frames | [Satellite geometry](guide/satellite-geometry.md) |
| know the library's conventions (units, rotations, enums, errors, threads) | [Conventions](guide/conventions.md) |

## Background

- [Accuracy](accuracy.md): what each part is validated against, and measured accuracy of the
  propagators and the consistency checks, with the programs that produce the numbers.
- [Performance](performance.md): timings, and how to get the most out of the propagators.
- [References](references.md): the papers, standards and books the implementations follow.
- [Glossary](glossary.md): every acronym and term used here, spelled out.

The API reference is generated from the headers with Doxygen (`cmake --build build --target docs`);
every page here is also part of that site. The code in these pages is copied from the programs in
[`examples/`](https://github.com/wgvisots/astro/blob/main/examples), which are built and run as tests, and a test checks that the copies
are current.
