# Conventions

## Units and angles

SI units throughout: metres, seconds, m/s, m/s^2, m^3/s^2 for GM. Angles are radians; convert at the
edges with `deg2rad`, `rad2deg`, `arcsec2rad`, and parse or format with `parseHms`, `parseDms`,
`formatHms`, `formatDms`.

<!-- snippet: examples/conventions.cpp#units -->
```cpp
// SI units and radians everywhere; convert at the edges.
const double lat = deg2rad(38.8895);
std::cout << formatDms(lat) << '\n';
```

## Rotations

Rotation matrices are passive (frame) rotations, as in the Explanatory Supplement
[Urban & Seidelmann (2013)](https://aa.usno.navy.mil/publications/exp_supp) and SOFA: `v_new = R3(theta) * v_old` re-expresses a fixed vector in
axes turned by theta. `Quat` converts to and from matrices (Shepperd's method
[Shepperd (1978)](https://doi.org/10.2514/3.55767b)) and to and from all twelve Euler sequences.

## Times and frames are explicit

A `Time` carries its time scale and converts on demand. States carry no frame tag, so the API
asks for the frame where it matters (`NumericalPropagator::j2j4(..., Frame::GCRS)`,
`TrackerOptions::frame`).

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
