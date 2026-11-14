# Satellite geometry

## From two-line elements to positions

`Tle::parse` reads a two-line element set (checksums verified, Alpha-5 catalogue numbers supported);
`Sgp4` propagates it [Vallado et al. (2006)](https://celestrak.org/publications/AIAA/2006-6753/). SGP4 works in TEME
[Seago & Vallado (2000)](https://doi.org/10.2514/6.2000-4025): convert before mixing with other data.

<!-- snippet: examples/sgp4_passes.cpp#sgp4 -->
```cpp
// ISS element set (the example from the Wikipedia article on two-line elements).
const Tle         tle = Tle::parse("1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927",
                           "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537",
                           "ISS (ZARYA)");
const Sgp4        sgp4(tle);
const ConstantEop eop;
const Time        t    = sgp4.epoch() + 3600.0;
const StateVector teme = sgp4.propagate(t); // TEME, m and m/s
const StateVector gcrs = temeToGcrs(teme, t, eop);
std::cout << tle.name << " at " << t.iso8601(0) << "\n  TEME " << teme.r << " m\n  GCRS " << gcrs.r << " m\n";
```

Things to know about `Sgp4`:

- **Time.** `propagate(dt)` takes seconds since the element set's epoch, as every propagator in the
  library does. The reference code and most other SGP4 implementations take minutes: for those, use
  `propagateMinutes(minutes)`. `propagate(Time)` uses the elapsed time, leap seconds included; a
  tool that subtracts Julian dates instead is 1 s (7.7 km in LEO) different across a leap second.
- **Epoch.** The day fraction of an element set is a fraction of 86 400 s, also on a day with a
  leap second, as in the reference code.
- **Constants.** The default is WGS 72 (`Sgp4Gravity::WGS72`), with which element sets are made.
- **Mode.** The default is `Sgp4OpsMode::Improved`, the mode of the verification output in this
  repository; `Afspc` reproduces the older operational behaviour.
- **Errors.** `propagate` throws `Sgp4Error` with the reference code's error number (1 to 6; 6 is a
  decayed satellite), or 7 for a time that is not finite. The overload
  `propagate(dt, StateVector&)` (and `propagateMinutes(minutes, StateVector&)`) returns the number
  instead of throwing.
- **Element sets.** `Tle::parseMany` reads a file of sets (with or without name lines);
  `line1()` and `line2()` write a set back and throw if a value does not fit its columns.

## Passes over a ground station

`findPasses` takes the satellite as a `Trajectory` (`orbit/trajectory.h`): a state as a function of
time that carries its frame. It needs the Earth-fixed ITRS, and throws `astro::Error` for any other
frame. `sgp4Trajectory(sgp4, eop, Frame::ITRS)` gives one from SGP4; any other source works as
`{Frame::ITRS, function}`, for example a `NumericalPropagator` whose states are converted with
`gcrsToItrs`. It samples the elevation and refines rise, set and maximum with Brent's method:

<!-- snippet: examples/sgp4_passes.cpp#passes -->
```cpp
// Passes above 10 deg over a station during the day after the epoch. findPasses needs the
// satellite in the Earth-fixed ITRS: sgp4Trajectory converts the SGP4 output (TEME) to it.
const Geodetic   station{degToRad(38.8895), degToRad(-77.0353), 100.0};
const Trajectory iss    = sgp4Trajectory(sgp4, eop, Frame::ITRS);
const auto       passes = findPasses(iss, station, sgp4.epoch(), sgp4.epoch() + 86400.0, degToRad(10.0));
for (const Pass& p : passes) {
    std::cout << "  AOS " << p.aos.iso8601(0) << "  max el " << std::fixed << std::setprecision(1)
              << radToDeg(p.maxElevation) << " deg  LOS " << p.los.iso8601(0) << '\n';
}
```

## Look angles, range rate and Doppler

<!-- snippet: examples/sgp4_passes.cpp#look -->
```cpp
// Look angles and the Doppler shift of a 437 MHz downlink at the first pass's closest approach.
if (!passes.empty()) {
    const LookAngles look = ecefToLookAngles(iss.state(passes[0].tca), station);
    std::cout << "  at maximum elevation: range " << look.range / 1e3 << " km, range rate " << look.rangeRate
              << " m/s, Doppler " << dopplerShift(437e6, look.rangeRate) << " Hz\n";
}
```

`ecefToLookAngles` (`astro/geodesy/local.h`) gives azimuth, elevation, range and range rate from a
site to a satellite's ITRS state (given only a position, it leaves the range rate NaN);
`dopplerShift` (first order) and `dopplerShiftRelativistic` turn range rate into a frequency shift.

## Eclipses and beta angle

<!-- snippet: examples/sun_moon.cpp#shadow -->
```cpp
// A satellite 7000 km from the Earth's center, directly behind it as seen from the Sun.
const ShadowResult s = earthShadow(-normalized(sun) * 7000e3, sun);
std::cout << "satellite: " << s.state << ", illumination " << s.illumination << '\n';
```

`earthShadow` uses a cylindrical or conical (umbra, penumbra, annular, with the visible fraction of the
solar disk) model [Montenbruck & Gill (2000), sec. 3.4](https://doi.org/10.1007/978-3-642-58351-3); `findEclipses` finds the
entry and exit times over a span; `betaAngle` is the Sun's elevation above the orbit plane.
`findEclipses(sat, sun, start, end)` takes two geocentric trajectories in the same frame (any frame:
it uses positions only), such as `sgp4Trajectory(sgp4, eop, Frame::GCRS)` and
`bodyTrajectory(eph, Body::Sun)`, which is in `Frame::GCRS`; different frames throw `astro::Error`.

## Ground track, footprint, line of sight

`subSatellitePoint` and `groundTrack` give the point below the satellite (`groundTrack` takes a
`Trajectory` in `Frame::ITRS`, as `findPasses` does); `coverageHalfAngle` and
`maxSlantRange` give the coverage for an altitude and minimum elevation, and `coverageFootprint` the circle
of that half angle as latitude/longitude points; `lineOfSight` says whether the Earth's ellipsoid
blocks the line between two points. All are in `satellite/ground.h`.

## Orbital frames and relative motion

`ricMatrix` (radial, in-track, cross-track; also called RSW), `lvlhMatrix`, `vncMatrix` and
`ntwMatrix` rotate inertial vectors into frames attached to a state; `rotateCovariance` rotates a
6x6 covariance. `relativeStateRic` and `relativeRicToDeputy` give the relative state of a deputy
in the chief's rotating RIC frame, the coordinates in which the Clohessy-Wiltshire equations
[Clohessy & Wiltshire (1960)](https://doi.org/10.2514/8.8704) are usually written today (check the axis order when
comparing with a particular text). The frame definitions are those of the CCSDS navigation
conventions [CCSDS 500.0-G-4](https://ccsds.org/Pubs/500x0g4.pdf).
