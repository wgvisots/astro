# Satellite geometry

## From two-line elements to positions

`Tle::parse` reads a two-line element set (checksums verified, Alpha-5 catalogue numbers supported);
`Sgp4` propagates it [[ref:vallado2006|Vallado et al. (2006)]]. SGP4 works in TEME
[[ref:teme|TEME definition]]: convert before mixing with other data.

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

## Passes over a ground station

`findPasses` takes any function giving the satellite's ITRS state at a time, so it works with SGP4,
`NumericalPropagator` or an ephemeris. It samples the elevation and refines rise, set and maximum with
Brent's method:

<!-- snippet: examples/sgp4_passes.cpp#passes -->
```cpp
// Passes above 10 deg over a station during the day after the epoch.
const Geodetic station{deg2rad(38.8895), deg2rad(-77.0353), 100.0};
const auto     passes =
    findPasses(sgp4ItrsFunction(sgp4, eop), station, sgp4.epoch(), sgp4.epoch() + 86400.0, deg2rad(10.0));
for (const Pass& p : passes) {
    std::cout << "  AOS " << p.aos.iso8601(0) << "  max el " << std::fixed << std::setprecision(1)
              << rad2deg(p.maxElevation) << " deg  LOS " << p.los.iso8601(0) << '\n';
}
```

## Look angles, range rate and Doppler

<!-- snippet: examples/sgp4_passes.cpp#look -->
```cpp
// Look angles and the Doppler shift of a 437 MHz downlink at the first pass's closest approach.
if (!passes.empty()) {
    const LookAngles look = lookAngles(temeToItrs(sgp4.propagate(passes[0].tca), passes[0].tca, eop), station);
    std::cout << "  at TCA: range " << look.range / 1e3 << " km, range rate " << look.rangeRate << " m/s, Doppler "
              << dopplerShift(437e6, look.rangeRate) << " Hz\n";
}
```

`lookAngles` gives azimuth, elevation, range and range rate from a site to a satellite's ITRS state;
`dopplerShift` (first order) and `dopplerShiftRelativistic` turn range rate into a frequency shift.

## Eclipses and beta angle

<!-- snippet: examples/sun_moon.cpp#shadow -->
```cpp
// A satellite 7000 km from the Earth's center, directly behind it as seen from the Sun.
const ShadowResult s = shadow(-normalized(sun) * 7000e3, sun);
std::cout << "satellite: " << s.state << ", illumination " << s.illumination << '\n';
```

`shadow` uses a cylindrical or conical (umbra, penumbra, annular, with the visible fraction of the
solar disk) model [[ref:montenbruck|Montenbruck & Gill (2000), sec. 3.4]]; `findEclipses` finds the
entry and exit times over a span; `betaAngle` is the Sun's elevation above the orbit plane.

## Ground track, footprint, line of sight

`subSatellitePoint` and `groundTrack` give the point below the satellite; `coverageHalfAngle` and
`maxSlantRange` give the coverage for an altitude and minimum elevation, and `footprint` the circle
of that half angle as latitude/longitude points; `lineOfSight` says whether the Earth's ellipsoid
blocks the line between two points. All are in `satellite/ground.h`.

## Orbital frames and relative motion

`ricMatrix` (radial, in-track, cross-track; also called RSW), `lvlhMatrix`, `vncMatrix` and
`ntwMatrix` rotate inertial vectors into frames attached to a state; `rotateCovariance` rotates a
6x6 covariance. `relativeStateRic` and `deputyFromRelativeRic` give the relative state of a deputy
in the chief's rotating RIC frame, the coordinates of the Clohessy-Wiltshire equations
[[ref:cw1960|Clohessy & Wiltshire (1960)]].
