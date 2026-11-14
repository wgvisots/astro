# Ephemerides and astrometry

## Sun, Moon and planets

Every ephemeris implements `Ephemeris::state(target, center, t)` (and `position`), returning metres
and m/s in GCRS axes. Bodies use their NAIF ids (`Body::Moon` is 301, `Body::Earth` 399, ...).

| Ephemeris | Accuracy | Needs |
|---|---|---|
| `LowPrecisionEphemeris` | Sun better than 30" and 6e-5 in distance; Moon better than 6' and 0.14% (1950-2050, measured against DE440) | nothing |
| `SpkEphemeris` | that of the kernel (reads it to 1e-12 relative) | a JPL SPK file such as `de440s.bsp` |

The low-precision Sun follows Meeus [Meeus (1998)](../references.md) with a correction from the Earth-Moon
barycentre to the Earth; the Moon is the truncated series in Montenbruck & Gill
[Montenbruck & Gill (2000)](https://doi.org/10.1007/978-3-642-58351-3). `SpkEphemeris` reads JPL SPK kernels (DAF files, segment
types 2 and 3, either byte order) [NAIF SPK Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html) such as DE440
[Park et al. (2021)](https://doi.org/10.3847/1538-3881/abd414), chaining segments (Moon -> Earth-Moon barycentre -> Earth, ...) as
needed. It matches the reference `jplephem` reader to 1e-12 relative.

A kernel has only the bodies it has. DE440, for one, gives the barycentre of the Mars system
(`Body::MarsBarycenter`) but not Mars itself (`Body::Mars`); asking for the planet throws, with a
message that names the barycentre. One is never returned in place of the other.

<!-- snippet: examples/sun_moon.cpp#ephemeris -->
```cpp
// Low precision (Sun < 30", Moon < 6') needs no files; a JPL kernel gives the full DE440 accuracy.
std::shared_ptr<const Ephemeris> eph = std::make_shared<LowPrecisionEphemeris>();
if (argc > 1) {
    eph = std::make_shared<SpkEphemeris>(argv[1]);
}
const Time t   = Time::fromIso8601("2024-02-01T00:00:00Z");
const Vec3 sun = eph->position(Body::Sun, Body::Earth, t); // geocentric, m
std::cout << "Sun distance " << std::fixed << std::setprecision(6) << norm(sun) / constants::au << " au\n";
```

`scripts/fetch_iers_data.sh --de440s` downloads `de440s.bsp` (32 MB, 1849-2150). Times may be in
any scale; the reader converts to TDB.

## Rise, set, transit, twilight, Moon phase

<!-- snippet: examples/sun_moon.cpp#rise_set -->
```cpp
const Geodetic site{degToRad(51.4779), 0.0, 45.0}; // Greenwich
for (const SkyEvent& e : sunRiseSet(site, t, t + 86400.0, ConstantEop{}, *eph)) {
    const char* what = e.type == SkyEventType::Rise ? "sunrise "
                     : e.type == SkyEventType::Set  ? "sunset  "
                                                    : "noon    ";
    std::cout << what << e.time.iso8601(0) << '\n';
}
const MoonPhase phase = moonPhase(t, *eph);
std::cout << "Moon " << std::setprecision(0) << phase.illuminatedFraction * 100 << "% illuminated, "
          << (phase.waxing ? "waxing" : "waning") << '\n';
```

`sunRiseSet` uses the conventional -0.8333 deg altitude of the Sun's centre (refraction plus the
semi-diameter); `twilightEvents` takes `Twilight::Civil`, `Nautical` or `Astronomical` (-6, -12, -18 deg).
`moonRiseSet` uses the Moon's upper limb with 34' of refraction, as seen from the site (so the
parallax is included). `sunRiseSet` and `moonRiseSet` also return the transits.
`findElevationEvents` finds rises, sets and transits of any elevation function: it samples it and
refines each event with Brent's method to about 1 ms.

## Apparent and observed places

For stars (`CatalogStar`: ICRS position, proper motion, parallax, radial velocity) and solar-system
bodies, the library computes the full chain SOFA uses [SOFA](https://www.iausofa.org/): space motion, parallax,
light deflection by the Sun, annual and diurnal aberration, precession-nutation, Earth rotation,
polar motion, and optionally refraction:

- `starApparentGcrs` gives the apparent geocentric direction;
- `starObserved` and `bodyObserved` give a `LookAngles` at a site: azimuth, elevation, and range
  for bodies (0 for stars), with refraction when an `Atmosphere` is passed. They compute no
  velocity, so the range rate is NaN; so does `gcrsToTopocentric`, the geometric direction.

For stars the result agrees with ERFA's `atco13` to 1 mas without refraction, and above 10 deg
with it ([Accuracy](../accuracy.md)). `bodyObserved` takes the same steps with the light time in
place of the space motion; the tests check it for consistency between ephemerides, not against an
outside reference. Two small differences from SOFA: the proper-motion interval leaves out the
light-time (Roemer) term, worth up to 0.15 mas for the fastest stars, and the limit that keeps the
light deflection finite at the Sun's limb is fixed.

Refraction: `refractTrueToObserved` / `refractObservedToTrue` use the A tan z + B tan^3 z model
with constants from pressure, temperature, humidity and wavelength (as SOFA's `refco`);
`refractionBennett` [Bennett (1982)](https://doi.org/10.1017/S0373463300022037) and `refractionSaemundsson`
[Saemundsson (1986), as given in Meeus (1998)](../references.md) are the simple formulas good to about 0.1' down to the horizon.

The A/B model is meant for zenith distances up to about 75-80 deg. Lower down it is not the
conventional refraction: at a true elevation of 0 it gives about 11' where the usual figure is
34', and below about 5 deg its two directions are not exact inverses of each other (27" apart
at 0-2 deg). Use the simple formulas near the horizon.

Coordinates: `vectorToRaDec`, `angularSeparation`, `positionAngle`, ecliptic conversions with the
IAU 2006 obliquity, and ICRS <-> galactic with the Hipparcos matrix [ESA (1997)](https://ui.adsabs.harvard.edu/abs/1997ESASP1200.....E).
Sidereal time: `localMeanSiderealTime`, `localApparentSiderealTime`, `hourAngle`.
