# Ephemerides and astrometry

## Sun, Moon and planets

Every ephemeris implements `Ephemeris::state(target, center, t)` (and `position`), returning metres
and m/s in GCRS axes. Bodies use their NAIF ids (`Body::Moon` is 301, `Body::Earth` 399, ...).

| Ephemeris | Accuracy | Needs |
|---|---|---|
| `LowPrecisionEphemeris` | Sun better than 30" and 6e-5 in distance; Moon better than 6' and 0.14% (1950-2050, measured against DE440) | nothing |
| `SpkEphemeris` | that of the kernel (reads it to 1e-12 relative) | a JPL SPK file such as `de440s.bsp` |

The low-precision Sun follows Meeus [[ref:meeus|Meeus (1998)]] with a correction from the Earth-Moon
barycentre to the Earth; the Moon is the truncated series in Montenbruck & Gill
[[ref:montenbruck|Montenbruck & Gill (2000)]]. `SpkEphemeris` reads JPL SPK kernels (DAF files, segment
types 2 and 3, either byte order) [[ref:naif-spk|NAIF SPK Required Reading]] such as DE440
[[ref:de440|Park et al. (2021)]], chaining segments (Moon -> Earth-Moon barycentre -> Earth, ...) as
needed. It matches the reference `jplephem` reader to 1e-12 relative.

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
const Geodetic site{deg2rad(51.4779), 0.0, 45.0}; // Greenwich
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
semi-diameter); `twilight` takes `Twilight::Civil`, `Nautical` or `Astronomical` (-6, -12, -18 deg).
`findElevationEvents` finds rises, sets and transits of any elevation function: it samples it and
refines each event with Brent's method to about 1 ms.

## Apparent and observed places

For stars (`CatalogStar`: ICRS position, proper motion, parallax, radial velocity) and solar-system
bodies, the library computes the full chain SOFA uses [[ref:sofa|SOFA]]: space motion, parallax,
light deflection by the Sun, annual and diurnal aberration, precession-nutation, Earth rotation,
polar motion, and optionally refraction:

- `starApparentGcrs` gives the apparent geocentric direction;
- `starObserved` and `bodyObserved` give azimuth, elevation (and range for bodies) at a site,
  with refraction when an `Atmosphere` is passed.

They agree with ERFA's `atco13` to 1 mas without refraction and above 10 deg with it
([Accuracy](../accuracy.md)).

Refraction: `refractTrueToObserved` / `refractObservedToTrue` use the A tan z + B tan^3 z model
with constants from pressure, temperature, humidity and wavelength (as SOFA's `refco`);
`refractionBennett` [[ref:bennett1982|Bennett (1982)]] and `refractionSaemundsson`
[[ref:saemundsson1986|Saemundsson (1986)]] are the simple formulas good to about 0.1' down to the horizon.

Coordinates: `raDecFromVector`, `angularSeparation`, `positionAngle`, ecliptic conversions with the
IAU 2006 obliquity, and ICRS <-> galactic with the Hipparcos matrix [[ref:hipparcos|ESA (1997)]].
Sidereal time: `localMeanSiderealTime`, `localApparentSiderealTime`, `hourAngle`.
