# Geodesy

All coordinates are on an `Ellipsoid` (WGS84 by default [NGA WGS 84](../references.md); GRS80 is also
provided). Latitudes and longitudes are in radians, heights in metres above the ellipsoid.

## Latitude, longitude, height <-> ECEF

<!-- snippet: examples/geodesy.cpp#ecef -->
```cpp
// WGS84 latitude, longitude (radians) and height above the ellipsoid (m) -> ECEF (m), and back.
const Geodetic site{deg2rad(38.8895), deg2rad(-77.0353), 100.0};
const Vec3     ecef = geodeticToEcef(site);
const Geodetic back = ecefToGeodetic(ecef);
std::cout << "ECEF " << ecef << " m\n";
std::cout << "lat " << rad2deg(back.lat) << " deg, lon " << rad2deg(back.lon) << " deg, h " << back.h << " m\n";
```

`ecefToGeodetic` uses Vermeille's closed form [Vermeille (2011)](https://doi.org/10.1007/s00190-010-0419-x), which is exact
everywhere outside the evolute of the ellipsoid (a small region about 43 km around the centre); inside it a bracketed
root search takes over, so every point has an answer. Round trips agree to 1e-8 m, and the result
matches ERFA's `gc2gd` to 2e-12 rad and 10 um near the surface ([Accuracy](../accuracy.md)).

Geocentric and reduced (parametric) latitudes convert with `geodeticToGeocentricLatitude`,
`geodeticToReducedLatitude` and their inverses.

## Local frames: ENU, NED, azimuth / elevation / range

`ecefToEnu` / `enuToEcef` (and the NED versions) express a point relative to a site; `enuToAer`
gives azimuth (clockwise from north), elevation and range:

<!-- snippet: examples/geodesy.cpp#aer -->
```cpp
// Where is a point 10 km east, 5 km north and 3 km up, as seen from the site?
const Vec3 target = enuToEcef({10e3, 5e3, 3e3}, site);
const Aer  aer    = enuToAer(ecefToEnu(target, site));
std::cout << "az " << rad2deg(aer.az) << " deg, el " << rad2deg(aer.el) << " deg, range " << aer.range << " m\n";
```

For a satellite's look angles and range rate from a ground site, see
[Satellite geometry](satellite-geometry.md).

## Distances on the ellipsoid

`geodesicInverse` and `geodesicDirect` solve the inverse (distance and azimuths between two
points) and direct (destination from a start, azimuth and distance) problems with Vincenty's
formulae [Vincenty (1975)](https://doi.org/10.1179/sre.1975.23.176.88), accurate to about 0.5 mm:

<!-- snippet: examples/geodesy.cpp#geodesic -->
```cpp
// Ellipsoidal distance and initial azimuth from Washington to Paris (Vincenty).
const GeodesicInverseResult g = geodesicInverse(site.lat, site.lon, deg2rad(48.8566), deg2rad(2.3522));
std::cout << "distance " << g.distance / 1e3 << " km, initial azimuth " << rad2deg(g.azimuth1) << " deg\n";
```

For nearly antipodal points Vincenty's inverse iteration can fail to converge; it then throws
`astro::Error` rather than return a wrong answer. `haversineDistance` is the spherical
approximation.

## Datum transformations

`helmertTransform` applies the IERS linearised 7-parameter similarity transformation
(`HelmertParams`: translation, rotation, scale) between terrestrial frames
[IERS Conventions (2010), ch. 4](https://www.iers.org/TN36). `HelmertTimeDependent` holds the 14-parameter form (the
parameters at an epoch plus their rates); its `at(year)` gives the 7 parameters for a date, and a
`helmertTransform` overload moves a position and velocity. Published parameter sets use either
the position-vector or the coordinate-frame sign convention for the rotations; set
`HelmertParams::convention` to match.
