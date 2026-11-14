# Impact prediction

`predictImpact` follows a ballistic trajectory from an ECI (Earth-Centered Inertial) state until it
falls through the WGS84 ellipsoid, or a surface at a given height above it. It reports when and
where that happens and, if you pass the state's 6x6 covariance, the impact error ellipse.
Acronyms are spelled out in the [Glossary](../glossary.md).

## Predicting an impact

<!-- snippet: examples/impact.cpp#impact -->
```cpp
// 1-sigma: 100 m per axis in position, 0.1 m/s per axis in velocity (GCRS, m^2 and m^2/s^2).
Mat6 P0{};
for (std::size_t i = 0; i < 3; ++i) {
    P0[i][i]         = 100.0 * 100.0;
    P0[i + 3][i + 3] = 0.1 * 0.1;
}
const std::optional<Impact> impact = predictImpact(epoch, s0, P0, Frame::GCRS, eop);
if (!impact) {
    std::cout << "no impact within ImpactOptions::maxFlightTime\n";
    return 1;
}
const ErrorEllipse e1  = impact->uncertainty->ellipse; // 1-sigma: 39 % of impacts
const ErrorEllipse e95 = e1.scaled(0.95);              // 95 % of impacts
```

<!-- snippet: examples/impact.cpp#print -->
```cpp
std::cout << std::fixed << std::setprecision(3);
std::cout << "impact time      " << impact->time.iso8601() << '\n';
std::cout << "flight time      " << impact->flightTime << " s\n";
std::cout << "latitude         " << rad2deg(impact->point.lat) << " deg\n";
std::cout << "longitude        " << rad2deg(impact->point.lon) << " deg\n";
std::cout << "impact angle     " << rad2deg(impact->angle) << " deg below horizontal\n";
std::cout << "speed            " << impact->speed << " m/s (relative to the ground)\n";
std::cout << std::setprecision(1);
std::cout << "1-sigma ellipse  " << e1.semiMajor << " x " << e1.semiMinor << " m, rotation " << rad2deg(e1.rotation)
          << " deg (clockwise from north)\n";
std::cout << "95 % ellipse     " << e95.semiMajor << " x " << e95.semiMinor << " m\n";
std::cout << std::setprecision(3) << "time sigma       " << impact->uncertainty->timeSigma << " s\n";
```

Output of `examples/impact.cpp` (a launch from 30 N, 100 W at 5 km/s, 35 degrees above the horizon):

```
impact time      2026-03-01T12:15:36.097Z
flight time      936.097 s
latitude         46.833 deg
longitude        -67.926 deg
impact angle     35.008 deg below horizontal
speed            5000.011 m/s (relative to the ground)
1-sigma ellipse  295.7 x 92.2 m, rotation 69.9 deg (clockwise from north)
95 % ellipse     723.7 x 225.6 m
time sigma       0.087 s
```

The state may be in any celestial frame: GCRS, J2000 (EME2000), MOD, TOD, TEME or CIRS. Of-date
frames are taken at their own time, and `Impact::state` comes back in the same frame at the impact
time. `predictImpact` returns `std::nullopt` when nothing comes down within
`ImpactOptions::maxFlightTime` (2 h by default), for example for an orbit.

## What you get

| Field | Meaning |
|---|---|
| `time`, `flightTime` | impact time, and seconds since the initial epoch |
| `state` | position and velocity at impact, in the input frame |
| `itrs` | position and ground-relative velocity at impact, in ITRS (ECEF) |
| `point` | geodetic latitude, longitude (rad) and height (= `surfaceHeight`) |
| `speed` | speed relative to the ground (m/s) |
| `angle` | impact angle: flight path below the local horizontal (rad) |
| `heading` | direction of travel over the ground, clockwise from north (rad) |
| `uncertainty` | present when a covariance was given; see below |

## The impact error ellipse

`Impact::uncertainty` holds:

| Field | Meaning |
|---|---|
| `ellipse.semiMajor`, `ellipse.semiMinor` | 1-sigma semi-axes on the ground (m) |
| `ellipse.rotation` | direction of the semi-major axis, clockwise from north, 0 to 180 degrees (rad) |
| `eastVariance`, `northVariance`, `eastNorthCovariance` | the same, as a 2x2 covariance (m^2) |
| `timeSigma` | standard deviation of the impact time (s) |
| `covariance` | the full 6x6 state covariance at the nominal impact time, in the input frame |

The ellipse lies in the local horizontal (east/north) plane at the impact point. A 1-sigma ellipse
contains 39.3 % of impacts. `ellipse.scaled(p)` gives the ellipse that contains a fraction p, by
the factor k = sqrt(-2 ln(1 - p)) (chi-square with 2 degrees of freedom):

| Contains | 39.3 % | 50 % | 90 % | 95 % | 99 % |
|---|---|---|---|---|---|
| k | 1 | 1.1774 | 2.1460 | 2.4477 | 3.0349 |

How the ellipse is computed:

1. The covariance is propagated to the nominal impact time with the variational equations, as in
   [Covariance and tracking](covariance-tracking.md). `ImpactOptions::accelerationNoisePsd` adds
   white acceleration noise for forces the model leaves out.
2. A position error dr at that time moves the impact point along the ground-relative velocity v
   until the trajectory is back on the surface: dp = (I - v u^T / (u . v)) R dr, where R rotates
   GCRS into ITRS and u is the local vertical. Velocity errors enter through step 1.
3. The east and north parts of dp give the 2x2 covariance, and its eigenvalues and eigenvectors
   give the axes and rotation.

This is a linear (first-order) result. A Monte Carlo study (section "Impact prediction" of [Accuracy](../accuracy.md))
shows it holds well even for 1 km / 1 m/s errors: 94-95 % of sampled impacts fall inside the
linear 95 % ellipse. At shallow impact angles the along-track axis grows as 1 / tan(angle), and the
linear result is less reliable for grazing trajectories.

## How it works

| Step | What happens |
|---|---|
| Frame | the state is rotated into GCRS (`gcrsToCelestialMatrix`) |
| Force model | point mass + J2 + J3 + J4, with the zonal field aligned to the true pole at the epoch (`ImpactOptions::J`); extra forces through `ImpactOptions::additionalForces` |
| Integration | RKF7(8), relTol 1e-9 (`ImpactOptions::integrator`) |
| Height | each sample is rotated into ITRS with an `ItrsRotationSeries` and converted to geodetic height ([Vermeille (2011)](https://doi.org/10.1007/s00190-010-0419-x)) |
| Impact | `NumericalPropagator::propagateToEvent` samples the height four times per integrator step on the degree-7 Hermite interpolant, and refines the first downward crossing with Brent's method to 1e-7 s ([Shampine & Thompson (2000)](https://doi.org/10.1016/S0898-1221%2800%2900045-6)) |

Starting on the surface (a launch) is fine: only a downward crossing counts, and the start itself is
never an event.

`propagateToEvent` is general: any function of time and state works, for example height crossings,
apogee (radial velocity changing sign) or node crossings (z changing sign).

## Earth orientation data matter

The impact point is found in the Earth-fixed frame, so the Earth's orientation at the impact time
matters:

| Left out | Error in the impact point |
|---|---|
| UT1 - UTC (up to 0.9 s) | up to ~420 m east-west |
| polar motion (up to ~0.5 arcsec) | up to ~15 m |
| IAU 2000B instead of 2000A (the default, ~1 mas) | ~3 cm |

Load the IERS EOP (Earth Orientation Parameters) with `EopTable::loadFinals2000A`
(see "External data" in [Getting started](getting-started.md)). `ConstantEop{}` (all zero) is only for
examples and tests.

## Many predictions: share the Earth rotations

Each call builds an `ItrsRotationSeries` over `[epoch, epoch + maxFlightTime]` unless you pass one.
For many predictions in the same time window, build it once:

<!-- snippet: examples/impact.cpp#shared -->
```cpp
// Many predictions in the same time window: build the Earth rotations once and share them.
const ItrsRotationSeries rotations(epoch, 3 * 3600.0, eop, ImpactOptions{}.frameOptions);
ImpactOptions            options;
options.rotations = &rotations;
for (const double speed : {4000.0, 4500.0, 5000.0}) {
    const Vec3        v = vEnu * (speed / V);
    const StateVector s = itrsToGcrs(StateVector{geodeticToEcef(site), enu.transpose() * v}, epoch, eop);
    if (const auto hit = predictImpact(epoch, s, Frame::GCRS, eop, options)) {
        std::cout << std::setprecision(0) << speed << std::setprecision(3) << " m/s -> latitude "
                  << rad2deg(hit->point.lat) << ", longitude " << rad2deg(hit->point.lon) << " deg\n";
    }
}
```

| Shot | Flight | Per call | Shared rotations | Shared, with covariance |
|---|---|---|---|---|
| 3 km/s, 35 deg | 413 s | 25 us | 9 us | 14 us |
| 5 km/s, 35 deg | 936 s | 29 us | 14 us | 26 us |
| 7 km/s, 35 deg | 2353 s | 38 us | 21 us | 50 us |

From `bench/impact_bench.cpp` ([Performance](../performance.md)). A series and the propagators are
immutable, so threads can share them.

## Limits

- **No drag yet.** The flight is ballistic in vacuum. Below about 100 km altitude drag usually
  dominates, so for re-entry the impact point can be off by kilometres. A drag `ForceModel` can be
  passed through `ImpactOptions::additionalForces`; a built-in atmosphere model is planned.
- **The ellipsoid is not the ground.** Mean sea level (the geoid) differs from the WGS84 ellipsoid
  by up to about 100 m, and terrain adds more. Set `ImpactOptions::surfaceHeight` to the local
  ground height above the ellipsoid. A height error dh moves the impact point by about
  dh / tan(impact angle) horizontally (170 m for 100 m at 30 degrees).
- **Linear uncertainty**, as described above.
- **Event sampling.** A crossing that goes in and out of the surface within a quarter of an
  integrator step could be missed. For an impact this would need a trajectory that grazes the
  surface and climbs back out.
