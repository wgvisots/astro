# Propagation

## Choosing a propagator

| Propagator | Model | Accuracy (LEO) | Speed (4 h at 10 s) | Use it for |
|---|---|---|---|---|
| `propagateKepler` | point mass | km per orbit (no J2) | ~150 ns per call | quick two-body estimates, Lambert transfers |
| `SecularPropagator` | J2, J2^2, J4 mean-element drift (Brouwer) | **tens to hundreds of km** from a state vector; ~2-6 km with fitted mean elements | 40 ns per output | mean-element studies: node and perigee drift, sun-synchronous and repeat-track design |
| `SecularPropagatorCorrected` | the same positions; velocity = d(position)/dt | the same | ~2 % more | the same, when velocity must agree with position (range rate, Doppler, closest approach) |
| `NumericalPropagator` | point mass + zonal harmonics (+ third bodies, custom forces) | set by tolerance: 0.2 m over 4 h at relTol 1e-9, 0.2 mm at the default | 55 us per propagation at relTol 1e-9 | accurate J2/J4 trajectories, covariance |
| `Sgp4` | SGP4/SDP4 from two-line elements | that of the TLE (typically ~km) | 0.3 us per call (1.3 us deep space) | catalogue element sets |

Numbers are from [Accuracy](../accuracy.md) and [Performance](../performance.md).

## Numerical J2/J4 on an output grid

For many equally spaced outputs, use `propagateGrid(start, step, count)`. The integrator takes its own
error-controlled steps (RKF7(8) [Fehlberg (1968)](https://ntrs.nasa.gov/citations/19680027281), about 90-210 s in LEO) instead
of stopping at every output time, and the outputs in between come from a degree-7 Hermite
polynomial that matches position, velocity, acceleration and jerk at both ends of each step
[Hairer et al. (1993), sec. II.6](https://doi.org/10.1007/978-3-540-78862-1). The cost therefore barely depends on the output spacing,
and the interpolation adds no measurable error.

<!-- snippet: examples/propagation.cpp#numerical -->
```cpp
// Numerical J2 + J4, 4 h at 10 s. Say which frame the state is in: the zonal field is applied
// about the true pole (in GCRS ~0.15 deg off z). relTol 1e-9 keeps LEO within ~0.2 m over 4 h.
AdaptiveOptions opt;
opt.relTol                          = 1e-9;
opt.absTol                          = 1e-6;
const NumericalPropagator      prop = NumericalPropagator::j2j4(epoch, s0, Frame::GCRS, opt);
const std::vector<StateVector> grid = prop.propagateGrid(0.0, 10.0, 1441); // t = 0, 10, ..., 14400 s
std::cout << "numerical, 4 h: " << grid.back().r << " m\n";
```

**Tell it the frame.** J2 acts about the Earth's true pole, which is ~0.15 deg from GCRS z; the
`Frame` argument aligns the field (see [Reference frames](frames.md)). The overload without a frame
assumes z already is the true pole (TOD, TEME, CIRS).

**Pick the tolerance from your error budget.** Worst position error over LEO, SSO and an eccentric
orbit, against a tight-tolerance reference:

| `relTol` / `absTol` | Steps per 4 h | Error after 4 h | Error after 1 day |
|---|---|---|---|
| 1e-12 / 1e-9 (default) | 158 | 0.2 mm | 3 mm |
| 1e-10 / 1e-7 | 90 | 2 cm | 0.4 m |
| 1e-9 / 1e-6 | 68 | 0.2 m | 5 m |
| 1e-8 / 1e-5 | 52 | 2.3 m | 60 m |

These are integration errors against the same J2/J4 model; the model itself omits drag, solar
pressure and higher harmonics, which usually dominate.

`GravityField` takes any list of zonal coefficients J2...Jn (default J2, J3, J4) and a rotation to the
pole frame; `ThirdBodyForce` adds the Sun or Moon from any `Ephemeris`; derive from `ForceModel` for
your own forces. `propagate(dt)` and `propagate(dts)` stop exactly at the requested times instead.

## Two-body

<!-- snippet: examples/propagation.cpp#two_body -->
```cpp
const StateVector kepler = propagateKepler(s0, 14400.0); // universal variables, point mass only
std::cout << "two-body differs from J2/J4 by " << norm(kepler.r - grid.back().r) / 1e3 << " km\n";
```

`propagateKepler` uses universal variables with the Laguerre-Conway iteration
[Conway (1986)](https://doi.org/10.1007/BF01230852), valid for elliptic, parabolic and hyperbolic orbits. `lambertIzzo`
solves Lambert's problem, including multi-revolution branches [Izzo (2015)](https://doi.org/10.1007/s10569-014-9587-y). Classical
and equinoctial elements [Broucke & Cefola (1972)](https://doi.org/10.1007/BF01228432) convert to and from states in
`orbit/elements.h`.

## Analytic J2/J4 (secular), and its limits

`SecularPropagator` advances mean elements with the J2, J2^2 and J4 secular rates of Brouwer's theory
[Brouwer (1959)](https://doi.org/10.1086/107958). It is set up once per orbit and then costs about 40 ns per output.

<!-- snippet: examples/propagation.cpp#secular -->
```cpp
// Analytic J2/J4 mean-element drift. Fast (~40 ns per output) but the state's osculating
// elements are taken as mean elements, so expect kilometres of error: see the propagation guide.
const SecularPropagator        secular(s0); // z must be the pole: fine for this illustration only
const std::vector<StateVector> sec = secular.propagateGrid(0.0, 10.0, 1441);
std::cout << "secular differs from numerical by " << norm(sec.back().r - grid.back().r) / 1e3 << " km\n";
```

**It needs mean elements, and a state vector gives osculating ones.** J2 makes the osculating
semi-major axis oscillate by kilometres twice per orbit, so treating a state's elements as mean
gets the period wrong and the error grows every orbit:

| Orbit | Worst error over 4 h, elements from a state | Fitted mean elements (2 orbits of states) |
|---|---|---|
| LEO 400 km | 117 km | 1.9 km |
| SSO 800 km, e = 0.01 | 136 km | 5.9 km |
| 1000 x 20000 km | 219 km | 3.8 km |

Even fitted, the residual stays at the 1-4 km level: those are the short-period oscillations a
mean-element theory leaves out. For metre-level work use `NumericalPropagator`. The state's frame
must have z along the pole (TOD, TEME, CIRS), or be rotated with `poleOfDateMatrix` first.

### Two velocities: `SecularPropagator` or `SecularPropagatorCorrected`

The two give **identical positions** and differ only in the velocity:

| | `SecularPropagator` | `SecularPropagatorCorrected` |
|---|---|---|
| Velocity | two-body velocity of the mean elements at each time (as Vallado's PKEPLER) | time derivative of the position: v = (Mdot / n) v_mean + raanDot (z x r) + argpDot (w x r) |
| Agrees with dr/dt of its own positions | no: off by 1.2-9.0 m/s (the drift of node, perigee and mean anomaly, ~1e-3 relative) | yes: to 0.002 m/s (the finite-difference check's own error) |
| Velocity error against numerical J2/J4, fitted mean elements, 4 h | 6.3 / 8.3 / 3.9 m/s (LEO / SSO / HEO) | 3.6 / 8.4 / 4.0 m/s |
| Reproduces the input state at dt = 0 | yes, position and velocity | position only; the velocity differs by the drift term (~6-10 m/s in LEO) |
| A propagated state converts back to the mean elements | yes, exactly | no |
| Cost | 39 ns per grid output | ~2 % more |

**Neither velocity is the true one.** Against the numerical J2/J4 trajectory both miss by metres
per second, set by the short-period terms a mean-element theory leaves out, and the corrected one
is better only in LEO. What the correction buys is *consistency* between velocity and position.

Choose **`SecularPropagatorCorrected`** when velocity and position are used together:

- range rate and Doppler, or any rate of a geometric quantity;
- closest-approach and event searches that find a zero of a rate (`Trajectory` uses it for the
  secular models: with the original velocity, nearly co-orbital objects closing at a few m/s had
  their time of closest approach off by up to 90 s);
- interpolating outputs with Hermite polynomials, or differencing or integrating them.

Choose **`SecularPropagator`** when the state stands for the mean orbit:

- mean-element studies, where `stateToElements` of an output should give the propagated mean
  elements back;
- starting from a state that must be reproduced exactly at the epoch, such as a filter's estimate;
- matching other tools that build the state from propagated elements (PKEPLER and its relatives).

Both take the same constructors and `SecularPropagatorCorrected` is a `SecularPropagator`, so it
can be passed anywhere the original is accepted. `velocityIsPositionRate()` says which one an
instance is. The free functions `propagateSecular` keep the original velocity.

Starting the corrected propagator from mean elements chosen so that its velocity *does* match an
input state at the epoch was tried and rejected: the mean semi-major axis it implies is worse. The
4 h position error grew in five of six test orbits (LEO 117 -> 195 km, polar 143 -> 571 km,
near-equatorial 22 -> 896 km) and improved only for the eccentric 1000 x 20000 km orbit
(219 -> 154 km); see [Accuracy](../accuracy.md#analytic-j2j4-secular).

Mission-design helpers: `secularRates`, `sunSynchronousInclination`, `repeatGroundTrackSemiMajorAxis`,
`kCriticalInclination`.

## 12-component states: position, velocity, acceleration, jerk

`StateVector12` carries acceleration and jerk too. Under gravity alone, position and velocity
determine the trajectory, and the model determines acceleration and jerk; `HigherOrderMode` says what
input acceleration and jerk that disagree with the model mean:

- `ModelOnly` (default): propagate from position and velocity; the differences are reported
  (`accelerationResidual`, `jerkResidual`) but ignored.
- `UnmodeledForce`: the differences act as an extra force da0 + dj0 t (thrust, drag) during the
  propagation, so the propagated state reproduces the input exactly at t = 0.

<!-- snippet: examples/propagation.cpp#state12 -->
```cpp
// A 12-component state: the input acceleration and jerk differ from gravity by a small thrust,
// which UnmodeledForce keeps acting (da0 + dj0 t) during the propagation.
StateVector12 s12     = prop.withDerivatives(0.0, s0); // r, v and the model's a and j
s12.a                += normalized(s0.v) * 1e-3;       // plus 1 mm/s^2 along-track
const auto thrusting  = NumericalPropagator::j2j4(epoch, s12, Frame::GCRS, HigherOrderMode::UnmodeledForce, opt);
std::cout << std::setprecision(6) << norm(thrusting.accelerationResidual()) << " m/s^2 unmodeled; after 1 h "
          << std::setprecision(1) << norm(thrusting.propagate(3600.0).r - prop.propagate(3600.0).r) / 1e3
          << " km from the unthrusted orbit\n";
std::cout << "mode: " << HigherOrderMode::UnmodeledForce << '\n'; // enums print as Enum::Value
```

`propagateKepler`, `propagateSecular` and `NumericalPropagator` all accept `StateVector12`; the
outputs carry the model's acceleration and jerk at the new time (the jerk is exact for gravity).

## SGP4

`Sgp4` implements SGP4/SDP4 as revised by Vallado et al. [Vallado et al. (2006)](https://celestrak.org/publications/AIAA/2006-6753/)
(after Spacetrack Report #3 [Hoots & Roehrich (1980)](https://celestrak.org/NORAD/documentation/spacetrk.pdf)) and reproduces their
verification set to 0.12 mm. It outputs TEME; see [Satellite geometry](satellite-geometry.md).

## Closest approach of two objects

`closestApproach(first, second, start, end)` finds the smallest distance between two objects over a
time window. `closeApproaches(..., maxDistance)` lists every local minimum of the distance under a
threshold: a conjunction list. Each object is a `Trajectory`, built from any of these:

| Built from | Propagated with |
|---|---|
| `{epoch, StateVector}` | `TrajectoryModel::TwoBody` (default, `propagateKepler`), `SecularJ2` or `SecularJ2J4` (`SecularPropagator`) |
| `{epoch, KeplerianElements}` | the same; osculating elements for two-body, mean elements for the secular models |
| `{epoch, EquinoctialElements}` | the same, via `equinoctialToKeplerian` |
| any `std::function<StateVector(const Time&)>` | your own: a `NumericalPropagator`, an `Sgp4`, an ephemeris table, ... |

So two element sets, two states, or one of each pass straight in:

<!-- snippet: examples/propagation.cpp#close_approach -->
```cpp
// Closest approach within a day of a state vector and an element set (either may be either, or
// any function of time, e.g. [&prop](const Time& t) { return prop.propagate(t); }). Both must be
// in the same inertial frame; the default model is two-body.
const KeplerianElements debris{6780e3, 0.001, deg2rad(97.6), 0.25, 0.0, 1.22};
const Time              end = epoch + 86400.0;
const CloseApproach     ca  = closestApproach({epoch, s0}, {epoch, debris}, epoch, end);
std::cout << "closest approach " << ca.distance / 1e3 << " km at " << ca.time.iso8601(0) << ", "
          << ca.relativeSpeed / 1e3 << " km/s\n";
// Every approach within 10 km: a conjunction list.
for (const CloseApproach& c : closeApproaches({epoch, s0}, {epoch, debris}, epoch, end, 10e3)) {
    std::cout << "  " << c.time.iso8601(0) << "  " << c.distance / 1e3 << " km\n";
}
```

**How it searches.** The range rate g(t) = (r2 - r1) . (v2 - v1) is negative while the objects
close and positive while they separate. It is sampled with a step of 0.05 |r| / |v| (the shorter of
the two objects' values at each sample: ~125 samples per circular orbit, 45 s in LEO, 11 min in
GEO), so eccentric orbits are sampled finely near perigee and coarsely near apogee. Each change of
sign from negative to positive brackets a minimum, which Brent's method refines to
`timeTolerance` (1 us). The minima and maxima of the distance between two orbits about the same
body are typically a large fraction of a revolution apart, so the step leaves a wide margin: over 300 random orbit
pairs (LEO to 60 000 km, e up to 0.9) a 1-day search found every minimum a 0.5 s brute-force scan
did, at the default step and at twice it. `closestApproach` also compares the window's two ends
(flagged `atBoundary`) for objects still closing at the end or already separating at the start.

**Both objects must be in the same inertial frame.** Convert first if they are not (for example
an `Sgp4` gives TEME; see [Reference frames](frames.md)). The secular models also need the Earth's
pole as z (TOD, TEME or CIRS), as for `SecularPropagator`.

**Secular velocity.** The secular models use `SecularPropagatorCorrected`, whose velocity is the
time derivative of its position (see [Two velocities](#two-velocities-secularpropagator-or-secularpropagatorcorrected)).
The search finds zeros of the range rate, so that matters: with `SecularPropagator`'s velocity,
nearly co-orbital objects (relative speeds of m/s) had their time of closest approach moved by
tens of seconds.

**Cost.** One propagation of each object per sample: about 2000 samples per day of window for
two LEO objects, plus ~10 per minimum found. A wrapped `NumericalPropagator::propagate(t)`
integrates from its epoch on every call, so keep its epoch near the window, or sample it once
with `propagateGrid` and interpolate in your function.

The result is geometric: it does not include the objects' sizes, covariances or a collision
probability.
