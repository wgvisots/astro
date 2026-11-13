# Propagation

## Choosing a propagator

| Propagator | Model | Accuracy (LEO) | Speed (4 h at 10 s) | Use it for |
|---|---|---|---|---|
| `propagateKepler` | point mass | km per orbit (no J2) | ~150 ns per call | quick two-body estimates, Lambert transfers |
| `SecularPropagator` | J2, J2^2, J4 mean-element drift (Brouwer) | **tens to hundreds of km** from a state vector; ~2-6 km with fitted mean elements | 40 ns per output | mean-element studies: node and perigee drift, sun-synchronous and repeat-track design |
| `SecularPropagatorCorrected` | the same positions; velocity = d(position)/dt | the same | ~2 % more | the same, when velocity must agree with position (range rate, Doppler, closest approach) |
| `NumericalPropagator` | point mass + zonal harmonics (+ third bodies, custom forces) | set by tolerance: 0.2 m over 4 h at relTol 1e-9, 0.2 mm at the default | 55 us per propagation at relTol 1e-9 | accurate J2/J4 trajectories, covariance |
| `Sgp4` | SGP4/SDP4 from two-line elements | that of the TLE (typically ~km) | 0.3 us per call (1.3 us deep space) | catalogue element sets |

Numbers are from [Accuracy](../accuracy.md) and [Performance](../performance.md). All of them except
`Sgp4` also propagate a covariance: see [Covariance and tracking](covariance-tracking.md).

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

**A grid may lie anywhere relative to the epoch.** It can start before the epoch, span it, or step
towards it. The integration only ever moves away from the epoch, so every output has the accuracy of
integrating straight to its own time, and an output at the epoch is the initial state exactly.

**A force that depends on time alone needs a step limit.** The step-size control compares
evaluations at the two ends of a step [Fehlberg (1968)](https://ntrs.nasa.gov/citations/19680027281).
A part of the force that depends only on time, such as a thrust profile, cancels in that comparison
exactly, so the integrator cannot tell when it is stepping over it.

- If such a force varies smoothly, set `AdaptiveOptions::maxStep` to a fraction of the time over
  which it changes.
- If it jumps (an engine switching on or off), propagate to the jump and start a new
  `NumericalPropagator` from the state there.
- Forces that depend on position or velocity (gravity, drag, third bodies) are not affected.

## Two-body

<!-- snippet: examples/propagation.cpp#two_body -->
```cpp
const StateVector kepler = propagateKepler(s0, 14400.0); // universal variables, point mass only
std::cout << "two-body differs from J2/J4 by " << norm(kepler.r - grid.back().r) / 1e3 << " km\n";
```

`propagateKepler` uses universal variables with the Laguerre-Conway iteration
[Conway (1986)](https://doi.org/10.1007/BF01230852), valid for elliptic, parabolic and hyperbolic orbits.
`propagateKeplerStm` adds its state transition matrix in closed form, and
`propagateKeplerWithCovariance` a propagated covariance. Classical
and equinoctial elements [Broucke & Cefola (1972)](https://doi.org/10.1007/BF01228432) convert to and from states in
`orbit/elements.h`.

### Lambert's problem

Lambert's problem asks for the two-body orbit that joins two positions in a given time.

<!-- snippet: examples/propagation.cpp#lambert -->
```cpp
// Lambert's problem: the orbit that joins two positions in a given time. These two are 40 minutes
// apart on the two-body orbit above, so the solution gives back its velocity.
const StateVector     later    = propagateKepler(s0, 2400.0);
const LambertSolution transfer = lambert(s0.r, later.r, 2400.0);
std::cout << "Lambert recovers the velocity to " << norm(transfer.v1 - s0.v) * 1e6 << " um/s\n";
// With 5 hours to spare, other orbits make complete revolutions on the way: two for each number.
const std::vector<LambertSolution> all = lambertIzzo(s0.r, later.r, 18000.0, constants::gmEarth, true, 3);
std::cout << all.size() << " orbits join the same points in 5 h with up to 3 revolutions\n";
```

- `lambert` returns the solution without a complete revolution.
- `lambertIzzo` returns that one first, then two solutions for each number of complete revolutions
  the time allows, up to `maxRevolutions`.
- `prograde` picks the direction of motion: counter-clockwise seen from +z, or clockwise.
- The method is Izzo's [Izzo (2015)](https://doi.org/10.1007/s10569-014-9587-y), written from the paper: a
  third-order iteration on the variable of Lancaster and Blanchard
  [Lancaster & Blanchard (1969)](https://ntrs.nasa.gov/citations/19690027552). One problem takes about 0.3 us
  ([Performance](../performance.md)).
- The arrival point is reproduced to about 1e-12 of its radius. Three kinds of problem are
  sensitive by nature, whatever the solver: positions within about 1e-6 rad of opposite, flights
  of very many time units on one arc, and the long way round in a very short time
  ([Accuracy](../accuracy.md)).
- Positions that are exactly collinear with the centre leave the orbit plane undefined, and throw.

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
| Reproduces the input state at dt = 0 | yes, position and velocity | position only; the velocity differs by the drift term (5-20 m/s in LEO: 5 at 51.6 deg, 9 sun-synchronous, 20 near-equatorial) |
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
| `{epoch, StateVector}` | `TrajectoryModel::TwoBody` (default, `propagateKepler`), `SecularJ2` or `SecularJ2J4` (`SecularPropagatorCorrected`) |
| `{epoch, KeplerianElements}` | the same; osculating elements for two-body, mean elements for the secular models |
| `{epoch, EquinoctialElements}` | the same, via `equinoctialToKeplerian` |
| any callable `StateVector(const Time&)`, e.g. a lambda | your own: a `NumericalPropagator`, an `Sgp4`, an ephemeris table, ... |

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
sign from negative to positive brackets a minimum, which Brent's method
([Brent 1973](../references.md#orbits-and-propagation)) refines to `timeTolerance` (1 us). Searching
for the zeros of the range rate is the approach of [Alfano & Negron (1993)](../references.md#conjunctions-and-collision-probability), who
fit cubic polynomials through its samples; here each bracketed zero is refined on the trajectories
themselves. Analytic
pre-filters such as [Hoots, Crawford & Roehrich (1984)](https://doi.org/10.1007/BF01234152) can
thin a large catalogue before it. The minima and maxima of the distance between two orbits about
the same body are typically a large fraction of a revolution apart, so the step leaves a wide
margin: over 300 random orbit pairs (LEO to 66 700 km, e up to 0.9) a 1-day search found every
minimum a 0.5 s brute-force scan did, at the default step and at twice it (see
[Accuracy](../accuracy.md#closest-approach-search)). `closestApproach` also compares the window's two ends
(flagged `atBoundary`) for objects still closing at the end or already separating at the start.

**Both objects must be in the same inertial frame.** Convert first if they are not (for example
an `Sgp4` gives TEME; see [Reference frames](frames.md)). The secular models also need the Earth's
pole as z (TOD, TEME or CIRS), as for `SecularPropagator`.

**Secular velocity.** The secular models use `SecularPropagatorCorrected`, whose velocity is the
time derivative of its position (see [Two velocities](propagation.md#two-velocities-secularpropagator-or-secularpropagatorcorrected)).
The search finds zeros of the range rate, so that matters: with `SecularPropagator`'s velocity,
nearly co-orbital objects (relative speeds of m/s) had their time of closest approach moved by
tens of seconds.

**Cost.** One propagation of each object per sample: about 2000 samples per day of window for
two LEO objects, plus ~10 per minimum found. A wrapped `NumericalPropagator::propagate(t)`
integrates from its epoch on every call, so keep its epoch near the window, or sample it once
with `propagateGrid` and interpolate in your function.

### Collision probability

`collisionProbability(approach, covariance1, covariance2, hardBodyRadius)` turns an approach into
the probability that the objects actually collide, given their sizes and position uncertainties:

<!-- snippet: examples/propagation.cpp#collision_probability -->
```cpp
// Probability of collision at that approach, for a 5 m and a 1 m object. Position covariances at
// the time of the approach: 400 m radial, 4 km in-track, 400 m cross-track (1 sigma) in each
// object's RIC frame, rotated into the inertial frame. (From an epoch covariance,
// NumericalPropagator::propagateGridWithCovariance gives them at any time.)
// clang-format off
const Mat3 ric{400.0 * 400.0,             0.0,           0.0,
                         0.0, 4000.0 * 4000.0,           0.0,
                         0.0,             0.0, 400.0 * 400.0};
// clang-format on
const Mat3                 m1 = ricMatrix(ca.first), m2 = ricMatrix(ca.second);
const CollisionProbability pc = collisionProbability(ca, m1.transpose() * ric * m1, m2.transpose() * ric * m2, 6.0);
std::cout << "collision probability " << std::scientific << std::setprecision(2) << pc.probability << std::fixed
          << std::setprecision(1) << ", miss " << pc.mahalanobis << " sigma\n";
```

It is the standard short-term encounter model of [Foster & Estes (1992)](../references.md#conjunctions-and-collision-probability),
[Akella & Alfriend (2000)](https://doi.org/10.2514/2.4611) and [Alfano (2005a)](https://doi.org/10.1007/BF03546397);
[Chan (2008)](https://doi.org/10.2514/4.989186) gives the background:

| Step | What happens |
|---|---|
| Objects | spheres whose radii add up to the hard-body radius (HBR) |
| Uncertainty | independent Gaussian position errors at TCA: the relative position has covariance C1 + C2 |
| Encounter plane | the plane through the miss vector perpendicular to the relative velocity; errors along the relative velocity only shift the time of the encounter, so they drop out |
| Probability | the 2-D Gaussian of the miss in that plane integrated over the HBR circle (`probabilityInCircle`): exact along the covariance's minor axis, adaptive Gauss-Kronrod along the major axis, ~1e-12 relative in 7-18 us |

The result also gives the miss distance, the projected 2x2 covariance, its sigmas and the miss in
standard deviations (`mahalanobis`).

**Validated against NASA CARA (Conjunction Assessment Risk Analysis) reference cases**
([CARA Analysis Tools](https://github.com/nasa/CARA_Analysis_Tools)): to 1e-12 relative where CARA's own values are that
good, and to an independent integral where they are not (a near-singular covariance, and large radii). Also
against [Alfano's (2009)](../references.md#conjunctions-and-collision-probability) test cases (see [Accuracy](../accuracy.md#collision-probability)).
One convention differs from CARA's `PcCircle`: CARA takes the full distance |r2 - r1| as the miss,
assuming the states are exactly at TCA; `collisionProbability` projects r2 - r1 onto the encounter
plane, which is the closest distance of the straight-line encounter. They agree for a true closest
approach from `closestApproach`; for states 2.1 deg off TCA (CARA's Omitron 2 case) they differ by 0.017 %.

**When it applies.** The model assumes the encounter is short: straight-line relative motion and
constant covariances while the objects cross each other's uncertainty, as at relative speeds of
km/s. Alfano's slow test cases (0.2-0.5 m/s) still agree with his Monte Carlo to within 2 %,
because their relative motion stays nearly straight, but for slow, curving, nearly co-orbital
encounters use a Monte Carlo or a method for nonlinear relative motion instead, such as
[Patera (2003)](https://doi.org/10.2514/2.5127). Use a true closest approach (not
`atBoundary`).

**More uncertainty does not always mean more risk.** For a given miss, the probability peaks when
the covariance is about the size of the miss and falls as it grows (probability dilution;
[Alfano 2005b](https://doi.org/10.1007/BF03546350), [Balch, Martin & Ferson 2019](https://doi.org/10.1098/rspa.2018.0565)). In the
example, growing the covariance from 400 m / 4 km (radial / in-track) to 800 m / 6 km raises Pc from
1.9e-8 to 4.8e-7, but on to 1 km / 8 km barely changes it (4.84e-7): it has levelled off near its
peak. A small Pc from a large covariance can mean the data cannot tell, not that the approach is safe.
Balch, Martin & Ferson show that this is not peculiar to this model: with poor data, any
epistemic probability of collision assigns high confidence to "no collision" whether or not it is
true (false confidence).

### From states and covariances

`assessConjunction` and `assessConjunctions` do all three steps from each object's state and
covariance at its epoch: find the approach, propagate both covariances to it, and compute the
probability. Each object is an `UncertainTrajectory`:

| Built from | Covariance propagated with |
|---|---|
| `{epoch, state, covariance}` | `propagateKeplerWithCovariance` (default `TrajectoryModel::TwoBody`) |
| `{epoch, state, covariance, TrajectoryModel::SecularJ2J4}` (or `SecularJ2`) | `SecularPropagatorCorrected::propagateWithCovariance` |
| a `Trajectory` and any function of time returning a `StateCovariance` | your own: a `NumericalPropagator`'s `propagateGridWithCovariance`, `EstimateTracker::predict`, ... |

<!-- snippet: examples/propagation.cpp#conjunction_assessment -->
```cpp
// The same from each object's state and covariance at the epoch, in one call: the covariances
// are propagated to every approach within 10 km with the trajectory's model (two-body here;
// TrajectoryModel::SecularJ2J4 needs a frame with z along the pole). 1 sigma per axis: 100 m and
// 10 cm/s for the first object, 1 km and 1 m/s for the debris.
Mat6 P1{}, P2{};
for (std::size_t i = 0; i < 3; ++i) {
    P1[i][i]         = 100.0 * 100.0;
    P1[i + 3][i + 3] = 0.1 * 0.1;
    P2[i][i]         = 1000.0 * 1000.0;
    P2[i + 3][i + 3] = 1.0 * 1.0;
}
const UncertainTrajectory first{epoch, s0, P1}, second{epoch, elementsToState(debris), P2};
for (const ConjunctionAssessment& c : assessConjunctions(first, second, epoch, end, 6.0, 10e3)) {
    std::cout << "  " << c.approach.time.iso8601(0) << "  " << c.approach.distance / 1e3 << " km, Pc "
              << std::scientific << std::setprecision(2) << c.collision.probability << std::fixed
              << std::setprecision(1) << " (" << c.collision.mahalanobis << " sigma)\n";
}
```

```text
  2026-03-01T02:00:09Z  6.2 km, Pc 1.48e-14 (5.9 sigma)
  2026-03-01T02:46:27Z  2.4 km, Pc 1.01e-07 (1.5 sigma)
```

Each `ConjunctionAssessment` holds the `CloseApproach`, both covariances at its time, and the
`CollisionProbability`. `assessConjunction` gives the closest approach only; that is not always the
most probable collision, so to rank a window's approaches use `assessConjunctions` and sort by
`collision.probability`.

- **The search samples the trajectory; the covariance is computed once per approach.** The cost is
  that of `closeApproaches` plus 0.3 us (two-body) or 6.4 us (secular) per approach and object.
- **Validated end to end.** In the tests, 100 000 draws of both initial states, each propagated and
  searched, pass within the hard-body radius as often as `assessConjunction` predicts (Pc 0.3,
  checked to 2 %).
- **For a NumericalPropagator**, give the trajectory and the covariance separately, so that the search
  does not integrate the covariance at every sample:

```cpp
const auto prop = NumericalPropagator::j2j4(epoch, s0, Frame::TOD);
const UncertainTrajectory object{
    Trajectory([&prop](const Time& t) { return prop.propagate(t); }),
    [&prop, P0, epoch](const Time& t) { return prop.propagateGridWithCovariance(t - epoch, 1.0, 1, P0)[0]; }};
```

  Each `propagate(t)` integrates from the epoch, so keep the window near it (see the cost note above).

