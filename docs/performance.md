# Performance

All timings: one core of an Intel Core i7-1260P (laptop, power-saving governor), GCC 13, `-O3`,
Release build. Laptop timings vary by 10-30% from run to run; compare ratios rather than single
numbers, and measure on your own hardware with the programs in [`bench/`](https://github.com/wgvisots/astro/blob/main/bench).

## Propagating one orbit to many times

LEO, output every Y seconds:

| Method | X = 1 day, Y = 60 s | X = 1 day, Y = 10 s | X = 4 h, Y = 10 s |
|---|---|---|---|
| `propagateSecular(state, dt)` called per output | 201 ns/output | 227 ns/output | 201 ns/output |
| `SecularPropagator::propagateGrid` | 40 ns/output | 43 ns/output | 40 ns/output |
| `SecularPropagatorCorrected::propagateGrid` | 39 ns/output | 39 ns/output | 39 ns/output |
| `SecularPropagator::propagateGrid12` | 69 ns/output | 68 ns/output | 68 ns/output |
| `NumericalPropagator::propagate(dts)`, stopping at each output | 0.68 ms | 4.4 ms | 0.68 ms |
| `NumericalPropagator::propagateGrid`, default tolerance | 0.51 ms | 0.56 ms | 0.10 ms |
| `NumericalPropagator::propagateGrid`, relTol 1e-9 | | | 0.056 ms |
| `NumericalPropagator::propagateGrid12`, relTol 1e-9 | | | 0.13 ms |

`SecularPropagatorCorrected` was timed in a later run, alongside `SecularPropagator::propagateGrid` at
38.7 ns/output: the velocity correction adds about 2 %. Times for the numerical propagator are per propagation (all outputs). Two-body `propagateKepler`
costs about 150 ns per call; SGP4 about 0.3 us per call (1.3 us for deep-space SDP4).

What makes the grids fast:

- The numerical grid lets the integrator take its natural steps (about 96 s in LEO at the default
  tolerance) and interpolates the outputs, so going from Y = 60 s to Y = 10 s costs almost nothing;
  stopping at each output time would cap the step at Y.
- `SecularPropagator` does the per-orbit work once, and steps the angles of an evenly spaced grid
  with angle-addition formulas and a warm-started Kepler solve, so there are no library sin/cos
  calls per output (recomputed exactly every 64 outputs).
- Gravity evaluations do not allocate, share one pass for point mass and zonal terms, and avoid
  `std::hypot`.

## Covariance, checks and frames

4 h at Y = 10 s (1441 outputs), LEO:

| relTol | States only | With covariance | With covariance and process noise |
|---|---|---|---|
| 1e-12 | 111 us | 533 us | 751 us |
| 1e-9 | 54 us | 308 us | 426 us |

| Operation | Time |
|---|---|
| `propagateKeplerWithCovariance` (closed-form Phi) | 0.31 us |
| `lambertIzzo`, no complete revolution | 0.27 us |
| `lambertIzzo`, up to 3 revolutions (all solutions of one problem) | 1.3 us |
| `SecularPropagatorCorrected::propagateWithCovariance` (24 differenced propagations) | 6.4 us |
| `SecularPropagatorCorrected::propagateGridWithCovariance`, 4 h at 10 s (1441 outputs) | 1.7 ms (1.2 us per output) |
| `checkConsistency` / `EstimateTracker::update`, gap 10 min / 30 min / 4 h | 15 / 26 / 137 us |
| `FrameTransformer`, IAU 2006/2000A | 125 us |
| `FrameTransformer`, IAU 2006/2000B | 6 us |
| `poleOfDateMatrix(Frame::GCRS, t)` | 6 us |

## Impact prediction

From `bench/impact_bench.cpp`: `predictImpact` for launches from the ellipsoid (35 degrees above
the horizon), with the Earth rotations built per call or shared between calls:

| Shot | Flight | Per call | Shared rotations | Shared, with covariance |
|---|---|---|---|---|
| 3 km/s | 413 s | 32 us | 12 us | 18 us |
| 5 km/s | 936 s | 31 us | 15 us | 29 us |
| 7 km/s | 2353 s | 41 us | 24 us | 50 us |

| Operation | Time |
|---|---|
| `ItrsRotationSeries`, 2 h window (3 nodes), IAU 2006/2000B | 17 us |
| `ItrsRotationSeries`, 2 h window (3 nodes), IAU 2006/2000A | 400 us |

## Getting the most out of it

- Choose the integrator tolerance from your error budget: relTol 1e-9 (0.2 m over 4 h in LEO) is
  twice as fast as the default; see [Accuracy](accuracy.md).
- Use the grid calls for equally spaced outputs, and the span overloads (`propagateGrid(start, step,
  std::span<StateVector>)`) to reuse your output buffers.
- Construct a `FrameTransformer` once per epoch, not per vector; for many times in a window, use
  an `ItrsRotationSeries` (GCRS <-> ITRS) and share it, e.g. through `ImpactOptions::rotations`.
- Propagate many orbits in parallel: the propagators are immutable, so threads can share them, and
  separate orbits scale with the number of cores.

## Closest approach and collision probability

From `astro_conjunction_bench`, 1-day window:

| Call | Time |
|---|---|
| `closestApproach`, two LEO objects, two-body | 0.66 ms (about 2000 samples) |
| `closestApproach`, two LEO objects, secular J2/J4 | 0.40 ms |
| `closestApproach`, LEO and GTO, two-body | 0.72 ms |
| `closeApproaches` (every minimum), two LEO objects, two-body | 0.65 ms |
| `probabilityInCircle`, typical (miss ~ sigma, radius << sigma) | 6.7 us |
| `probabilityInCircle`, 1000:1 covariance | 12.6 us |
| `probabilityInCircle`, radius >> sigma | 17.5 us |
| `probabilityInCircle`, 123 km x 19 m covariance, 30 km radius | 9-12 us |
| `collisionProbability`, 200 m miss (projection + `probabilityInCircle`) | 6.8 us |
| `assessConjunction`, two LEO objects, two-body (search, covariances, probability) | 0.66 ms |
| `assessConjunction`, two LEO objects, secular J2/J4 | 0.42 ms |

Near-singular covariances at the circle's edge take 10-20 us in the cases of the study and up to
0.12 ms in a random search; the refinement is capped at ~1 ms whatever the input.

Over the 300 random orbit pairs of `astro_study_conjunctions` (LEO to 66 700 km, e up to 0.9), a
1-day search averages 0.21 ms at the default step: higher orbits need fewer samples.

## Reproducing

```sh
cmake -B build -DCMAKE_BUILD_TYPE=Release -DASTRO_BUILD_BENCHMARKS=ON
cmake --build build -j
build/bench/astro_propagation_bench 86400 60 100   # X seconds, Y seconds, number of states
build/bench/astro_propagation_bench 14400 10 100
build/bench/astro_covariance_bench
build/bench/astro_conjunction_bench                # optional argument: window seconds
```
