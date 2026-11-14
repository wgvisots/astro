# Covariance and tracking

## Propagating a covariance

`propagateGridWithCovariance` carries a 6x6 position/velocity covariance along the output grid. The
state transition matrix Phi is integrated alongside the state from the force models' exact partial
derivatives (the variational equations [Montenbruck & Gill (2000), ch. 7](https://doi.org/10.1007/978-3-642-58351-3)), and

    P(t) = Phi(t) P0 Phi(t)^T + Q(t)

<!-- snippet: examples/covariance.cpp#covariance -->
```cpp
// 10 m and 1 cm/s (1-sigma) per axis at the epoch, in the state's frame (GCRS).
Mat6 P0{};
for (std::size_t i = 0; i < 3; ++i) {
    P0[i][i]         = 10.0 * 10.0;
    P0[i + 3][i + 3] = 0.01 * 0.01;
}
AdaptiveOptions opt;
opt.relTol      = 1e-9;
opt.absTol      = 1e-6;
const auto prop = NumericalPropagator::j2j4(epoch, s0, Frame::GCRS, opt);
// States and covariances every 10 min for 4 h, with a little process noise (m^2/s^3).
const std::vector<StateCovariance> out = prop.propagateGridWithCovariance(0.0, 600.0, 25, P0, 1e-10);
```

- The integrator's error control looks at the state only, so the steps, and the states, are exactly
  those of `propagateGrid`. Between steps Phi uses quintic Hermite interpolation. Against stopping
  exactly at each output, over 4 h at 10 s: 1e-8 relative at the default tolerance and 1e-6 at
  relTol 1e-9 for near-circular orbits, ten times that for a 1000 x 20000 km orbit.
- `Q(t)` is optional process noise: white, isotropic acceleration noise with power spectral density
  q (m^2/s^3), standing in for what the model leaves out. Over a short interval it adds about
  q t^3 / 3 to each position variance.
- `propagateGridStm` returns Phi itself.
- Cost: about 6x the states alone (8x with process noise); 0.3 ms for 4 h at 10 s at relTol 1e-9.

Uncertainty is easiest to read in radial / in-track / cross-track axes:

<!-- snippet: examples/covariance.cpp#ric -->
```cpp
std::cout << " t (h)   sigma R (m)   sigma I (m)   sigma C (m)\n" << std::fixed << std::setprecision(1);
for (std::size_t k = 0; k < out.size(); k += 6) {
    const Mat6 ric = rotateCovariance(out[k].covariance, ricMatrix(out[k].state));
    std::cout << std::setw(5) << k / 6.0 << std::setw(14) << std::sqrt(ric[0][0]) << std::setw(14)
              << std::sqrt(ric[1][1]) << std::setw(14) << std::sqrt(ric[2][2]) << '\n';
}
```

In-track uncertainty grows fastest: a velocity error changes the period.

### With the analytic propagators

The two-body and secular propagators carry a covariance too, without integrating anything:

| Propagator | Call | How Phi is found | Cost |
|---|---|---|---|
| two-body | `propagateKeplerWithCovariance(s0, P0, dt)`; `propagateKeplerStm` for Phi | closed form: the universal-variable Lagrange coefficients differentiated (the matrix that [Shepperd (1985)](https://doi.org/10.1007/BF01227666) and [Battin (1999)](https://doi.org/10.2514/4.861543) derive in other forms); exact to rounding, any conic | 0.31 us |
| secular J2/J4 | `SecularPropagator::propagateWithCovariance(dt, P0)`, `propagateGridWithCovariance`; `propagateStm` for Phi | fourth-order central differences of the whole state-to-state map; 1e-11 to 2e-5 relative | 6.4 us per time, 1.2 us per grid output |
| numerical | `propagateGridWithCovariance` (above) | variational equations, integrated | 0.3 ms for 4 h at 10 s |

<!-- snippet: examples/covariance.cpp#analytic_covariance -->
```cpp
// The same initial covariance through the analytic propagators, at 4 h: two-body (closed-form
// state transition matrix) and secular J2/J4. The secular model needs z along the pole, so the
// state and covariance are rotated from GCRS first (a fixed rotation, which commutes with the
// propagation; RIC sigmas do not depend on the frame).
const StateCovariance            kepler = propagateKeplerWithCovariance(s0, P0, 4.0 * 3600.0);
const Mat3                       toPole = poleOfDateMatrix(Frame::GCRS, epoch);
const SecularPropagatorCorrected secular(StateVector{toPole * s0.r, toPole * s0.v});
const StateCovariance            mean = secular.propagateWithCovariance(4.0 * 3600.0, rotateCovariance(P0, toPole));
const auto                       row  = [](const char* name, const StateCovariance& s) {
    const Mat6 ric = rotateCovariance(s.covariance, ricMatrix(s.state));
    std::cout << std::left << std::setw(26) << name << std::right << std::setw(9) << std::sqrt(ric[0][0])
              << std::setw(14) << std::sqrt(ric[1][1]) << std::setw(14) << std::sqrt(ric[2][2]) << '\n';
};
std::cout << "at 4 h                     sigma R (m)   sigma I (m)   sigma C (m)\n";
row("numerical J2/J4 (+ noise)", out.back());
row("secular J2/J4", mean);
row("two-body", kepler);
```

```text
at 4 h                     sigma R (m)   sigma I (m)   sigma C (m)
numerical J2/J4 (+ noise)      44.5         675.0           9.7
secular J2/J4                  44.0         675.5           9.7
two-body                       44.0         675.1           9.7
```

- **The covariance is mostly Keplerian shear.** Against the numerical J2/J4 covariance, both analytic
  ones match the in-track sigma to 0.6 % at 4 h on all six test orbits, and two-body does as well as
  secular. Over a day the small radial sigma differs by up to 37 % (it oscillates, and the analytic
  trajectories are at a different phase), and the HEO in-track sigma by 10 %
  ([Accuracy](../accuracy.md#analytic-state-transition-matrices-and-covariances)). The state, not the
  covariance, is where the models really differ: tens of km for the secular model, more for two-body.
- **No process noise.** The analytic covariances show only how the model maps the initial
  uncertainty. For process noise, or forces beyond J2/J4, use `NumericalPropagator`.
- **For a whole grid, the numerical propagator is faster.** Differencing costs 25 secular grids, 1.7 ms
  for 4 h at 10 s, against 0.3 ms for the numerical covariance at relTol 1e-9, which is also more
  accurate. The secular covariance pays off at single times (6.4 us) and where it must go with the
  secular trajectory; the two-body one is cheap anywhere.
- **Which state.** P0 is the covariance of the state the propagator was built from
  (`epochState()`; for mean elements, their state). `SecularPropagatorCorrected`'s velocity at the
  epoch is not the input velocity, so its Phi(0) is not the identity.
- **Frames.** As for the states: any inertial frame for two-body, z along the pole for the secular
  models (rotate with `poleOfDateMatrix` and `rotateCovariance`, as above).
- `assessConjunction` uses these to propagate both objects' covariances to their closest approach
  (see [Propagation](propagation.md#from-states-and-covariances)).

## Following an upstream filter

If you receive a stream of states and covariances from one upstream filter (an orbit determination,
a tracking system), each delivery already contains everything before it. Successive deliveries are
correlated, so **do not fuse them** as if they were independent measurements: a Kalman update would
count the same data again and again, and the covariance would collapse far below the truth. (To
combine estimates from different sources with unknown correlation, covariance intersection is
the standard tool [Julier & Uhlmann (1997)](https://doi.org/10.1109/ACC.1997.609105).)

`EstimateTracker` does what is valid: the newest delivery replaces the current estimate ("latest
wins"), and each delivery is first checked against the previous one propagated (with its covariance)
to the new epoch.

<!-- snippet: examples/tracker.cpp#tracker -->
```cpp
TrackerOptions options;
options.frame = Frame::GCRS;      // frame of the deliveries (also TOD, TEME, J2000, ...)
EstimateTracker tracker(options); // J2 + J4, relTol 1e-9, element consistency test

(void)tracker.update({t0, s0, P}); // the first delivery has nothing to be checked against
if (auto check = tracker.update({t0 + 1800.0, truth(t0, s0, 1800.0), P})) {
    std::cout << "delivery 2 (no burn):\n";
    print(*check);
}
```

A delivery after a burn:

<!-- snippet: examples/tracker.cpp#flagged -->
```cpp
if (auto check = tracker.update({t0 + 3600.0, truth(t0 + 3300.0, burn, 300.0), P}); check && !check->consistent) {
    std::cout << "delivery 3 (after a burn):\n";
    print(*check);
}
// Latest wins: predictions now start from delivery 3.
const auto ahead = tracker.predictGrid(0.0, 60.0, 61); // states + covariances for the next hour
std::cout << "predicted 1 h ahead: " << ahead.back().state.r << " m\n";
```

## The consistency check

`ConsistencyCheck` reports:

| Field | Meaning |
|---|---|
| `positionResidualRic`, `velocityResidualRic` | new minus predicted, in radial / in-track / cross-track (m, m/s) |
| `positionSigmaRic`, `velocitySigmaRic` | the predicted 1-sigma in the same axes |
| `semiMajorAxisResidual`, `semiMajorAxisSigma` | change of the osculating semi-major axis, a clear sign of a burn |
| `elementsMahalanobis2`, `elementsPValue` | 6-D distance in equinoctial elements, chi-square p-value |
| `positionMahalanobis2`, `positionPValue` | 3-D distance in position |
| `mahalanobis2`, `pValue` | 6-D distance in Cartesian position and velocity |
| `consistent` | the selected test's p-value is at least `alarmProbability` (default 1e-3) |

The distances use the predicted covariance, so for a consistent upstream filter the p-values are
conservative [Bar-Shalom et al. (2001)](https://doi.org/10.1002/0471221279). `TrackerOptions::test` selects the statistic:

| `ConsistencyTest` | False alarms (nominal 1e-3), gaps to 4 h | Sees a burn just before the delivery |
|---|---|---|
| `Elements` (default) | nominal in every case (1 or 2 alarms in 1500 samples) | yes |
| `Position` | nominal (0 to 5 alarms in 1500 samples) | no: blind until the burn has moved the position |
| `PositionVelocity` (Cartesian) | up to 0.06 (20 m / 2 cm/s priors) and 0.34 (100 m / 0.1 m/s) at 4 h | yes |

Linear covariance propagation in Cartesian coordinates stops describing the error once it stretches
along the orbit (the "banana" shape [Junkins et al. (1996)](../references.md)); in orbital elements the
error grows mostly in one coordinate, the mean longitude, and stays close to Gaussian much longer
[Aristoff, Horwood & Alfriend (2021)](https://doi.org/10.1007/s10569-021-10004-0); the Monte Carlo comparison here follows
[Sabol et al. (2010)](../references.md). The full Monte Carlo tables are in [Accuracy](../accuracy.md).

## Frames, tolerance and process noise

- `TrackerOptions::frame` is the frame of the deliveries (GCRS by default; J2000, MOD, TOD, TEME or
  CIRS). Of-date frames are taken at each delivery's own epoch, and predictions are expressed in the
  frame of their own time, matching how later deliveries will arrive (to 1.4 mm over 4 h).
- `TrackerOptions::body` is the gravity model: `CentralBody::earthJ2J4()` (point mass, J2 and J4) by
  default; any `CentralBody` may be given (see [Propagation](propagation.md)).
- `TrackerOptions::integrator` defaults to relTol 1e-9 (0.2 m over 4 h in LEO).
- `TrackerOptions::accelerationNoisePsd` defaults to 0. With drag acting in low LEO, good deliveries
  will start to fail over long gaps; raise it until the in-track residuals of routine deliveries sit
  within a few predicted sigmas and the alarm rate is near `alarmProbability`.
- Const member functions (`current`, `predict`, `predictGrid`) may be called from several threads;
  `update` needs exclusive access.
