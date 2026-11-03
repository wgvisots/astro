# Covariance and tracking

## Propagating a covariance

`propagateGridWithCovariance` carries a 6x6 position/velocity covariance along the output grid. The
state transition matrix Phi is integrated alongside the state from the force models' exact partial
derivatives (the variational equations [[ref:montenbruck|Montenbruck & Gill (2000), ch. 7]]), and

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
  those of `propagateGrid`. Between steps Phi uses quintic Hermite interpolation (relative error
  2e-9 at the default tolerance, 4e-7 at relTol 1e-9).
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

## Following an upstream filter

If you receive a stream of states and covariances from one upstream filter (an orbit determination,
a tracking system), each delivery already contains everything before it. Successive deliveries are
correlated, so **do not fuse them** as if they were independent measurements: a Kalman update would
count the same data again and again, and the covariance would collapse far below the truth. (To
combine estimates from different sources with unknown correlation, covariance intersection is
the standard tool [[ref:julier1997|Julier & Uhlmann (1997)]].)

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
conservative [[ref:barshalom|Bar-Shalom et al. (2001)]]. `TrackerOptions::test` selects the statistic:

| `ConsistencyTest` | False alarms (nominal 1e-3), gaps to 4 h | Sees a burn just before the delivery |
|---|---|---|
| `Elements` (default) | 0.0007-0.0013 in every case | yes |
| `Position` | 0-0.003 | no: blind until the burn has moved the position |
| `PositionVelocity` (Cartesian) | up to 0.06 (20 m / 2 cm/s priors) and 0.34 (100 m / 0.1 m/s) at 4 h | yes |

Linear covariance propagation in Cartesian coordinates stops describing the error once it stretches
along the orbit (the "banana" shape [[ref:junkins1996|Junkins et al. (1996)]]); in orbital elements the
error grows mostly in one coordinate, the mean longitude, and stays close to Gaussian much longer
[[ref:sabol2010|Sabol et al. (2010)]]. The full Monte Carlo tables are in [Accuracy](../accuracy.md).

## Frames, tolerance and process noise

- `TrackerOptions::frame` is the frame of the deliveries (GCRS by default; J2000, MOD, TOD, TEME or
  CIRS). Of-date frames are taken at each delivery's own epoch, and predictions are expressed in the
  frame of their own time, matching how later deliveries will arrive (to 1.4 mm over 4 h).
- `TrackerOptions::integrator` defaults to relTol 1e-9 (0.2 m over 4 h in LEO).
- `TrackerOptions::accelerationNoisePsd` defaults to 0. With drag acting in low LEO, good deliveries
  will start to fail over long gaps; raise it until the in-track residuals of routine deliveries sit
  within a few predicted sigmas and the alarm rate is near `alarmProbability`.
- Const member functions (`current`, `predict`, `predictGrid`) may be called from several threads;
  `update` needs exclusive access.
