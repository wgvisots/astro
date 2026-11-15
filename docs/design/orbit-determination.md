# Orbit determination (design)

| | |
|---|---|
| Status | Planned for after 1.0 (1.1 or later). Nothing here is implemented yet. |
| Agreed | 2026-11-15 |
| Module | a new one, in this repository |

This page records what the feature must do, the design agreed so far, the questions still open, and
the order of work. It changes as the work proceeds; when a phase is done, the user guide takes over
its description.

## What it must do

1. **Fit an orbit to position and velocity reports.** The reports are scattered and irregular in
   time, but certainly of the same object. The result is the orbit (a state at an epoch, and its
   orbital elements) and an error estimate (a covariance).
2. **Determine an orbit from lines of sight.** Each observation gives a time, the observer's position
   and velocity, and a line of sight to the satellite. The result is the same as in 1.

Both must work for **low Earth orbit and geostationary orbit**, with observations that may be
**several orbits apart**.

## Why in this repository

Most of what orbit determination needs is already here and public:

| Needed | Already in astro |
|---|---|
| A propagator, and the state transition matrix (how a later state depends on the initial one) | `NumericalPropagator::propagateStm`, `propagateGridStm`, `propagateKeplerStm`, `StateTransition` |
| Covariance propagation, element conversions, consistency statistics | `propagateWithCovariance`, `Mat6`, `mahalanobis2`, `EstimateTracker` |
| Frames for the observer and the satellite | `FrameTransformer`, `EopTable` |
| Light time and aberration for lines of sight | the astrometry module |
| A two-point boundary solver, for initial orbits from lines of sight | `lambertIzzo` |

It needs no new dependency: the systems to solve are 6×6, or a little larger with estimated
parameters.

## Before 1.0: keep the state transition matrix public

The fit's core is the state transition matrix. The pre-1.0 review of the public API
([inventory and decisions](api-review-1.0.md), done in 0.4.0) kept these public:

- `NumericalPropagator::propagateStm` and `propagateGridStm`
- `propagateKeplerStm`
- `StateTransition` (`astro/orbit/state_covariance.h`)
- `ForceModel::partials`, which the new force models will implement (the review changed it to
  return a `ForcePartials`)

Then orbit determination can arrive in a 1.x release as a pure addition, as the
[compatibility policy](../compatibility.md) requires. Each of these declarations carries a comment
that points here.

## Inputs

### Position and velocity reports

Each report is labelled **independent** or **don't know**.

| Label | Meaning | How it is used |
|---|---|---|
| Independent | a measurement whose error is independent of every other report's | a measurement in the least-squares fit |
| Don't know, with a covariance | possibly the output of a filter, so correlated with earlier reports | **covariance intersection** ([Julier & Uhlmann 1997](../references.md)) with the fit, after a consistency check against it; inconsistent reports are flagged, not merged |
| Don't know, without a covariance | cannot be weighted | a starting guess or a sanity check only |

A filter's output already contains its earlier inputs, so fusing such outputs as independent would
give a covariance that is too small: the reason `EstimateTracker` does not fuse them
([Covariance and tracking](../guide/covariance-tracking.md)). Covariance intersection stays
consistent for any correlation, at the cost of a larger covariance. The result says which reports
were used, and how.

### Lines of sight

The provider's documentation is vague, so the inputs are described here as received:

- the position and velocity of the centre of the sensor's focal plane, in "ECI";
- a unit vector in "ECI" from the centre of the focal plane to the object;
- possibly, rates in azimuth and elevation (rad/s).

Unknowns, and their size:

| Unknown | Effect if wrong | Handling |
|---|---|---|
| Which "ECI": J2000, GCRS, TEME, true of date | TEME and true of date differ from GCRS by about 0.3 degrees today; J2000 by milliarcseconds | an option |
| Aberration: none, from the observer's orbital velocity, or also annual (if the direction was calibrated against stars) | about 5 arcseconds for an observer in low orbit; up to about 20 arcseconds annual | an option |
| Light time | about 2 arcseconds for a geostationary object | an option, on by default |
| Time tag: start, middle or end of the exposure | 750 m along the track for 0.1 s on a low-orbit object | an estimated offset |
| Angle noise | unknown | estimated from the residuals of the fit |
| Frame of the azimuth and elevation rates | unknown for a sensor in orbit | not used at first; once known, they help the initial orbit from short tracks |

These are settled by calibration (phase D), not by guessing.

## Force models needed

A least-squares fit spans gaps of several orbits only if the force model holds over them.

| Force | Why | In astro today |
|---|---|---|
| Gravity field to degree and order *n* ([EGM2008, Pavlis et al. 2012](../references.md)) | tesseral terms: J22 drives the longitude drift of geostationary orbits; metre-level effects in low orbit within a few revolutions | zonal terms only |
| Atmospheric drag | dominant in low orbit | none |
| Solar radiation pressure, with the Earth's shadow | important in geostationary orbit | the shadow geometry only (`eclipse.h`) |
| Sun and Moon | both regimes | `ThirdBodyForce` |
| Estimated parameters: drag and radiation-pressure coefficients | area-to-mass ratios are unknown | the transition matrix is 6×6; a sensitivity to parameters is needed, as an addition |

Atmosphere models, in order:

1. **Harris-Priester** (as in Montenbruck & Gill 2000, section 3.5): simple, no space-weather input.
2. **NRLMSIS 2.0** (Emmert et al. 2021), with a loader for the space-weather file (F10.7 and Ap, from
   CelesTrak), kept up to date by the application like the EOP file.
3. **JB2008** (Bowman et al. 2008), more accurate in low orbit, if its licence allows.

Each is written from its paper, and checked against the reference code's published outputs. The
licence of each is checked before work starts.

## Methods

### Phase B: least-squares fit to position and velocity reports

A batch weighted least-squares fit (Tapley, Schutz & Born 2004, chapter 4) solves for one state at a
chosen epoch, plus the estimated parameters:

1. Start from any one report.
2. For each report, predict it from the current estimate, and take the residual. For a position and
   velocity report, the dependence of the prediction on the epoch state is the state transition
   matrix itself.
3. Solve the normal equations, weighted by each report's covariance, and update the estimate.
4. Repeat until the update is small; reject reports whose residuals are too large (with a stated
   threshold, and a list of what was rejected).

The result:

- the epoch state and its covariance, the inverse of the normal matrix;
- the covariance scaled by the weighted residuals: if the reports' covariances are right, the
  weighted root-mean-square residual is near 1;
- the residual of each report, the rejected reports, and the number of iterations;
- the osculating elements at the epoch, and their covariance;
- a warning when the normal matrix is badly conditioned (too little geometry in the data).

### Phase C: lines of sight

1. **Initial orbit from three lines of sight, with no prior.** Gooding's method (Gooding 1997): robust,
   handles a moving observer, and solves through a Lambert problem (`lambertIzzo`). Gauss's method
   (Vallado 2013, section 7.3) as the simple fallback. With observations orbits apart, the number of
   revolutions between them is ambiguous, so the initial orbit tries several and the fit chooses.
2. **The fit of phase B with an angle measurement:** two angles per observation, with light time and
   the chosen aberration, in the chosen frame.

From a short arc the range is poorly observable, so the covariance is long along the line of sight;
an observer that moves adds triangulation. The fit reports weak geometry rather than return a
confident wrong answer.

### Phase D: calibration against a known orbit

The sensor observes a satellite whose orbit is published to the centimetre level: a GPS satellite
(IGS precise orbits) or a laser-ranging satellite (ILRS). Fitting the unknowns listed under
*Lines of sight* (in *Inputs*, above) against that truth shows which frame, which aberration and which
time tag the provider uses. This is a study in `studies/`, like the others in
[Accuracy and validation](../accuracy.md), so the conventions are measured and written down.

## Validation

As for the rest of the library ([Accuracy and validation](../accuracy.md)):

- **Simulated data:** observations from known orbits, with known noise. The truth must fall inside
  the reported covariance at the expected rate (the chi-square test of `EstimateTracker`), over many
  trials.
- **Published examples:** worked examples of the initial-orbit methods.
- **Force models:** each against its reference outputs.
- **Real data:** the calibration study of phase D.

## Order of work

| Phase | Content | Size |
|---|---|---|
| A. Forces | gravity field of degree and order *n*; drag (Harris-Priester, then NRLMSIS and space weather); radiation pressure with shadow; sensitivity to estimated parameters | 3 to 4 pull requests; useful for propagation on its own |
| B. Position and velocity fit | least squares; covariance intersection for "don't know" reports; outlier rejection; consistency statistics; element covariance | 1 to 2 pull requests |
| C. Lines of sight | measurement model and its options; Gooding and Gauss with revolution hypotheses; angles in the fit | 2 pull requests |
| D. Calibration | the conventions of the lines of sight, against a precise orbit | 1 pull request, once real data exists |

## Open questions

- The provider's conventions for the lines of sight (see *Lines of sight*, under *Inputs*), until
  phase D answers them.
- The frame of the azimuth and elevation rates.
- Whether "don't know" reports without a covariance arrive often enough to need more than a sanity
  check.
- The licences of NRLMSIS 2.0 and JB2008.

## Sources

Already in [References](../references.md): Julier & Uhlmann (1997), Pavlis et al. (2012), Montenbruck
& Gill (2000), Vallado (2013), Izzo (2015). New for this work, to be added there when used:

- **Tapley, B. D., Schutz, B. E. & Born, G. H. (2004).** *Statistical Orbit Determination*. Elsevier
  Academic Press. Batch least squares, the state transition matrix in estimation.
- **Gooding, R. H. (1997).** A new procedure for the solution of the classical problem of minimal
  orbit determination from three lines of sight. *Celest. Mech. Dyn. Astr.* 66, 387-423.
- **Emmert, J. T. et al. (2021).** NRLMSIS 2.0: A whole-atmosphere empirical model of temperature
  and neutral species densities. *Earth and Space Science* 8, e2020EA001321.
  <https://doi.org/10.1029/2020EA001321>.
- **Bowman, B. R. et al. (2008).** A new empirical thermospheric density model JB2008 using new solar
  and geomagnetic indices. AIAA 2008-6438.
