# Reference frames

The celestial ("ECI") and terrestrial ("ECEF") frames, and the transformations between them,
follow the IERS Conventions (2010) [[ref:iers2010|IERS Conventions (2010), ch. 5]] and match SOFA/ERFA
[[ref:sofa|SOFA]] [[ref:erfa|ERFA]] to about a microarcsecond ([Accuracy](../accuracy.md)).

| `Frame` | What it is |
|---|---|
| `GCRS` | Geocentric Celestial Reference System: the modern "ECI" |
| `J2000` | mean equator and equinox of J2000.0 (EME2000); GCRS plus a 20 mas frame bias |
| `MOD` | mean equator and equinox of date (precession only) |
| `TOD` | true equator and equinox of date (precession and nutation) |
| `TEME` | true equator, mean equinox: the frame of SGP4 and two-line elements [[ref:vallado2006|Vallado et al. (2006)]] |
| `CIRS` | Celestial Intermediate Reference System (true equator, CIO origin) |
| `TIRS` | Terrestrial Intermediate Reference System (Earth-fixed, before polar motion) |
| `PEF` | pseudo Earth-fixed (equinox-based counterpart of TIRS) |
| `ITRS` | International Terrestrial Reference System: the modern "ECEF" |

The main chain is the CIO-based transformation

    ITRS = W^T(xp, yp, s') * R3(ERA) * C(X + dX, Y + dY, s) * GCRS

with the IAU 2006 precession [[ref:hilton2006|Hilton et al. (2006)]] and IAU 2000A nutation
[[ref:capitaine2006|Capitaine & Wallace (2006)]], ERA the Earth rotation angle, and W the polar
motion. The equinox-based frames (MOD, TOD, PEF) and TEME are available for interoperability.
`PrecessionNutationModel::IAU2006_2000B` [[ref:mccarthy2003|McCarthy & Luzum (2003)]] is about 20x
faster (in `FrameTransformer`) and accurate to about 1 mas.

## Earth orientation data

Earth-fixed transformations need UT1 - UTC, polar motion (xp, yp) and, for full accuracy, the
celestial pole offsets (dX, dY) [[ref:iers-eop|IERS EOP]]. Pass an `EopProvider`:

<!-- snippet: examples/frames.cpp#eop -->
```cpp
// Earth orientation: an IERS finals2000A file for full accuracy, or ConstantEop{} (UT1 = UTC, no
// polar motion: good to ~1 s of Earth rotation, ~400 m at the equator).
std::unique_ptr<EopProvider> eop = std::make_unique<ConstantEop>();
if (argc > 1) {
    eop = std::make_unique<EopTable>(EopTable::loadFinals2000A(argv[1]));
}
```

`EopTable` reads `finals2000A` (with predictions about a year ahead) and the IERS 20 C04 series,
and interpolates with 4-point Lagrange polynomials across leap seconds. `ConstantEop{}` (all
zero) is good to about 400 m at the equator: UT1 - UTC can reach 0.9 s.

## ECEF <-> ECI

<!-- snippet: examples/frames.cpp#ecef_eci -->
```cpp
std::cout << std::fixed << std::setprecision(3);
const Time        t = Time::fromIso8601("2024-03-20T12:00:00Z");
const StateVector ecef{geodeticToEcef({deg2rad(38.8895), deg2rad(-77.0353), 100.0}), {}};
const StateVector eci  = itrsToGcrs(ecef, t, *eop); // includes the omega x r velocity
const StateVector back = gcrsToItrs(eci, t, *eop);
std::cout << "GCRS r " << eci.r << " m, v " << eci.v << " m/s\n";
std::cout << "round trip error " << std::scientific << norm(back.r - ecef.r) << std::fixed << " m\n";
```

State transformations include the Earth-rotation term (omega x r) in the velocity, with omega
including the length-of-day variation.

For many vectors at one epoch, build a `FrameTransformer` once: its constructor evaluates the
precession-nutation series (about 125 us with IAU 2000A, 6 us with 2000B), after which each
transformation is a matrix product:

<!-- snippet: examples/frames.cpp#transformer -->
```cpp
// Many vectors at one epoch: build the rotations once.
const FrameTransformer ft(t, *eop);
const Vec3             teme = ft.transform(eci.r, Frame::GCRS, Frame::TEME);
const Vec3             tod  = ft.transform(eci.r, Frame::GCRS, Frame::TOD);
std::cout << "TEME " << teme << "\nTOD  " << tod << '\n';
```

## The true pole, and why it matters for propagation

Zonal gravity (J2, J4) acts about the Earth's true spin axis, the Celestial Intermediate Pole. In
GCRS and J2000 that axis is about 0.13-0.16 deg from z during the 2020s, drifting with precession;
applying J2 about GCRS z instead would cost about 56 m after 1 h and 270 m after 4 h in LEO
([Accuracy](../accuracy.md)). `poleOfDateMatrix` gives the rotation to a frame whose z is the true pole:

<!-- snippet: examples/frames.cpp#pole -->
```cpp
// The Earth's true pole is not GCRS z (0.134 deg apart in 2024, growing ~0.0056 deg a year with
// precession). Zonal gravity must be applied about the true pole.
const Mat3 toPole = poleOfDateMatrix(Frame::GCRS, t);
std::cout << std::setprecision(7) << "true pole in GCRS: " << Vec3{toPole(2, 0), toPole(2, 1), toPole(2, 2)}
          << ", tilt " << std::setprecision(4) << rad2deg(std::acos(toPole(2, 2))) << " deg\n";
```

The propagators take the frame of their state (`NumericalPropagator::j2j4(epoch, state, Frame::GCRS)`,
`TrackerOptions::frame`) and apply it for you. TOD, TEME and CIRS already have the true pole as z.

## Pitfalls

- "ECI" is ambiguous. GCRS and J2000 differ by a 20 mas frame bias (0.7 m at LEO radius, 4 m at
  GEO); TEME and TOD differ from GCRS by about 0.13-0.16 deg in the 2020s (15-19 km at LEO radius).
  Know which one your data uses.
- Of-date frames (MOD, TOD, TEME, CIRS) turn slowly: a fixed vector expressed in TOD of one epoch
  and TOD of 4 h later differs by 0.85 m at LEO radius. `EstimateTracker` expresses predictions in
  the frame of their own time; `NumericalPropagator` treats its frame as fixed at the epoch.
- SGP4 outputs TEME. Convert with `temeToGcrs` / `temeToItrs` before mixing it with other data.
