# Glossary

Acronyms and terms used in this documentation and in the API, grouped by topic. Within each group
they are in alphabetical order.

## Time

| Term | Stands for | What it means here |
|---|---|---|
| BDT | BeiDou Time | BeiDou navigation system time; TAI - 33 s (`TimeScale::BDT`) |
| GPS | Global Positioning System | GPS time is TAI - 19 s, counted in weeks and seconds of week (`TimeScale::GPS`) |
| GST | Galileo System Time | Galileo navigation system time, aligned with GPS time (`TimeScale::GST`). Not Greenwich Sidereal Time, which this library calls GMST or GAST |
| ISO 8601 | International Organization for Standardization standard 8601 | text times such as `2026-03-01T12:00:00Z` (`Time::fromIso8601`, `Time::iso8601`) |
| J2000.0 | Julian epoch 2000.0 | 2000 January 1, 12:00 TT (JD 2451545.0) |
| JD | Julian Date | days since 4713 BC January 1, 12:00; stored in two parts for precision (`JulianDate`) |
| Leap second | | a second inserted into UTC to keep it within 0.9 s of UT1 (`LeapSecondTable`) |
| MJD | Modified Julian Date | JD - 2400000.5; days since 1858 November 17, 0:00 |
| TAI | International Atomic Time (Temps Atomique International) | the continuous atomic time scale; UTC differs from it by whole leap seconds |
| TCB | Barycentric Coordinate Time | relativistic coordinate time of the solar-system barycenter |
| TCG | Geocentric Coordinate Time | relativistic coordinate time of the geocenter |
| TDB | Barycentric Dynamical Time | the time argument of planetary ephemerides; within 2 ms of TT |
| TOW | Time Of Week | seconds since the start of a GNSS week |
| TT | Terrestrial Time | TAI + 32.184 s; the time argument of precession-nutation |
| UT1 | Universal Time 1 | time from the Earth's actual rotation angle; differs from UTC by up to 0.9 s (see EOP) |
| UTC | Coordinated Universal Time | civil time: atomic seconds, kept near UT1 by leap seconds |

## Reference frames and Earth orientation

| Term | Stands for | What it means here |
|---|---|---|
| CIO | Celestial Intermediate Origin | the origin of right ascension on the CIP equator in the IERS 2010 transformation |
| CIP | Celestial Intermediate Pole | the Earth's rotation pole used by the IAU 2006/2000 models |
| CIRS | Celestial Intermediate Reference System | GCRS rotated by precession-nutation: CIP pole, CIO origin (`Frame::CIRS`) |
| dX, dY | | small measured corrections to the modelled CIP position, part of the EOP |
| ECEF | Earth-Centered, Earth-Fixed | rotates with the Earth; in this library the ITRS |
| ECI | Earth-Centered Inertial | does not rotate with the Earth; in this library the GCRS by default, or J2000, TOD, TEME, ... |
| EME2000 | Earth Mean Equator and equinox of J2000 | same as the J2000 frame (`Frame::J2000`) |
| EOP | Earth Orientation Parameters | measured UT1 - UTC, polar motion, length of day and dX, dY; published daily by the IERS (`EopTable`) |
| ERA | Earth Rotation Angle | the angle the Earth has turned, a linear function of UT1 |
| Frame bias | | the small fixed rotation between J2000 (EME2000) and GCRS (~23 mas) |
| GAST | Greenwich Apparent Sidereal Time | the Earth's rotation angle measured from the true equinox |
| GCRS | Geocentric Celestial Reference System | the IAU's geocentric inertial frame; the library's default ECI (`Frame::GCRS`) |
| GMST | Greenwich Mean Sidereal Time | the Earth's rotation angle measured from the mean equinox |
| IAU 2006/2000A, 2000B | | the IAU precession (2006) and nutation (2000A, full; 2000B, truncated to ~1 mas) models (`PrecessionNutationModel`) |
| ICRS / ICRF | International Celestial Reference System / Frame | the barycentric celestial reference, realized by distant quasars; GCRS is its geocentric counterpart |
| ITRS / ITRF | International Terrestrial Reference System / Frame | the Earth-fixed frame; WGS84 agrees with it to a few centimetres (`Frame::ITRS`) |
| LAST / LMST | Local Apparent / Mean Sidereal Time | sidereal time at a given longitude |
| LOD | Length Of Day | excess of the day over 86 400 s, in seconds (part of the EOP) |
| MOD | Mean Of Date | mean equator and equinox at a given time: precession only (`Frame::MOD`) |
| Nutation | | the short-period wobble of the Earth's axis (18.6 years and shorter) |
| PEF | Pseudo Earth-Fixed | equinox-based Earth-fixed frame before polar motion (`Frame::PEF`) |
| Polar motion | | the motion of the rotation axis relative to the Earth's crust (xp, yp; part of the EOP) |
| Precession | | the slow (26 000-year) turning of the Earth's axis |
| TEME | True Equator, Mean Equinox | the frame of SGP4 output (`Frame::TEME`) |
| TIO | Terrestrial Intermediate Origin | the longitude origin on the CIP equator in the Earth-fixed frame |
| TIRS | Terrestrial Intermediate Reference System | Earth-fixed frame before polar motion (`Frame::TIRS`) |
| TOD | True Of Date | true equator and equinox at a given time: precession and nutation (`Frame::TOD`) |

## Geodesy and local frames

| Term | Stands for | What it means here |
|---|---|---|
| AER | Azimuth, Elevation, Range | look angles from a site; azimuth clockwise from north (`ecefToAer`) |
| Ellipsoid height | | height above the reference ellipsoid along its normal (`Geodetic::h`); not height above sea level |
| ENU | East, North, Up | local horizontal frame at a site (`ecefToEnuMatrix`) |
| Geoid | | the equipotential surface that mean sea level follows; up to ~100 m above or below the WGS84 ellipsoid |
| GRS80 | Geodetic Reference System 1980 | reference ellipsoid, within 0.1 mm of WGS84 (`Ellipsoid::grs80`) |
| NED | North, East, Down | local frame common in aviation (`ecefToNedMatrix`) |
| WGS72 / WGS84 | World Geodetic System 1972 / 1984 | WGS84 is the GPS reference ellipsoid and frame (`Ellipsoid::wgs84`); WGS72 constants are used by SGP4 |

## Orbits and propagation

| Term | Stands for | What it means here |
|---|---|---|
| Closest approach | | the smallest distance between two objects, where their range rate is zero (`closestApproach`) |
| Conjunction | | a close approach of two objects, usually below a screening distance (`closeApproaches`) |
| Event | | a zero crossing of a function of time and state, such as height at impact (`propagateToEvent`) |
| Flight-path angle | | angle of the velocity above (launch) or below (impact) the local horizontal |
| GEO | Geostationary Earth Orbit | circular equatorial orbit with a one-sidereal-day period |
| GM, mu | Gravitational constant times Mass | gravitational parameter of a body (m^3/s^2) |
| HEO | Highly Elliptical Orbit | e.g. Molniya |
| J2, J3, J4 | | zonal harmonics of the Earth's gravity field: J2 for the equatorial bulge, J3 for the pear shape, J4 next |
| Lambert's problem | | find the two-body orbit that joins two positions in a given time (`lambertIzzo`) |
| LEO | Low Earth Orbit | altitudes up to about 2000 km |
| NORAD | North American Aerospace Defense Command | its catalogue numbers identify objects in TLEs |
| ODE | Ordinary Differential Equation | the equations of motion the integrators solve |
| RAAN | Right Ascension of the Ascending Node | orientation of the orbit plane |
| RK4 | Runge-Kutta, 4th order | fixed-step integrator (`integrateRk4`) |
| RKF7(8) | Runge-Kutta-Fehlberg 7(8) | adaptive-step integrator used by `NumericalPropagator` |
| SGP4 / SDP4 | Simplified General / Deep-space Perturbations 4 | the analytic propagators for TLEs (`Sgp4`) |
| SSO | Sun-Synchronous Orbit | orbit whose plane turns with the Sun, through J2 |
| TCA | Time of Closest Approach | when the distance between two objects is smallest (`CloseApproach::time`) |
| TLE | Two-Line Element set | the NORAD element format for SGP4 (`Tle`) |

## Covariance, tracking and impact

| Term | Stands for | What it means here |
|---|---|---|
| 1-sigma | | one standard deviation; a 1-sigma error ellipse contains 39.3 % of 2-D Gaussian errors |
| CARA | Conjunction Assessment Risk Analysis | NASA's conjunction assessment team; its open-source tools give the collision probability test cases |
| CEP | Circular Error Probable | radius of the circle with 50 % of impacts. Not computed directly: `ellipse.scaled(0.5)` gives the 50 % ellipse |
| Chi-square | | distribution of a sum of squared Gaussian errors; gives the p-values of the consistency checks |
| Covariance | | the 6x6 matrix of position and velocity error variances and correlations (`Mat6`) |
| Encounter plane | | the plane perpendicular to the relative velocity at closest approach, where the collision probability is integrated (`CollisionProbability::xAxis`, `yAxis`) |
| Error ellipse | | the curve of equal probability on the ground; semi-major, semi-minor and rotation (`ErrorEllipse`) |
| HBR | Hard-Body Radius | the radius of a sphere enclosing both objects: the sum of their radii (`collisionProbability`) |
| LVLH | Local Vertical, Local Horizontal | orbit-following frame, CCSDS convention (`lvlhMatrix`) |
| Mahalanobis distance | | a difference measured in units of its covariance: d^T P^-1 d |
| NTW, VNC | Normal-Tangential-W, Velocity-Normal-Co-normal | orbit-following frames (`ntwMatrix`, `vncMatrix`) |
| Pc | Probability of collision | chance that two objects pass within the hard-body radius (`collisionProbability`) |
| Probability dilution | | a larger covariance can lower Pc for the same miss: a small Pc may only mean the data are poor |
| Process noise | | extra uncertainty added during propagation for unmodelled forces (`accelerationNoisePsd`) |
| PSD | Power Spectral Density | strength of white process noise (m^2/s^3 for acceleration noise) |
| RIC / RSW | Radial, In-track, Cross-track | orbit-following frame for residuals and covariances (`ricMatrix`) |
| Rotation (of an ellipse) | | direction of the semi-major axis, clockwise from north |
| Shear (Keplerian) | | the growth of in-track error because a different energy means a different period; dominates propagated covariances |
| STM | State Transition Matrix | Phi: how a small change of the initial state changes the state later (`propagateKeplerStm`, `propagateStm`, `propagateGridStm`) |
| Symplectic | | Phi^T J Phi = J with J = [[0, I], [-I, 0]]: the property of every Phi of a Hamiltonian (e.g. purely gravitational) motion; a check on a computed Phi |

## Satellite geometry

| Term | Stands for | What it means here |
|---|---|---|
| AOS / LOS | Acquisition / Loss Of Signal | when a satellite rises above / sets below a site's elevation mask |
| Beta angle | | angle between the orbit plane and the Sun direction (`betaAngle`) |
| LOS | Line Of Sight | also used for an unobstructed straight line between two points (`lineOfSight`) |
| Umbra / penumbra | | full / partial Earth shadow (`shadow`, `findEclipses`) |

## Data, standards and software

| Term | Stands for | What it means here |
|---|---|---|
| C04 | IERS EOP C04 series (14 C04, 20 C04) | the IERS combined series of final EOP values (`EopTable::loadC04`) |
| DAF | Double precision Array File | NAIF's binary file format underlying SPK files |
| DCM | Direction Cosine Matrix | a rotation matrix |
| DE440, DE441, DE440s | Development Ephemeris 440 (441, short 440s) | JPL planetary ephemerides, read by `SpkEphemeris` |
| EGM2008 | Earth Gravitational Model 2008 | source of the J2-J6 values and the reference radius |
| ERFA | Essential Routines for Fundamental Astronomy | BSD-licensed copy of SOFA; the reference values in the tests come from it |
| IAU | International Astronomical Union | defines the time scales and precession-nutation models |
| IERS | International Earth Rotation and Reference Systems Service | publishes the EOP, leap seconds and the IERS Conventions (2010) |
| JPL | Jet Propulsion Laboratory | NASA laboratory that publishes the DE ephemerides |
| NAIF | Navigation and Ancillary Information Facility | JPL group behind the SPICE toolkit and SPK format |
| SI | International System of Units | metres, seconds, kilograms; all inputs and outputs use SI units and radians |
| SOFA | Standards Of Fundamental Astronomy | the IAU's reference software for these algorithms |
| SPK | Spacecraft and Planet Kernel | NAIF ephemeris file format (`.bsp`) |
| SSB | Solar System Barycenter | center of mass of the solar system; the root of ephemeris chains |
