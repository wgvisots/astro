# References

The standards, papers and books the implementations follow, grouped by topic, with what each is
used for. DOIs and links were checked when this list was written (November 2026); the two marked
"link not checked" could not be reached from the checking machine. Third-party tables and test data
included in the repository are credited in
[THIRD_PARTY_NOTICES.md](https://github.com/wgvisots/astro/blob/main/THIRD_PARTY_NOTICES.md).

## Standards and reference software

- **Petit, G. & Luzum, B. (eds.) (2010).** *IERS Conventions (2010)*. IERS Technical Note 36,
  Verlag des BKG. <https://www.iers.org/TN36> (updates: <https://iers-conventions.obspm.fr/>).
  Time scales, CIO-based GCRS <-> ITRS chain, X/Y/s and nutation series, Helmert transformations,
  EOP interpolation.
- **IAU SOFA Board.** *Standards of Fundamental Astronomy*. <https://www.iausofa.org/>. Cookbooks:
  [SOFA Time Scale and Calendar Tools](https://www.iausofa.org/s/sofa_ts_c.pdf) (UTC quasi-JD
  convention, time scales) and [SOFA Tools for Earth Attitude](https://www.iausofa.org/s/sofa_pn_c.pdf)
  (precession-nutation, frames).
- **ERFA**, *Essential Routines for Fundamental Astronomy* (a BSD-licensed SOFA derivative),
  <https://github.com/liberfa/erfa>, and its Python wrapper **pyerfa**,
  <https://github.com/liberfa/pyerfa>. The reference values in the tests.
- **Urban, S. E. & Seidelmann, P. K. (eds.) (2013).** *Explanatory Supplement to the Astronomical
  Almanac*, 3rd ed. University Science Books. <https://aa.usno.navy.mil/publications/exp_supp>.
  The project's primary reference; rotation conventions, time scales, frames, astrometry.
- **Kaplan, G. H. (2005).** *The IAU Resolutions on Astronomical Reference Systems, Time Scales, and
  Earth Rotation Models*. USNO Circular 179. <https://aa.usno.navy.mil/publications/Circular_179>.

## Time

- **Fairhead, L. & Bretagnon, P. (1990).** An analytical formula for the time transformation
  TB-TT. *A&A* 229, 240-247. <https://ui.adsabs.harvard.edu/abs/1990A%26A...229..240F>. TDB - TT.
- **IERS Bulletin C** (leap seconds): <https://datacenter.iers.org/data/latestVersion/bulletinC.txt>;
  `Leap_Second.dat`: <https://hpiers.obspm.fr/iers/bul/bulc/Leap_Second.dat>.
- **IANA time zone database**, `leap-seconds.list`: <https://data.iana.org/time-zones/data/leap-seconds.list>.
- **IS-GPS-200N** (2022), *NAVSTAR GPS Space Segment / Navigation User Segment Interfaces*.
  <https://www.gps.gov/interface-control-documents-icds-interface-specifications-iss>. GPS time and weeks.
- **European Union (2025).** *Galileo Open Service Signal-In-Space ICD*, Issue 2.2.
  <https://www.gsc-europa.eu/electronic-library/programme-reference-documents/galileo-in-force/open-service>. Galileo System Time.
- **China Satellite Navigation Office (2019).** *BeiDou SIS ICD, Open Service Signal B1I*,
  Version 3.0. <http://en.beidou.gov.cn/SYSTEMS/ICD/> (link not checked). BeiDou time.

## Reference frames and Earth orientation

- **Capitaine, N. & Wallace, P. T. (2006).** High precision methods for locating the celestial
  intermediate pole and origin. *A&A* 450, 855-872. <https://doi.org/10.1051/0004-6361:20054550>.
- **Hilton, J. L. et al. (2006).** Report of the IAU Division I Working Group on Precession and
  the Ecliptic. *Celest. Mech. Dyn. Astr.* 94, 351-367. <https://doi.org/10.1007/s10569-006-0001-2>.
  IAU 2006 precession, Fukushima-Williams angles.
- **McCarthy, D. D. & Luzum, B. J. (2003).** An abridged model of the precession-nutation of the
  celestial pole. *Celest. Mech. Dyn. Astr.* 85, 37-49. <https://doi.org/10.1023/A:1021762727016>. IAU 2000B.
- **Seago, J. & Vallado, D. (2000).** Coordinate frames of the U.S. Space Object Catalogs. AIAA
  2000-4025. <https://doi.org/10.2514/6.2000-4025>; and Vallado et al. (2006), Appendix C. TEME.
- **IERS Earth orientation data**: <https://www.iers.org/iers/en/dataproducts/earthorientationdata/eop>;
  EOP 20 C04: <https://hpiers.obspm.fr/iers/eop/eopc04/>; `finals2000A`:
  <https://maia.usno.navy.mil/ser7/finals2000A.all> (format:
  <https://maia.usno.navy.mil/ser7/readme.finals2000A>).

## Geodesy

- **Vermeille, H. (2011).** An analytical method to transform geocentric into geodetic
  coordinates. *J. Geodesy* 85, 105-117. <https://doi.org/10.1007/s00190-010-0419-x>.
- **Vincenty, T. (1975).** Direct and inverse solutions of geodesics on the ellipsoid with
  application of nested equations. *Survey Review* 23(176), 88-93. <https://doi.org/10.1179/sre.1975.23.176.88>.
- **NGA (2014).** *Department of Defense World Geodetic System 1984: Its Definition and
  Relationships with Local Geodetic Systems*. NGA.STND.0036_1.0.0_WGS84.
  <https://nsgreg.nga.mil/doc/view?i=4085> (link not checked). WGS84 constants.
- **Pavlis, N. K. et al. (2012).** The development and evaluation of the Earth Gravitational Model
  2008 (EGM2008). *J. Geophys. Res.* 117, B04406. <https://doi.org/10.1029/2011JB008916>. Zonal harmonics.

## Ephemerides and astrometry

- **Park, R. S. et al. (2021).** The JPL Planetary and Lunar Ephemerides DE440 and DE441.
  *Astron. J.* 161, 105. <https://doi.org/10.3847/1538-3881/abd414>.
- **NAIF.** *SPK Required Reading* <https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html>
  and *DAF Required Reading* <https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/daf.html>. The SPK reader.
- **Meeus, J. (1998).** *Astronomical Algorithms*, 2nd ed. Willmann-Bell, ISBN 0-943396-61-1.
  Low-precision Sun, rise/set; also the source used for Saemundsson's refraction formula
  (Saemundsson, *Sky & Telescope* 72, 70, 1986).
- **Bennett, G. G. (1982).** The calculation of astronomical refraction in marine navigation.
  *J. Navigation* 35, 255-259. <https://doi.org/10.1017/S0373463300022037>.
- **ESA (1997).** *The Hipparcos and Tycho Catalogues*, ESA SP-1200, Vol. 1, Sect. 1.5.3.
  <https://ui.adsabs.harvard.edu/abs/1997ESASP1200.....E>. ICRS <-> galactic matrix.

## Orbits and propagation

- **Battin, R. H. (1999).** *An Introduction to the Mathematics and Methods of Astrodynamics*,
  revised ed. AIAA. <https://doi.org/10.2514/4.861543>. Universal variables and the two-body state
  transition matrix.
- **Shepperd, S. W. (1985).** Universal Keplerian state transition matrix. *Celestial Mechanics*
  35(2), 129-144. <https://doi.org/10.1007/BF01227666>. The closed-form two-body Phi
  (`propagateKeplerStm`).
- **Montenbruck, O. & Gill, E. (2000).** *Satellite Orbits: Models, Methods and Applications*.
  Springer. <https://doi.org/10.1007/978-3-642-58351-3>. Low-precision Moon, shadow model,
  variational equations and state transition matrix.
- **Vallado, D. A. (2013).** *Fundamentals of Astrodynamics and Applications*, 4th ed. Microcosm
  Press (5th ed. 2022). <https://celestrak.org/software/vallado-sw.php>. Worked examples in the tests.
- **Curtis, H. D. (2020).** *Orbital Mechanics for Engineering Students*, 4th ed. Butterworth-Heinemann.
  <https://shop.elsevier.com/books/orbital-mechanics-for-engineering-students/curtis/978-0-12-824025-0>.
  Lambert example in the tests.
- **Vallado, D. A., Crawford, P., Hujsak, R. & Kelso, T. S. (2006).** Revisiting Spacetrack Report
  #3. AIAA 2006-6753. <https://doi.org/10.2514/6.2006-6753>; paper, code and the SGP4-VER test set:
  <https://celestrak.org/publications/AIAA/2006-6753/>. SGP4/SDP4.
- **Hoots, F. R. & Roehrich, R. L. (1980).** *Models for Propagation of NORAD Element Sets*.
  Spacetrack Report No. 3. <https://celestrak.org/NORAD/documentation/spacetrk.pdf>.
- **Brouwer, D. (1959).** Solution of the problem of artificial satellite theory without drag.
  *Astron. J.* 64, 378. <https://doi.org/10.1086/107958>. J2, J2^2, J4 secular rates.
- **Izzo, D. (2015).** Revisiting Lambert's problem. *Celest. Mech. Dyn. Astr.* 121, 1-15.
  <https://doi.org/10.1007/s10569-014-9587-y>.
- **Conway, B. A. (1986).** An improved algorithm due to Laguerre for the solution of Kepler's
  equation. *Celest. Mech.* 39, 199-211. <https://doi.org/10.1007/BF01230852>.
- **Danby, J. M. A. (1987).** The solution of Kepler's equation, III. *Celest. Mech.* 40, 303-312.
  <https://doi.org/10.1007/BF01235847>; and *Fundamentals of Celestial Mechanics*, 2nd ed.,
  Willmann-Bell, 1988. Starting value for Kepler's equation.
- **Broucke, R. A. & Cefola, P. J. (1972).** On the equinoctial orbit elements. *Celest. Mech.* 5,
  303-310. <https://doi.org/10.1007/BF01228432>.
- **Fehlberg, E. (1968).** *Classical fifth-, sixth-, seventh-, and eighth-order Runge-Kutta
  formulas with stepsize control*. NASA TR R-287. <https://ntrs.nasa.gov/citations/19680027281>. RKF7(8).
- **Hairer, E., Norsett, S. P. & Wanner, G. (1993).** *Solving Ordinary Differential Equations I:
  Nonstiff Problems*, 2nd ed. Springer. <https://doi.org/10.1007/978-3-540-78862-1>. Dense output
  and Hermite interpolation.
- **Shampine, L. F. & Thompson, S. (2000).** Event location for ordinary differential equations.
  *Computers & Mathematics with Applications* 39(5-6), 43-54.
  <https://doi.org/10.1016/S0898-1221%2800%2900045-6>. Locating events (impact) on the dense output.
- **Brent, R. P. (1973).** *Algorithms for Minimization without Derivatives*. Prentice-Hall
  (Dover reprint, 2002). Brent's method for event times (impacts, passes, rise and set) and times
  of closest approach.
- **Shepperd, S. W. (1978).** Quaternion from rotation matrix. *J. Guidance and Control* 1(3),
  223-224. <https://doi.org/10.2514/3.55767b>.
- **Clohessy, W. H. & Wiltshire, R. S. (1960).** Terminal guidance system for satellite
  rendezvous. *J. Aerospace Sci.* 27(9), 653-658. <https://doi.org/10.2514/8.8704>. Relative motion in RIC.

## Conjunctions and collision probability

- **Foster, J. L. & Estes, H. S. (1992).** *A Parametric Analysis of Orbital Debris Collision
  Probability and Maneuver Rate for Space Vehicles*. NASA/JSC-25898, NASA Johnson Space Center
  (not on NTRS; no link). The short-term encounter model: the 2-D Gaussian over the hard-body
  circle in the encounter plane.
- **Hoots, F. R., Crawford, L. L. & Roehrich, R. L. (1984).** An analytic method to determine
  future close approaches between satellites. *Celestial Mechanics* 33(2), 143-158.
  <https://doi.org/10.1007/BF01234152>. Geometric pre-filters for screening a catalogue.
- **Alfano, S. & Negron, D. (1993).** Determining satellite close approaches. *J. Astronaut. Sci.*
  41(2), 217-225 (no DOI found). Close approaches as zeros of the range rate, found from cubic
  polynomials through its samples; here the samples bracket each zero for Brent's method.
- **Akella, M. R. & Alfriend, K. T. (2000).** Probability of collision between space objects.
  *J. Guidance, Control, and Dynamics* 23(5), 769-772. <https://doi.org/10.2514/2.4611>.
- **Patera, R. P. (2003).** Satellite collision probability for nonlinear relative motion.
  *J. Guidance, Control, and Dynamics* 26(5), 728-733. <https://doi.org/10.2514/2.5127>. For slow
  encounters whose relative motion is not a straight line, which the short-term model excludes.
- **Alfano, S. (2005a).** A numerical implementation of spherical object collision probability.
  *J. Astronaut. Sci.* 53(1), 103-109. <https://doi.org/10.1007/BF03546397>. The 1-D integral
  with an exact error-function inner integral.
- **Alfano, S. (2005b).** Relating position uncertainty to maximum conjunction probability.
  *J. Astronaut. Sci.* 53(2), 193-205. <https://doi.org/10.1007/BF03546350>. Probability dilution:
  the covariance size that maximizes Pc.
- **Chan, F. K. (2008).** *Spacecraft Collision Probability*. The Aerospace Press / AIAA.
  <https://doi.org/10.2514/4.989186>. Background, including probability dilution.
- **Alfano, S. (2009).** Satellite conjunction Monte Carlo analysis. AAS/AIAA Space Flight
  Mechanics Meeting, Savannah, GA (AAS 09-233; no DOI found). Test cases with Monte Carlo
  probabilities, used in the tests through CARA's copies.
- **Balch, M. S., Martin, R. & Ferson, S. (2019).** Satellite conjunction analysis and the false
  confidence theorem. *Proc. R. Soc. A* 475(2227), 20180565. <https://doi.org/10.1098/rspa.2018.0565>.
  Why a small Pc from a large covariance does not mean the approach is safe.
- **NASA CARA Analysis Tools.** <https://github.com/nasa/CARA_Analysis_Tools>. Reference
  inputs and probabilities for the tests (`PcCircle_UnitTest.m`, Alfano 2009 cases).
- **Rice, S. O. (1945).** Mathematical analysis of random noise (part III). *Bell System Technical
  Journal* 24(1), 46-156. <https://doi.org/10.1002/j.1538-7305.1945.tb00453.x>. The distribution
  of the distance of an offset isotropic Gaussian point, used as a reference in the tests.
- **Piessens, R., de Doncker-Kapenga, E., Ueberhuber, C. W. & Kahaner, D. K. (1983).**
  *QUADPACK: A Subroutine Package for Automatic Integration*. Springer. The Gauss-Kronrod 7-15 rule.

## Estimation and covariance

- **Bar-Shalom, Y., Li, X. R. & Kirubarajan, T. (2001).** *Estimation with Applications to Tracking
  and Navigation*. Wiley. <https://doi.org/10.1002/0471221279>. Chi-square consistency tests.
- **Junkins, J. L., Akella, M. R. & Alfriend, K. T. (1996).** Non-Gaussian error propagation in
  orbital mechanics. *J. Astronaut. Sci.* 44(4), 541-563. Cartesian vs element-space uncertainty.
- **Aristoff, J. M., Horwood, J. T. & Alfriend, K. T. (2021).** On a set of J2 equinoctial orbital
  elements and their use for uncertainty propagation. *Celest. Mech. Dyn. Astr.* 133, 9.
  <https://doi.org/10.1007/s10569-021-10004-0>. Why the element-space consistency test holds up.
- **Sabol, C. et al. (2010).** Linearized orbit covariance generation and propagation analysis via
  simple Monte Carlo simulations. 20th AAS/AIAA Space Flight Mechanics Meeting.
- **Julier, S. J. & Uhlmann, J. K. (1997).** A non-divergent estimation algorithm in the presence of
  unknown correlations. *Proc. American Control Conf.*, 2369-2373.
  <https://doi.org/10.1109/ACC.1997.609105>. Covariance intersection.
