/**********************************************************************
 * astro.h
 *
 * Copyright 2024 W. Greg Visotski
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 *         limitations under the License.
 **********************************************************************/

#pragma once

/// \file
/// Every public header of the library, for convenience: `#include "astro/astro.h"` gives the whole
/// interface in one line. It costs compile time, because every header is read for every file that
/// includes it. Including the specific headers a file needs is still the norm, and the rule inside
/// the library. Grouped by module, in the order of the module list in the documentation; a test
/// checks that no header under include/astro/ is missing.

#include "astro/export.h"
#include "astro/version.h"

// core
#include "astro/core/angle.h"
#include "astro/core/constants.h"
#include "astro/core/error.h"
#include "astro/core/mat6.h"
#include "astro/core/rotation.h"
#include "astro/core/vec.h"

// time
#include "astro/time/julian_date.h"
#include "astro/time/leap_seconds.h"
#include "astro/time/time.h"

// geodesy
#include "astro/geodesy/ellipsoid.h"
#include "astro/geodesy/geodesic.h"
#include "astro/geodesy/geodetic.h"
#include "astro/geodesy/helmert.h"
#include "astro/geodesy/local.h"

// frames
#include "astro/frames/earth_rotation.h"
#include "astro/frames/eop.h"
#include "astro/frames/precession_nutation.h"
#include "astro/frames/rotation_series.h"
#include "astro/frames/transform.h"

// ephemeris
#include "astro/ephemeris/ephemeris.h"
#include "astro/ephemeris/low_precision.h"
#include "astro/ephemeris/spk.h"

// astrometry
#include "astro/astrometry/apparent.h"
#include "astro/astrometry/refraction.h"
#include "astro/astrometry/rise_set.h"
#include "astro/astrometry/spherical.h"

// orbit
#include "astro/orbit/central_body.h"
#include "astro/orbit/close_approach.h"
#include "astro/orbit/elements.h"
#include "astro/orbit/gravity.h"
#include "astro/orbit/impact.h"
#include "astro/orbit/integration.h"
#include "astro/orbit/kepler.h"
#include "astro/orbit/numerical.h"
#include "astro/orbit/secular.h"
#include "astro/orbit/sgp4.h"
#include "astro/orbit/state12.h"
#include "astro/orbit/state_covariance.h"
#include "astro/orbit/tle.h"
#include "astro/orbit/trajectory.h"
#include "astro/orbit/two_body.h"

// satellite
#include "astro/satellite/eclipse.h"
#include "astro/satellite/estimate_tracker.h"
#include "astro/satellite/ground.h"
#include "astro/satellite/look.h"
#include "astro/satellite/orbital_frames.h"
