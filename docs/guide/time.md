# Time scales

`astro::Time` is an instant in a stated time scale, held as a two-part Julian date for
sub-nanosecond resolution over centuries. Converting with `to(scale)` goes through TAI.

| Scale | What it is | Relation |
|---|---|---|
| `UTC` | civil time, with leap seconds | TAI - UTC from the leap-second table (37 s since 2017) |
| `TAI` | International Atomic Time | the hub |
| `TT` | Terrestrial Time | TAI + 32.184 s |
| `TDB` | Barycentric Dynamical Time (ephemerides) | TT + periodic terms below 2 ms (Fairhead & Bretagnon series) |
| `TCG`, `TCB` | geocentric / barycentric coordinate time | linear rates L_G, L_B (IAU 2006) |
| `UT1` | Earth rotation angle time | UTC + (UT1 - UTC) from Earth orientation data |
| `GPS`, `GST` | GPS time, Galileo System Time | TAI - 19 s |
| `BDT` | BeiDou time | TAI - 33 s |

The definitions follow the IERS Conventions [IERS Conventions (2010), ch. 10](https://www.iers.org/TN36) and the
IAU resolutions summarised in [USNO Circular 179](https://aa.usno.navy.mil/publications/Circular_179); results are validated against
ERFA [ERFA](https://github.com/liberfa/erfa) to 1 ns over 1960-2035 ([Accuracy](../accuracy.md)).

## Converting

<!-- snippet: examples/time_scales.cpp#convert -->
```cpp
// A UTC instant inside the leap second at the end of 2016.
const Time utc = Time::fromIso8601("2016-12-31T23:59:60.5Z");
std::cout << utc.iso8601(3) << '\n'                     // 2016-12-31T23:59:60.500Z
          << utc.to(TimeScale::TAI).iso8601(3) << '\n'  // 2017-01-01T00:00:36.500 TAI
          << utc.to(TimeScale::TT).iso8601(3) << '\n'   // TAI + 32.184 s
          << utc.to(TimeScale::TDB).iso8601(6) << '\n'  // TT + periodic terms (< 2 ms)
          << utc.to(TimeScale::GPS).iso8601(3) << '\n'; // TAI - 19 s
```

## Leap seconds

UTC uses the SOFA "quasi-JD" convention [SOFA Time Scale and Calendar Tools](https://www.iausofa.org/s/sofa_ts_c.pdf): on a
day with a leap second, the UTC day is stretched, so `23:59:60.5` is a valid, distinct instant.
Arithmetic is in SI seconds and counts the leap second:

<!-- snippet: examples/time_scales.cpp#arithmetic -->
```cpp
// Arithmetic is in SI seconds, so it counts the leap second.
const Time before = Time::fromIso8601("2016-12-31T23:59:59Z");
const Time after  = Time::fromIso8601("2017-01-01T00:00:00Z");
std::cout << "elapsed: " << (after - before) << " s\n"; // 2 s
std::cout << (before + 1.0).iso8601(0) << '\n';         // 2016-12-31T23:59:60Z
```

The built-in table ends with the 2017 leap second. If the IERS announces a new one
([IERS Bulletin C](https://datacenter.iers.org/data/latestVersion/bulletinC.txt)), load an updated `Leap_Second.dat` or `leap-seconds.list`
[IANA leap-seconds.list](https://data.iana.org/time-zones/data/leap-seconds.list) once at start-up with `setLeapSecondTable`; the table is
swapped atomically, so other threads may keep converting. Before 1972, UTC ran at a different rate
from TAI ("rubber seconds"); those offsets are built in too.

## Unix time and GNSS weeks

Unix time counts days of exactly 86400 s, so it cannot represent a leap second;
`Time::fromUnix` / `unixSeconds` follow that convention (the instant inside a leap second maps to
the following second). GNSS week and time of week are available for GPS, Galileo and BeiDou
[IS-GPS-200](https://www.gps.gov/interface-control-documents-icds-interface-specifications-iss):

<!-- snippet: examples/time_scales.cpp#unix_gnss -->
```cpp
const Time        t   = Time::fromUnix(1700000000.0); // Unix time ignores leap seconds
const GnssWeekTow gps = t.gnssWeekTow(GnssSystem::GPS);
std::cout << t.iso8601(0) << "  GPS week " << gps.week << ", TOW " << gps.tow << " s\n";
std::cout << "back to Unix: " << std::setprecision(12)
          << Time::fromGnssWeekTow(GnssSystem::GPS, gps.week, gps.tow).unixSeconds() << '\n';
```

Broadcast messages carry the week modulo 1024 (or 8192); `resolveWeekRollover` recovers the full
week from a nearby reference time.

## UT1

UT1 follows the Earth's actual rotation and needs Earth orientation data. `ConstantEop{}`
(UT1 = UTC) is within 0.9 s; an IERS `finals2000A.all` file makes it exact to the data's accuracy
[IERS Earth orientation data](https://www.iers.org/iers/en/dataproducts/earthorientationdata/eop):

<!-- snippet: examples/time_scales.cpp#ut1 -->
```cpp
// UT1 needs Earth orientation data (UT1 - UTC); here a constant -0.1 s instead of an IERS file.
EopValues v;
v.dut1 = -0.1;
const ConstantEop eop(v);
std::cout << t.to(TimeScale::UT1, eop).iso8601(3) << '\n';
```

## Pitfalls

- A `Time` remembers its scale, and printing shows the scale it was created in.
- `t1 - t2` converts `t2` to the scale of `t1` and gives the difference in seconds of that scale
  (a UTC `t1` goes through TAI, so leap seconds are counted).
    - Among UTC, TAI, TT and the GNSS scales, which all tick at the same rate, the result is the
      elapsed time whatever the mix.
    - TDB, TCG and TCB tick at slightly different rates, so a difference that involves one of them
      depends on which operand is on the left, by up to 1.5e-8 of the interval.
    - A difference or comparison between UT1 and another scale needs Earth orientation data and
      throws: convert with an `EopProvider` first.
- UT1 without Earth orientation data is UTC, which is up to 0.9 s off: about 400 m of Earth
  rotation at the equator.
- TDB - TT is periodic (largest term 1.66 ms over a year); the topocentric correction is a
  further ~2 us and is available through `tdbMinusTt`.
