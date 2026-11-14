# Data files

The library never touches the network. What it needs from outside comes from files that your
application downloads, keeps up to date, and loads. This page lists each file:
- what it is for
- where to get it
- how often it changes
- what you lose without it
- how to tell when it runs out

## What is built in

| Built in | Covers | Good until |
|---|---|---|
| Leap-second table (`LeapSecondTable::builtin()`) | TAI - UTC from 1960, including the "rubber second" years before 1972, up to the 2017-01-01 leap second (37 s) | the next leap second the IERS announces |
| IAU 2006/2000A precession-nutation series | the celestial pole, for every frame conversion | always (a model, not data) |
| TDB - TT series (Fairhead and Bretagnon) | TDB and TCB | always |
| Low-precision Sun and Moon (`LowPrecisionEphemeris`) | Sun to 30", Moon to 6' (1950-2050, against DE440; see [Accuracy](../accuracy.md)) | always |
| Constants (WGS 84, IERS 2010, GM values) | geodesy, gravity | always |

Everything else is a file you supply.

## What you supply

| File | Needed for | Source | Size | Changes | Without it |
|---|---|---|---|---|---|
| `finals2000A.all` | UT1, polar motion and pole offsets: ECEF <-> ECI to the centimetre | [IERS Rapid Service](https://datacenter.iers.org/products/eop/rapid/standard/finals2000A.all) (mirror at USNO) | 3.8 MB | daily | `ConstantEop{}`: UT1 = UTC, up to 0.9 s off, which is up to about 420 m east-west on the ground. No polar motion either: up to about 15 m. |
| `eopc04.1962-now` | the same, final values for reprocessing past data | [IERS EOP C04](https://hpiers.obspm.fr/iers/eop/eopc04/eopc04.1962-now) | 5 MB | twice a week; lags real time by weeks, no predictions | use `finals2000A.all` |
| `Leap_Second.dat` or `leap-seconds.list` | UTC after a leap second announced since 2017 | [IERS](https://hpiers.obspm.fr/iers/bul/bulc/Leap_Second.dat), [IANA](https://data.iana.org/time-zones/tzdb/leap-seconds.list) | 2-5 kB | when a leap second is announced ([Bulletin C](https://datacenter.iers.org/data/latestVersion/bulletinC.txt), twice a year) | the built-in table, which is correct until the next leap second |
| `de440s.bsp` | Sun, Moon and planets to the accuracy of JPL's DE440 | [NAIF](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp) | 32 MB | never (covers 1849-2150) | the low-precision Sun and Moon; no planets |
| Two-line element sets | SGP4 | [CelesTrak](https://celestrak.org), [Space-Track](https://www.space-track.org) (account needed) | small | several times a day per object | nothing: each TLE describes one object |

A few notes on these files:

- **Larger JPL kernels.** `de440.bsp` (1550-2650) and `de441.bsp` (13 200 BC to AD 17 191) work in
  the same way. They are 114 MB and 3 GB.
- **Mars in DE kernels.** These kernels hold the Mars system barycentre, not Mars itself:
  ask for `Body::MarsBarycenter`.
- **TLE age.** A TLE is a fit to recent observations, and its error grows with every day after its
  epoch. Fetch fresh ones before an analysis, and look at the epoch before trusting one.

## Getting them

`scripts/fetch_iers_data.sh` downloads `finals2000A.all`, `Leap_Second.dat` and `leap-seconds.list`.
- `--c04` adds the C04 series and `--de440s` adds the JPL kernel; `--all` fetches both.
- Each file goes to a temporary name first. It is checked, and only then replaces the previous copy,
  so a failed or truncated download never destroys a good file.
- The exit status is not zero if any file could not be updated, so cron reports it.

```sh
scripts/fetch_iers_data.sh -o /var/lib/myapp/astro          # once, by hand
# crontab: every day at 03:17, quietly, skipping files fetched in the last day
17 3 * * * /opt/myapp/fetch_iers_data.sh -q -m 1 -o /var/lib/myapp/astro
```

The library has no default search path. Your application decides where the files live, for example
a directory named in its configuration or in an environment variable.

## Loading them

Load the files once, at start-up, and keep the objects:

<!-- snippet: examples/data_files.cpp#startup -->
```cpp
// Once at start-up, before any UTC conversion: a newer leap-second table, if there is one.
if (fs::exists(dir / "Leap_Second.dat")) {
    setLeapSecondTable(LeapSecondTable::load((dir / "Leap_Second.dat").string()));
}

// Earth orientation: the IERS file, or UT1 = UTC with no polar motion (about 400 m at the equator).
std::shared_ptr<const EopProvider> eop = std::make_shared<ConstantEop>();
if (fs::exists(dir / "finals2000A.all")) {
    eop = std::make_shared<EopTable>(EopTable::loadFinals2000A((dir / "finals2000A.all").string()));
}

// Sun, Moon and planets: a JPL kernel, or the built-in low-precision Sun and Moon.
std::shared_ptr<const Ephemeris> eph = std::make_shared<LowPrecisionEphemeris>();
if (fs::exists(dir / "de440s.bsp")) {
    eph = std::make_shared<SpkEphemeris>((dir / "de440s.bsp").string());
}
```

- **Leap seconds:** `setLeapSecondTable` replaces the table every UTC conversion uses. The swap is
  atomic, so other threads may keep converting while it happens. Do it before the first conversion,
  so that every result uses the same table.
- **Earth orientation:** an `EopTable` is immutable and safe to share between threads. Pass it to
  `FrameTransformer`, `itrsToGcrs`, `Time::to(TimeScale::UT1, eop)` and the other functions that take
  an `EopProvider`.
- **Ephemeris:** an `SpkEphemeris` keeps the file open and reads records on demand. It is
  thread-safe, so one instance can serve the whole program.

## Knowing when data runs out

Each kind of data behaves differently when it runs out:

| Data | When it runs out | What the library does |
|---|---|---|
| `EopTable` | after its last record (the predictions in `finals2000A.all` reach about a year ahead) | throws `astro::Error`. It never extrapolates, so a stale file is an error, not a quietly wrong answer. |
| Leap-second table | after a leap second it does not know about | keeps the last offset. UTC times after the leap second are 1 s off, with no error. |
| Leap-second file | after the expiry date in the file (`expiresMjd()`) | nothing. It is your signal to fetch a newer file. |
| `SpkEphemeris` | outside the kernel's coverage | throws `astro::Error` |
| TLE | no hard limit | nothing. Judge by the age of its epoch. |

The leap-second case is the one that fails silently, so check it at start-up, with the EOP
coverage:

<!-- snippet: examples/data_files.cpp#validity -->
```cpp
// How long is the loaded data good for? Check at start-up, not when a conversion fails.
const double today = Time::now().mjd(); // MJD (UTC)
if (const auto expires = leapSecondTable()->expiresMjd()) {
    std::cout << "leap-second table expires in " << std::fixed << std::setprecision(0) << (*expires - today)
              << " days\n";
} else {
    std::cout << "built-in leap-second table: no expiry date; it ends with the 2017 leap second\n";
}
if (const auto* table = dynamic_cast<const EopTable*>(eop.get())) {
    std::cout << "EOP data (with predictions) run out in " << std::fixed << std::setprecision(0)
              << (table->lastMjd() - today) << " days\n";
} else {
    std::cout << "no EOP file: UT1 = UTC, no polar motion\n";
}
```

## Reproducible results

EOP values are revised as better observations arrive. A prediction becomes a rapid value, and a
rapid value later becomes a final one, so a conversion run today and again next month can differ
at the centimetre level. If you need to reproduce a result exactly, keep a copy of the files you
used (or their dates and checksums) with the result.

## Formats accepted

| Loader | Formats |
|---|---|
| `EopTable::loadFinals2000A` / `parseFinals2000A` | IERS `finals2000A.all`, `.data`, `.daily` (Bulletin A columns) |
| `EopTable::loadC04` / `parseC04` | IERS C04, the "14" and the "20" layouts |
| `LeapSecondTable::load` / `parse` | IERS `Leap_Second.dat`, IETF/NIST `leap-seconds.list`, USNO `tai-utc.dat` |
| `SpkEphemeris` | JPL/NAIF SPK (`.bsp`), segment types 2 and 3, J2000/ICRF frame, either byte order |
| `Tle::parse` / `Tle::parseMany` | two-line and three-line element sets, Alpha-5 catalogue numbers |

All loaders are strict:
- A malformed line is an error that gives the line number. It is never skipped, because a skipped
  line could leave a wrong value in place without notice.
- Every loader also has a form that reads from memory (`parse...`), for files your application
  fetched itself.
