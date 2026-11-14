# Contributing to astro

Thank you for helping. This page says how to build and test the library, and what a change needs
before it is merged. Report a security problem privately, as [SECURITY.md](SECURITY.md) describes,
not in a public issue.

## Build and test

```sh
cmake -B build -DCMAKE_BUILD_TYPE=Release -DASTRO_WARNINGS_AS_ERRORS=ON
cmake --build build -j
ctest --test-dir build --output-on-failure
```

The CMake options are listed in [Getting started](docs/guide/getting-started.md). A few you may want:

| Option | What it builds |
|---|---|
| `ASTRO_BUILD_STUDIES=ON` | the accuracy studies behind [Accuracy](docs/accuracy.md) |
| `ASTRO_BUILD_BENCHMARKS=ON` | the timing programs behind [Performance](docs/performance.md) |
| `ASTRO_BUILD_FUZZERS=ON` | fuzz targets for the parsers (Clang; see `fuzz/CMakeLists.txt`) |
| `ASTRO_COVERAGE=ON` | a line-coverage build (GCC; report with `tools/coverage_report.py`) |

## What a change needs

### Correctness

- **A fix comes with a test that fails without it.** Write the test first, see it fail, then fix the
  cause (not the example that showed it).
- **A test compares with an independent truth**: a closed form, an identity the code does not use, a
  published value (say where it is published), or a calculation in higher precision written in the
  test. Never with the output of the code under test.
- **Each tolerance has a reason**, in a comment: what the code achieves, or the precision of the
  reference.
- **A test must be able to fail.** When you add one, break the code on purpose in the way the test is
  meant to catch, and check that it does.

### Interface

- Follow the conventions in [docs/guide/conventions.md](docs/guide/conventions.md): SI units and
  radians, the naming rules, the argument order of the propagators, and the error style.
- Every public declaration has a `///` comment that says what each parameter is (with its unit) and
  what throws. The documentation build treats Doxygen warnings as errors.
- A new formula cites its source in a comment, and the source goes in
  [docs/references.md](docs/references.md). Code that follows another project's code is credited in
  [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md), with its licence.

### Documentation

- If a number in `docs/` changes, rerun the study or benchmark that produced it and update it.
- The code in the guides must match the examples, which are built and run as tests:
  `python3 tools/check_doc_snippets.py --fix` copies the example code into the guides.
- Add a line to [CHANGELOG.md](CHANGELOG.md) under the next release, in the section that fits.
  Breaking changes say what users must change.

## Style

- **Formatting:** clang-format 18 with the repository's `.clang-format`:
  `clang-format -i $(git ls-files '*.h' '*.cpp')`.
- **Alignment:** the code is kept column-aligned on purpose, because it is easier to read. Put a
  hand-aligned table of constants between `// clang-format off` and `// clang-format on`.
- **Names:** `camelCase` functions and variables, `PascalCase` types, `kName` for constants local to
  a source file.
- **Comments:** plain words; spell out an acronym the first time it is used.

## Commits and pull requests

- **One change per pull request**, and one idea per commit.
- **Commit messages:** an imperative subject line under about 72 characters; then a body that says
  what was wrong or missing, why, what changed, and how it was tested.
- **Signing:** commits are signed.
- **CI must pass.** It checks:
  - GCC and Clang builds, with warnings as errors;
  - the tests under the thread, address and undefined-behaviour sanitizers;
  - a minute of fuzzing for each parser;
  - line coverage (it fails below the floor in the workflow);
  - the `.deb` and `.rpm` packages and a program built against them;
  - formatting;
  - the documentation build.

## Releases

Each release is tagged `vMAJOR.MINOR.PATCH` on `main`, with the version set in `CMakeLists.txt` and
the changes listed in `CHANGELOG.md` under that version and its date.
