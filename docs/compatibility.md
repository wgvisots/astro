# Compatibility

This page states what a version number of astro promises, what it does not, how parts of the API
are retired, and which platforms are supported. It is the policy every release follows.

## Version numbers

astro uses [Semantic Versioning 2.0.0](https://semver.org): `MAJOR.MINOR.PATCH`.

| Release | 0.x (now) | 1.x and later |
|---|---|---|
| Patch (`0.3.0` to `0.3.1`) | bug fixes only; source and binary compatible | bug fixes only; source and binary compatible |
| Minor (`0.3` to `0.4`) | **may break code**: the interface is still being settled | new features; **source compatible** (see below) |
| Major (`1.x` to `2.0`) | (does not apply) | may break code; removes what was deprecated |

`astro/version.h` gives the version to the preprocessor (`ASTRO_VERSION_MAJOR`, `_MINOR` and
`_PATCH`) and at run time (`astro::versionString()`).

## 0.x: the interface is in flux

Until 1.0, any minor release may rename, move or remove a part of the API, or change how it
behaves. Such changes are never silent:
- The [CHANGELOG](https://github.com/wgvisots/astro/blob/main/CHANGELOG.md) lists every breaking
  change under its release.
- Each listing shows the old and the new code side by side.
- Patch releases never break code.

These releases are for settling the interface before it is promised.

## 1.x: source compatibility

From 1.0, **code that compiles against a 1.x release, and uses only the public API as described
below, compiles against every later 1.y release and behaves the same**. The exceptions are:

- **Bug fixes.** Behaviour that contradicts the documentation is corrected to match it.
- **Numerical results within the documented accuracy.** An improved algorithm may change results
  in the last digits, but never beyond the tolerances in [Accuracy and validation](accuracy.md).
  Any change larger than rounding is listed in the CHANGELOG.
- **Input the documentation calls invalid** may start to throw `astro::Error` where it used to
  give a meaningless result.
- **Error messages.** The text after the function name may change (see `astro/core/error.h`).

### What the public API is

The public API is everything declared in the headers under `include/astro/`, except:
- names in a `detail` namespace
- anything a comment marks as internal

Headers under `src/` are private.

### What not to rely on

A library can only grow compatibly if its users leave it room. Code that does any of the following
may break in a minor release, and that will not count as a break of the promise:

| Do not | Because a minor release may |
|---|---|
| take the address of a library function, or form a pointer to a member function | add an overload or a defaulted parameter |
| forward-declare astro's types yourself | change a class into an alias, or add template parameters with defaults |
| switch over an astro enum with no `default:` and treat a missing case as an error | add an enumerator |
| rely on one astro header including another | tidy the includes (include what you use) |
| add declarations to namespace `astro` or specialise its templates (other than where documented) | add a name that clashes |
| parse the text of `what()` | reword a message |
| depend on `sizeof`, the layout or the private members of a class | change them (binary compatibility is not promised) |

The `"Enum::Enumerator"` text written by each enum's `operator<<` is stable. Tools may read it.

## Binary compatibility

astro promises **source** compatibility, not binary compatibility:

| Release | Binary (ABI) compatibility | Shared library name |
|---|---|---|
| Patch | kept: replace the library without rebuilding | unchanged, `libastro.so.MAJOR.MINOR` |
| Minor or major | not promised: rebuild your program | new, for example `libastro.so.1.3` after `libastro.so.1.2` |

The shared library's soname carries the major and minor number. That way, a program built against
one minor release never loads another by accident.

## Deprecation and removal

| | 0.x | 1.x |
|---|---|---|
| Removing or renaming | in any minor release, listed in the CHANGELOG with the replacement | never within 1.x |
| Retiring something | (removed directly) | marked `[[deprecated("use X instead")]]` in a minor release and listed under **Deprecated** in the CHANGELOG. It keeps working for the rest of 1.x and is removed in 2.0. |

Deprecated names give a compiler warning that names the replacement. Something is deprecated for
at least one minor release before the major release that removes it.

## Supported platforms

| | Linux | Windows |
|---|---|---|
| Architecture | x86-64 (Ubuntu 24.04 in CI) | x86-64 (Windows Server 2025 in CI) |
| Compilers | GCC 13 and Clang 18, in C++20 mode | MSVC from Visual Studio 2022 or newer, in C++20 mode |
| Library types | static and shared | static and DLL |
| Packages | `.deb` and `.rpm` (CPack), Conan recipe | Conan recipe |
| Build | CMake 3.20 or newer | CMake 3.20 or newer |
| Dependency | Boost 1.77 or newer (headers only, for Boost.Describe) | the same |

Every supported combination is built and tested in CI on every change. One difference: with MSVC,
`long double` is the same as `double`, so four accuracy checks whose reference values are computed
in `long double` are skipped there (the Stumpff functions, the integrated geodesic, and two
parabolic-orbit checks). The Linux jobs run them on the same code.

Other platforms and compilers, such as macOS, MinGW on Windows, or other Linux architectures, are
not tested. The code uses only standard C++20 and Boost, so it may well work. Fixes are welcome when
they keep the code simple, but there is no promise.

Raising a minimum version (compiler, CMake or Boost) happens only in a minor release, and the
CHANGELOG says so.

## Data file formats

The file loaders keep accepting every format listed in [Data files](guide/data-files.md) for the
whole of a major version. New formats may be added in a minor release.

## Which releases get fixes

Fixes go into the latest minor release, as a patch release. Older minor releases are not
maintained: upgrading within 1.x is meant to be a rebuild, with no code changes.
