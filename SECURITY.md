# Security policy

## Supported versions

astro has not had a release yet. Fixes go to the `main` branch.

## Reporting a vulnerability

Please report security problems privately, not in a public issue:

1. Open the repository's **Security** tab.
2. Choose **Report a vulnerability**.

This uses GitHub's private vulnerability reporting. Only the maintainer sees the report.

Useful details:

| Item        | Example                                                     |
|-------------|-------------------------------------------------------------|
| Affected    | function, file or loader (e.g. `EopTable::loadFinals2000A`) |
| Input       | the file or values that trigger it                          |
| Effect      | crash, out-of-bounds read, hang, wrong result               |
| Environment | compiler, version, OS                                       |

You should get a reply within 7 days.

## Scope

The library does no networking. The main attack surface is its file and text parsers:

- SPK ephemeris files (`.bsp`)
- EOP files (`finals2000A`, C04)
- leap-second tables
- TLEs

A crash or memory error on malformed input is in scope.
