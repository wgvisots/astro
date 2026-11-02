#!/usr/bin/env bash
###########################################################################
# fetch_iers_data.sh
#
# Copyright 2024 W. Greg Visotski
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
###########################################################################
#
# Download the external data files used by the astro library. The library
# itself never touches the network; run this script (e.g. from cron) and point
# the loaders at the output directory:
#
#   finals2000A.all     -> astro::EopTable::loadFinals2000A()
#   eopc04.1962-now     -> astro::EopTable::loadC04()             (--c04)
#   Leap_Second.dat     -> astro::LeapSecondTable::load()
#   leap-seconds.list   -> astro::LeapSecondTable::load()
#   de440s.bsp          -> astro::SpkEphemeris                    (--de440s, ~32 MB)
#
# Each file is downloaded to a temporary name, sanity-checked, and only then
# moved into place, so a failed download never clobbers a good file.

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: fetch_iers_data.sh [options]

Download IERS Earth orientation and leap-second data (and optionally the
JPL DE440s ephemeris) for the astro library.

Options:
  -o, --output DIR   Output directory (default: ./data)
      --c04          Also fetch the IERS EOP 20 C04 series (eopc04.1962-now, ~5 MB)
      --de440s       Also fetch the JPL DE440s planetary ephemeris (de440s.bsp, ~32 MB)
      --all          Fetch everything (same as --c04 --de440s)
  -m, --max-age DAYS Skip files younger than DAYS days (default: 0 = always fetch)
  -q, --quiet        Only print errors
  -h, --help         Show this help

Files always fetched: finals2000A.all, Leap_Second.dat, leap-seconds.list.
Exit status is non-zero if any requested file could not be updated.
EOF
}

OUT="./data"
WANT_C04=0
WANT_DE=0
MAX_AGE=0
QUIET=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        -o|--output) OUT="${2:?missing directory}"; shift 2 ;;
        --c04) WANT_C04=1; shift ;;
        --de440s) WANT_DE=1; shift ;;
        --all) WANT_C04=1; WANT_DE=1; shift ;;
        -m|--max-age) MAX_AGE="${2:?missing days}"; shift 2 ;;
        -q|--quiet) QUIET=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "fetch_iers_data.sh: unknown option '$1'" >&2; usage >&2; exit 2 ;;
    esac
done

if ! [[ "$MAX_AGE" =~ ^[0-9]+$ ]]; then
    echo "fetch_iers_data.sh: --max-age must be a non-negative integer" >&2
    exit 2
fi

log() { [[ $QUIET -eq 1 ]] || echo "$@"; }
err() { echo "fetch_iers_data.sh: $*" >&2; }

if command -v curl >/dev/null 2>&1; then
    download() { curl --fail --silent --show-error --location --retry 3 --retry-delay 5 \
                      --connect-timeout 30 --max-time 900 -o "$2" "$1"; }
elif command -v wget >/dev/null 2>&1; then
    download() { wget --quiet --tries=3 --timeout=60 -O "$2" "$1"; }
else
    err "need curl or wget"
    exit 1
fi

mkdir -p "$OUT"

# Sanity checks: each prints nothing and returns 0 if the file looks right.
check_finals()  { head -n 1 "$1" | grep -Eq '^[ 0-9]{6} +[0-9]{5}\.[0-9]{2} '; }
check_leapdat() { grep -q 'File expires on' "$1"; }
check_leaplst() { grep -q '^#@' "$1"; }
check_c04()     { grep -q 'C04' "$1" && grep -Eq '^ *[0-9]{4} +[0-9]+ +[0-9]+ +[0-9]+ +[0-9]{5}\.' "$1"; }
check_bsp()     { [[ "$(head -c 8 "$1")" == "DAF/SPK " ]]; }

# fetch NAME CHECK URL [URL...]   (tries each URL in turn)
FAILED=0
fetch() {
    local name="$1" check="$2"; shift 2
    local dest="$OUT/$name"
    if [[ $MAX_AGE -gt 0 && -f "$dest" ]] && [[ -n "$(find "$dest" -mtime "-$MAX_AGE" 2>/dev/null)" ]]; then
        log "  $name: up to date (younger than $MAX_AGE days), skipped"
        return 0
    fi
    local tmp
    tmp="$(mktemp "$OUT/.$name.XXXXXX")"
    local url
    for url in "$@"; do
        if download "$url" "$tmp" 2>/dev/null && [[ -s "$tmp" ]] && "$check" "$tmp"; then
            chmod 644 "$tmp"
            mv -f "$tmp" "$dest"
            log "  $name: $(wc -c < "$dest" | tr -d ' ') bytes from $url"
            return 0
        fi
        log "  $name: failed from $url"
    done
    rm -f "$tmp"
    err "could not update $name (kept any existing copy)"
    FAILED=1
    return 0
}

log "Downloading to $OUT"
fetch finals2000A.all check_finals \
    "https://datacenter.iers.org/products/eop/rapid/standard/finals2000A.all" \
    "https://maia.usno.navy.mil/ser7/finals2000A.all"
fetch Leap_Second.dat check_leapdat \
    "https://hpiers.obspm.fr/iers/bul/bulc/Leap_Second.dat" \
    "https://hpiers.obspm.fr/eoppc/bul/bulc/Leap_Second.dat"
fetch leap-seconds.list check_leaplst \
    "https://data.iana.org/time-zones/tzdb/leap-seconds.list" \
    "https://hpiers.obspm.fr/iers/bul/bulc/ntp/leap-seconds.list"
if [[ $WANT_C04 -eq 1 ]]; then
    fetch eopc04.1962-now check_c04 \
        "https://hpiers.obspm.fr/iers/eop/eopc04/eopc04.1962-now"
fi
if [[ $WANT_DE -eq 1 ]]; then
    fetch de440s.bsp check_bsp \
        "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp" \
        "https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/de440s.bsp"
fi

# Summaries useful for monitoring.
if [[ $QUIET -eq 0 ]]; then
    if [[ -f "$OUT/Leap_Second.dat" ]]; then
        echo "Leap_Second.dat: $(grep -o 'File expires on.*' "$OUT/Leap_Second.dat" | head -n 1), " \
             "TAI-UTC = $(awk '!/^#/ && NF >= 5 {v=$5} END {print v}' "$OUT/Leap_Second.dat") s"
    fi
    if [[ -f "$OUT/finals2000A.all" ]]; then
        echo "finals2000A.all: MJD $(head -n 1 "$OUT/finals2000A.all" | cut -c8-15 | tr -d ' ') to" \
             "$(awk 'length($0) >= 68 && substr($0,59,10) ~ /[0-9]/ {m=substr($0,8,8)} END {gsub(/ /,"",m); print m}' \
                "$OUT/finals2000A.all") (incl. predictions)"
    fi
fi

exit $FAILED
