#!/usr/bin/env python3
###########################################################################
# check_exports.py
#
# Copyright 2024-2026 W. Greg Visotski
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
"""Check that the shared library exports no internal symbol.

The library is built with hidden visibility and marks its public interface with ASTRO_EXPORT. This
script lists the symbols the shared library exports (`nm -DC --defined-only`) and fails if one of
them is internal: anything in an `astro::...::detail` namespace or in an anonymous namespace.

Usage: check_exports.py path/to/libastro.so [--nm nm]
Run as the ctest test `exports.internal` when BUILD_SHARED_LIBS is on.
"""

import argparse
import re
import subprocess
import sys

# A symbol in a detail namespace of astro (astro::detail::f, astro::orbit::detail::f, and template
# instances such as std::vector<astro::detail::X>), or in an anonymous namespace.
INTERNAL = re.compile(r"\bastro::(?:\w+::)*detail::|\(anonymous namespace\)")


def exported_symbols(library: str, nm: str) -> list[str]:
    """The demangled names of the symbols `library` defines in its dynamic symbol table."""
    out = subprocess.run([nm, "-DC", "--defined-only", library], check=True, capture_output=True, text=True).stdout
    names = []
    for line in out.splitlines():
        parts = line.split(maxsplit=2)  # address, type, name
        if len(parts) == 3:
            names.append(parts[2])
    return names


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("library", help="the shared library (libastro.so)")
    parser.add_argument("--nm", default="nm", help="the nm program (default: nm)")
    args = parser.parse_args()

    names = exported_symbols(args.library, args.nm)
    public = [n for n in names if n.startswith("astro::")]
    if not public:
        # An empty list means nm failed to read the library, or nothing is exported: either way
        # the check would prove nothing.
        print(f"check_exports: {args.library} exports no astro:: symbol", file=sys.stderr)
        return 1
    leaks = [n for n in names if INTERNAL.search(n)]
    for name in leaks:
        print(f"check_exports: internal symbol exported: {name}", file=sys.stderr)
    print(f"check_exports: {len(names)} symbols exported, {len(public)} in astro::, {len(leaks)} internal")
    return 1 if leaks else 0


if __name__ == "__main__":
    sys.exit(main())
