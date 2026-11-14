#!/usr/bin/env python3
###########################################################################
# check_astro_header.py
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
"""Check that astro/astro.h includes every public header.

Every header under include/astro/ (and the generated astro/version.h and astro/export.h) must be
included by include/astro/astro.h, and astro.h must include nothing that does not exist.

Usage: check_astro_header.py [source directory]   (default: the parent of tools/)
Run as the ctest test `headers.astro_h`.
"""

import pathlib
import re
import sys

GENERATED = {"astro/version.h", "astro/export.h"}  # made by CMake in the build tree


def main() -> int:
    root = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else pathlib.Path(__file__).resolve().parent.parent
    include = root / "include"
    umbrella = include / "astro" / "astro.h"

    headers = {p.relative_to(include).as_posix() for p in (include / "astro").rglob("*.h")}
    headers.discard("astro/astro.h")
    headers |= GENERATED
    included = set(re.findall(r'^#include "([^"]+)"', umbrella.read_text(), flags=re.MULTILINE))

    missing = sorted(headers - included)
    unknown = sorted(included - headers)
    for h in missing:
        print(f"check_astro_header: {h} is not included by astro/astro.h", file=sys.stderr)
    for h in unknown:
        print(f"check_astro_header: astro/astro.h includes {h}, which does not exist", file=sys.stderr)
    print(f"check_astro_header: {len(included)} headers included, {len(missing)} missing, {len(unknown)} unknown")
    return 1 if missing or unknown else 0


if __name__ == "__main__":
    sys.exit(main())
