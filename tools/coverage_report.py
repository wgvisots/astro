#!/usr/bin/env python3
###########################################################################
# coverage_report.py
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
"""Line coverage of the library by the test suite, from a build made with -DASTRO_COVERAGE=ON.

Usage:
    cmake -B build-cov -DCMAKE_BUILD_TYPE=Debug -DASTRO_COVERAGE=ON
    cmake --build build-cov -j && ctest --test-dir build-cov
    python3 tools/coverage_report.py build-cov [--min PERCENT] [--uncovered]

Prints a Markdown table: lines that can run, and lines that ran, for each module (its sources in
src/<module> and its headers in include/astro/<module>). With --min the exit status is 1 if the
total is below PERCENT. With --uncovered it also lists, for each file, the lines that never ran.

It needs only gcov (part of GCC) and the Python standard library.
"""

import argparse
import collections
import json
import pathlib
import subprocess
import sys


def gcov_files(gcda: pathlib.Path):
    """The per-file line counts gcov reports for one .gcda file."""
    out = subprocess.run(["gcov", "--json-format", "--stdout", gcda.name], cwd=gcda.parent, capture_output=True, text=True)
    if out.returncode != 0:
        sys.exit(f"gcov failed on {gcda}: {out.stderr.strip()}")
    for document in out.stdout.splitlines():
        if document.strip():
            data = json.loads(document)
            cwd = pathlib.Path(data.get("current_working_directory", gcda.parent))
            for entry in data["files"]:
                yield (cwd / entry["file"]).resolve(), entry["lines"]


def main(argv):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("build", type=pathlib.Path, help="build directory configured with -DASTRO_COVERAGE=ON")
    parser.add_argument("--min", type=float, default=None, help="fail if the total line coverage is below this percentage")
    parser.add_argument("--uncovered", action="store_true", help="list the lines that never ran")
    args = parser.parse_args(argv)

    root = pathlib.Path(__file__).resolve().parent.parent
    library = [root / "src", root / "include" / "astro"]
    gcdas = sorted(args.build.resolve().rglob("*.gcda"))
    if not gcdas:
        sys.exit(f"no .gcda files under {args.build}: configure with -DASTRO_COVERAGE=ON, build, and run the tests")

    # Executions of each line, summed over every object file that contains it (headers are compiled
    # into many).
    counts = collections.defaultdict(lambda: collections.defaultdict(int))
    for gcda in gcdas:
        for path, lines in gcov_files(gcda):
            if any(base in path.parents for base in library):
                for line in lines:
                    counts[path][line["line_number"]] += line["count"]

    modules = collections.defaultdict(lambda: [0, 0])
    for path, lines in counts.items():
        relative = path.relative_to(root)
        module = relative.parts[1] if relative.parts[0] == "src" else relative.parts[2]
        modules[module][0] += len(lines)
        modules[module][1] += sum(1 for count in lines.values() if count > 0)

    total = [sum(m[0] for m in modules.values()), sum(m[1] for m in modules.values())]
    print("| Module | Lines | Covered | Coverage |")
    print("|---|---|---|---|")
    for module in sorted(modules):
        lines, covered = modules[module]
        print(f"| {module} | {lines} | {covered} | {100.0 * covered / lines:.1f} % |")
    percent = 100.0 * total[1] / total[0]
    print(f"| **all** | {total[0]} | {total[1]} | **{percent:.1f} %** |")

    if args.uncovered:
        print()
        for path in sorted(counts):
            missed = sorted(number for number, count in counts[path].items() if count == 0)
            if missed:
                print(f"{path.relative_to(root)}: {', '.join(map(str, missed))}")

    if args.min is not None and percent < args.min:
        print(f"\nline coverage {percent:.1f} % is below the required {args.min:.1f} %", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
