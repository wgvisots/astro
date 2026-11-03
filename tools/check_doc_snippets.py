#!/usr/bin/env python3
###########################################################################
# check_doc_snippets.py
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
"""Keep the code in the Markdown docs identical to the compiled examples.

A code block in docs/**/*.md that follows a marker

    <!-- snippet: examples/time_scales.cpp#convert -->
    ```cpp
    ...
    ```

must equal the region of examples/time_scales.cpp between the lines `//! [convert]` (Doxygen's
snippet syntax), dedented. The examples are built and run as tests, so the documented code
compiles and runs.

    check_doc_snippets.py [--fix] [repo root]

exits 1 and prints each mismatch, or with --fix rewrites the blocks from the examples.
"""
import pathlib
import re
import sys
import textwrap

MARKER = re.compile(r"<!--\s*snippet:\s*(\S+)#(\S+)\s*-->\s*\n```(\w*)\n(.*?)```", re.S)


def region(path, tag):
    lines = path.read_text().splitlines()
    marks = [i for i, line in enumerate(lines) if line.strip() == f"//! [{tag}]"]
    if len(marks) != 2:
        raise SystemExit(f"{path}: expected two '//! [{tag}]' lines, found {len(marks)}")
    return textwrap.dedent("\n".join(lines[marks[0] + 1 : marks[1]])).strip("\n") + "\n"


def main(argv):
    fix = "--fix" in argv
    args = [a for a in argv if a != "--fix"]
    root = pathlib.Path(args[0] if args else pathlib.Path(__file__).resolve().parent.parent)
    bad = 0
    count = 0
    for md in sorted((root / "docs").rglob("*.md")):
        text = md.read_text()

        def replace(m):
            nonlocal bad, count
            count += 1
            want = region(root / m.group(1), m.group(2))
            if m.group(4) != want:
                if not fix:
                    bad += 1
                    print(f"{md.relative_to(root)}: block for {m.group(1)}#{m.group(2)} differs from the example")
                return m.group(0)[: m.start(4) - m.start(0)] + want + "```"
            return m.group(0)

        new = MARKER.sub(replace, text)
        if fix and new != text:
            md.write_text(new)
            print(f"updated {md.relative_to(root)}")
    if count == 0:
        raise SystemExit("no snippet markers found")
    if bad:
        print(f"{bad} of {count} snippets out of date; run tools/check_doc_snippets.py --fix")
        return 1
    print(f"{count} snippets match the examples")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
