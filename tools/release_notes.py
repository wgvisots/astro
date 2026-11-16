#!/usr/bin/env python3
###########################################################################
# release_notes.py
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
"""Write the notes of a GitHub release from CHANGELOG.md.

The version is the one in the project() call of CMakeLists.txt. The notes are that version's
CHANGELOG section, with links to files in the repository made absolute (a release page has no
repository directory to resolve them against), and the section's compare link at the end.

With --tag, the script also checks that the release is ready: the tag must be v<version>, and the
section must carry a date, not "unreleased". CI runs it with --tag on a release tag, and without
--tag on every other push, so a CHANGELOG the script cannot read is found before the tag is.

Usage: release_notes.py [--tag vX.Y.Z] [--repo owner/name] [--ref REF] [-o notes.md]
"""

import argparse
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def project_version(cmake: str) -> str:
    match = re.search(r"project\(\s*astro\s+VERSION\s+([0-9]+\.[0-9]+\.[0-9]+)", cmake)
    if not match:
        raise SystemExit("release_notes.py: no project(astro VERSION x.y.z) in CMakeLists.txt")
    return match.group(1)


def changelog_section(changelog: str, version: str) -> tuple[str, str, str | None]:
    """The heading's date (or other text), the section body, and the compare link, if listed."""
    heading = re.compile(r"^## \[" + re.escape(version) + r"\] - (.+)$", re.MULTILINE)
    match = heading.search(changelog)
    if not match:
        raise SystemExit(f"release_notes.py: CHANGELOG.md has no section '## [{version}] - <date>'")
    following = re.compile(r"^## \[", re.MULTILINE).search(changelog, match.end())
    body = changelog[match.end() : following.start() if following else len(changelog)]
    # The link definitions at the end of the file belong to no section; cut them off the last one.
    body = re.split(r"^\[[^\]]+\]: ", body, maxsplit=1, flags=re.MULTILINE)[0]
    link = re.search(r"^\[" + re.escape(version) + r"\]: (\S+)$", changelog, re.MULTILINE)
    return match.group(1).strip(), body.strip(), link.group(1) if link else None


def absolute_links(text: str, repo: str, ref: str) -> str:
    """Point Markdown links to repository files at https://github.com/<repo>/blob/<ref>/<path>."""

    def replace(match: re.Match) -> str:
        target = match.group(2)
        if re.match(r"[a-z]+:|#|/", target):
            return match.group(0)
        return f"{match.group(1)}(https://github.com/{repo}/blob/{ref}/{target})"

    return re.sub(r"(\]\s?)\(([^)\s]+)\)", replace, text)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--tag", help="the release tag; checks that the release is ready")
    parser.add_argument("--repo", default="wgvisots/astro", help="owner/name, for absolute links")
    parser.add_argument("--ref", help="the commit or tag that links point to (default: the tag, or main)")
    parser.add_argument("-o", "--output", help="output file (default: standard output)")
    args = parser.parse_args()

    version = project_version((ROOT / "CMakeLists.txt").read_text())
    date, body, compare = changelog_section((ROOT / "CHANGELOG.md").read_text(), version)

    if args.tag is not None:
        if args.tag != f"v{version}":
            raise SystemExit(
                f"release_notes.py: tag {args.tag} does not match the version in CMakeLists.txt ({version})"
            )
        if not re.fullmatch(r"\d{4}-\d{2}-\d{2}", date):
            raise SystemExit(
                f"release_notes.py: CHANGELOG.md section {version} is dated '{date}'; give it the release date"
            )

    ref = args.ref or args.tag or "main"
    notes = absolute_links(body, args.repo, ref)
    if compare:
        notes += f"\n\n**Full list of commits:** {compare}"
    notes += "\n"

    if args.output:
        Path(args.output).write_text(notes)
    else:
        sys.stdout.write(notes)
    return 0


if __name__ == "__main__":
    sys.exit(main())
