"""Verify that committed CIF fixtures match their manifest."""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path
from typing import TypedDict, cast

if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib

ROOT = Path(__file__).resolve().parents[1]
FIXTURE_ROOT = ROOT / "tests" / "fixtures"
MANIFEST_FILE = FIXTURE_ROOT / "manifest.toml"


class RawFixture(TypedDict):
    """Serialized fixture record."""

    id: str
    file: str
    sha256: str
    classification: str
    profiles: list[str]
    features: list[str]


class RawManifest(TypedDict):
    """Serialized corpus manifest."""

    format_version: int
    fixture: list[RawFixture]


def load_manifest() -> RawManifest:
    """Load the typed fixture manifest."""
    raw = tomllib.loads(MANIFEST_FILE.read_text(encoding="utf-8"))
    return cast(RawManifest, raw)


def fixture_digest(file: Path) -> str:
    """Return the fixture's lowercase SHA-256 digest."""
    return hashlib.sha256(file.read_bytes()).hexdigest()


def verify_corpus() -> tuple[str, ...]:
    """Return validation problems, or an empty tuple for a valid corpus."""
    manifest = load_manifest()
    problems: list[str] = []
    fixture_ids: set[str] = set()
    declared_files: set[Path] = set()

    if manifest["format_version"] != 1:
        problems.append(f"unsupported manifest version: {manifest['format_version']}")

    for entry in manifest["fixture"]:
        fixture_id = entry["id"]
        if fixture_id in fixture_ids:
            problems.append(f"duplicate fixture id: {fixture_id}")
        fixture_ids.add(fixture_id)

        file = FIXTURE_ROOT / entry["file"]
        declared_files.add(file)
        if not file.is_file():
            problems.append(f"missing fixture: {entry['file']}")
            continue

        actual_digest = fixture_digest(file)
        if actual_digest != entry["sha256"]:
            problems.append(
                f"digest mismatch for {entry['file']}: "
                f"expected {entry['sha256']}, found {actual_digest}"
            )

        if entry["classification"] not in {"valid", "invalid"}:
            problems.append(f"invalid classification for {fixture_id}")
        if not entry["profiles"]:
            problems.append(f"fixture has no profiles: {fixture_id}")
        if not entry["features"]:
            problems.append(f"fixture has no features: {fixture_id}")

    disk_files = set(FIXTURE_ROOT.glob("**/*.cif"))
    for undeclared_file in sorted(disk_files - declared_files):
        relative_file = undeclared_file.relative_to(FIXTURE_ROOT)
        problems.append(f"fixture missing from manifest: {relative_file}")

    return tuple(problems)


def main() -> int:
    """Run corpus verification from the command line."""
    problems = verify_corpus()
    if problems:
        for problem in problems:
            print(problem)
        return 1

    fixture_count = len(load_manifest()["fixture"])
    print(f"verified {fixture_count} fixtures against {MANIFEST_FILE.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
