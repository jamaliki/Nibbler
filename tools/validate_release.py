"""Validate Nibbler release metadata and distribution archives."""

from __future__ import annotations

import argparse
import hashlib
import re
import sys
import tarfile
import zipfile
from collections import Counter
from collections.abc import Iterable, Sequence
from email.parser import BytesParser
from email.policy import default
from pathlib import Path

if sys.version_info >= (3, 11):
    import tomllib
else:  # pragma: no cover - Python 3.10
    import tomli as tomllib

ROOT = Path(__file__).resolve().parents[1]
PROJECT_FILE = ROOT / "pyproject.toml"
CARGO_FILE = ROOT / "Cargo.toml"
CHANGELOG_FILE = ROOT / "CHANGELOG.md"
PACKAGE_NAME = "nibbler-cif"
WHEEL_STEM = "nibbler_cif"
REQUIRES_PYTHON = ">=3.10"


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse metadata-only or archive-validation options."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dist", type=Path, default=Path("dist"))
    parser.add_argument("--metadata-only", action="store_true")
    parser.add_argument("--tag")
    parser.add_argument("--expected-wheel-count", type=int)
    parser.add_argument("--python-tags", nargs="*", default=[])
    parser.add_argument("--platform-count", type=int)
    parser.add_argument("--write-checksums", type=Path)
    return parser.parse_args(argv)


def load_toml(file: Path) -> dict[str, object]:
    """Load one TOML document with Python 3.10 compatibility."""
    with file.open("rb") as stream:
        return tomllib.load(stream)


def project_version(tag: str | None = None) -> str:
    """Return the single version shared by Python, Rust, tags, and the changelog."""
    project = load_toml(PROJECT_FILE)["project"]
    cargo = load_toml(CARGO_FILE)["package"]
    if not isinstance(project, dict) or not isinstance(cargo, dict):
        raise RuntimeError("package metadata tables are malformed")
    python_version = project.get("version")
    rust_version = cargo.get("version")
    if not isinstance(python_version, str) or python_version != rust_version:
        raise RuntimeError(
            f"Python version {python_version!r} and Rust version "
            f"{rust_version!r} differ"
        )
    if re.fullmatch(r"[0-9]+\.[0-9]+\.[0-9]+(?:rc[0-9]+)?", python_version) is None:
        raise RuntimeError(f"release version is not supported: {python_version!r}")
    if tag is not None and tag.removeprefix("refs/tags/") != f"v{python_version}":
        raise RuntimeError(f"tag {tag!r} does not match version {python_version}")
    changelog_header = next(
        (
            line
            for line in CHANGELOG_FILE.read_text().splitlines()
            if line.startswith(f"## [{python_version}]")
        ),
        None,
    )
    if changelog_header is None:
        raise RuntimeError(f"CHANGELOG.md has no {python_version} release entry")
    if tag is not None and "Unreleased" in changelog_header:
        raise RuntimeError(f"CHANGELOG.md entry {changelog_header!r} is not dated")
    return python_version


def metadata(content: bytes, label: str, version: str) -> None:
    """Validate the shared Core Metadata fields in one archive."""
    message = BytesParser(policy=default).parsebytes(content)
    expected = {
        "Name": PACKAGE_NAME,
        "Version": version,
        "Requires-Python": REQUIRES_PYTHON,
    }
    for field, value in expected.items():
        if message[field] != value:
            raise RuntimeError(
                f"{label}: expected {field} {value!r}, found {message[field]!r}"
            )


def one_matching(names: Iterable[str], suffix: str, label: str) -> str:
    """Return exactly one archive member ending with a required suffix."""
    matches = [name for name in names if name.endswith(suffix)]
    if len(matches) != 1:
        raise RuntimeError(f"{label}: expected one {suffix}, found {len(matches)}")
    return matches[0]


def validate_wheel(file: Path, version: str) -> tuple[str, str]:
    """Validate one wheel and return its Python and platform tags."""
    pattern = re.compile(
        rf"^{WHEEL_STEM}-{re.escape(version)}-(cp[0-9]+)-(cp[0-9]+)-(.+)\.whl$"
    )
    match = pattern.fullmatch(file.name)
    if match is None or match.group(1) != match.group(2):
        raise RuntimeError(f"unexpected interpreter-specific wheel name: {file.name}")
    python_tag, _, platform_tag = match.groups()
    with zipfile.ZipFile(file) as archive:
        names = archive.namelist()
        metadata_name = one_matching(names, ".dist-info/METADATA", file.name)
        metadata(archive.read(metadata_name), file.name, version)
        required = (
            "nibbler/__init__.py",
            "nibbler/_core.pyi",
            "nibbler/py.typed",
        )
        missing = [name for name in required if name not in names]
        if missing:
            raise RuntimeError(f"{file.name}: missing {', '.join(missing)}")
        native = [
            name
            for name in names
            if name.startswith("nibbler/_core.")
            and (name.endswith(".so") or name.endswith(".pyd"))
        ]
        if len(native) != 1:
            raise RuntimeError(
                f"{file.name}: expected one native module, found {native}"
            )
        if not any(name.endswith(".dist-info/licenses/LICENSE") for name in names):
            raise RuntimeError(f"{file.name}: LICENSE is missing")
        forbidden = [
            name for name in names if "__pycache__" in name or name.endswith(".pyc")
        ]
        if forbidden:
            raise RuntimeError(f"{file.name}: generated Python caches are present")
    return python_tag, platform_tag


def validate_sdist(file: Path, version: str) -> None:
    """Validate one source distribution's metadata and build inputs."""
    with tarfile.open(file, "r:gz") as archive:
        names = archive.getnames()
        roots = {name.split("/", 1)[0] for name in names}
        if len(roots) != 1:
            raise RuntimeError(f"{file.name}: source archive has multiple roots")
        root = roots.pop()
        required = (
            "Cargo.toml",
            "Cargo.lock",
            "pyproject.toml",
            "LICENSE",
            "README.md",
            "CHANGELOG.md",
            "src/lib.rs",
            "python/nibbler/__init__.py",
            "python/nibbler/_core.pyi",
            "python/nibbler/py.typed",
            "schemas/compiled/pdbx.nbs",
            "schemas/compiled/modelcif.nbs",
        )
        missing = [name for name in required if f"{root}/{name}" not in names]
        if missing:
            raise RuntimeError(f"{file.name}: missing {', '.join(missing)}")
        member = archive.getmember(f"{root}/PKG-INFO")
        stream = archive.extractfile(member)
        if stream is None:
            raise RuntimeError(f"{file.name}: PKG-INFO is unreadable")
        metadata(stream.read(), file.name, version)


def write_checksums(files: Sequence[Path], destination: Path) -> None:
    """Write deterministic SHA-256 checksums for qualified archives."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        f"{hashlib.sha256(file.read_bytes()).hexdigest()}  {file.name}"
        for file in sorted(files)
    ]
    destination.write_text("\n".join(lines) + "\n")


def validate_archives(arguments: argparse.Namespace, version: str) -> None:
    """Validate archive counts, contents, tags, and optional checksums."""
    distribution = arguments.dist.resolve()
    wheels = sorted(distribution.glob("*.whl"))
    source_distributions = sorted(distribution.glob("*.tar.gz"))
    if (
        arguments.expected_wheel_count is not None
        and len(wheels) != arguments.expected_wheel_count
    ):
        raise RuntimeError(
            f"expected {arguments.expected_wheel_count} wheels, found {len(wheels)}"
        )
    if len(source_distributions) != 1:
        raise RuntimeError(
            f"expected one source distribution, found {len(source_distributions)}"
        )
    tags = [validate_wheel(file, version) for file in wheels]
    validate_sdist(source_distributions[0], version)
    if arguments.python_tags:
        counts = Counter(python_tag for python_tag, _ in tags)
        expected = Counter(
            {
                python_tag: arguments.platform_count
                for python_tag in arguments.python_tags
            }
        )
        if counts != expected:
            raise RuntimeError(
                f"wheel Python tags differ: expected {expected}, found {counts}"
            )
    if len({tag for tag in tags}) != len(tags):
        raise RuntimeError("duplicate Python/platform wheel tags")
    files = [*wheels, *source_distributions]
    if arguments.write_checksums is not None:
        write_checksums(files, arguments.write_checksums.resolve())
    print(
        f"qualified nibbler-cif {version}: {len(wheels)} wheel(s), "
        f"{len(source_distributions)} source distribution"
    )


def main(argv: Sequence[str] | None = None) -> None:
    """Validate synchronized metadata and, unless excluded, distribution archives."""
    arguments = parse_args(argv)
    version = project_version(arguments.tag)
    if arguments.metadata_only:
        print(version)
        return
    validate_archives(arguments, version)


if __name__ == "__main__":
    main()
