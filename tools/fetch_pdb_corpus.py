"""Fetch and verify the pinned real-PDB stress corpus."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import os
import sys
import tempfile
import urllib.request
from collections.abc import Mapping, Sequence
from pathlib import Path

if sys.version_info >= (3, 11):
    import tomllib
else:  # pragma: no cover - Python 3.10
    import tomli as tomllib

MANIFEST = Path("benchmarks/pdb_corpus.toml")
DEFAULT_DESTINATION = Path(".cache/pdb-stress")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=MANIFEST)
    parser.add_argument("--destination", type=Path, default=DEFAULT_DESTINATION)
    parser.add_argument("--ids", nargs="*", default=[])
    return parser.parse_args(argv)


def sha256(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def checked_content(content: bytes, size: int, digest: str, label: str) -> bytes:
    if len(content) != size:
        raise RuntimeError(f"{label}: expected {size} bytes, received {len(content)}")
    actual = sha256(content)
    if actual != digest:
        raise RuntimeError(f"{label}: expected sha256 {digest}, received {actual}")
    return content


def download(url: str) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": "Nibbler/0.1"})
    with urllib.request.urlopen(request, timeout=120) as response:
        content = response.read()
    if not isinstance(content, bytes):
        raise RuntimeError(f"{url}: download did not return bytes")
    return content


def atomic_write(file: Path, content: bytes) -> None:
    file.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=file.parent, prefix=f".{file.name}."
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as output:
            output.write(content)
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary, file)
    finally:
        temporary.unlink(missing_ok=True)


def load_structures(manifest: Path) -> list[Mapping[str, object]]:
    with manifest.open("rb") as stream:
        document = tomllib.load(stream)
    structures = document.get("structures")
    if not isinstance(structures, list):
        raise RuntimeError(f"{manifest}: missing [[structures]] entries")
    return structures


def string_field(structure: Mapping[str, object], name: str) -> str:
    value = structure.get(name)
    if not isinstance(value, str):
        raise RuntimeError(f"corpus field {name!r} must be a string")
    return value


def integer_field(structure: Mapping[str, object], name: str) -> int:
    value = structure.get(name)
    if isinstance(value, bool) or not isinstance(value, int):
        raise RuntimeError(f"corpus field {name!r} must be an integer")
    return value


def fetch_structure(structure: Mapping[str, object], destination: Path) -> None:
    pdb_id = string_field(structure, "id")
    cif_gz = checked_content(
        download(string_field(structure, "cif_gz_url")),
        integer_field(structure, "cif_gz_bytes"),
        string_field(structure, "cif_gz_sha256"),
        f"{pdb_id}.cif.gz",
    )
    cif = checked_content(
        gzip.decompress(cif_gz),
        integer_field(structure, "cif_bytes"),
        string_field(structure, "cif_sha256"),
        f"{pdb_id}.cif",
    )
    bcif = checked_content(
        download(string_field(structure, "bcif_url")),
        integer_field(structure, "bcif_bytes"),
        string_field(structure, "bcif_sha256"),
        f"{pdb_id}.bcif",
    )
    atomic_write(destination / f"{pdb_id}.cif.gz", cif_gz)
    atomic_write(destination / f"{pdb_id}.cif", cif)
    atomic_write(destination / f"{pdb_id}.bcif", bcif)
    print(f"verified {pdb_id}: {len(cif):,} CIF bytes, {len(bcif):,} BinaryCIF bytes")


def main(argv: Sequence[str] | None = None) -> int:
    arguments = parse_args(argv)
    selected = set(arguments.ids)
    structures = load_structures(arguments.manifest)
    known = {string_field(structure, "id") for structure in structures}
    unknown = selected - known
    if unknown:
        raise RuntimeError(f"unknown PDB IDs: {', '.join(sorted(unknown))}")
    for structure in structures:
        pdb_id = string_field(structure, "id")
        if not selected or pdb_id in selected:
            fetch_structure(structure, arguments.destination)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
