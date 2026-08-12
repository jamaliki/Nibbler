"""Fetch the pinned PDB corpus with ``python -m tools.fetch_pdb_corpus``."""

from __future__ import annotations

import argparse
import gzip
import os
import tempfile
import urllib.request
from collections.abc import Sequence
from pathlib import Path

from benchmarks.corpus import (
    DEFAULT_MANIFEST,
    Structure,
    load_structures,
    select_structures,
    verify_content,
)

DEFAULT_DESTINATION = Path(".cache/pdb-stress")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--destination", type=Path, default=DEFAULT_DESTINATION)
    parser.add_argument("--ids", nargs="*", default=[])
    return parser.parse_args(argv)


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


def fetch_structure(structure: Structure, destination: Path) -> None:
    pdb_id = structure.pdb_id
    cif_gz = verify_content(
        download(structure.cif_gz.url), structure.cif_gz, f"{pdb_id}.cif.gz"
    )
    cif = verify_content(
        gzip.decompress(cif_gz),
        structure.cif,
        f"{pdb_id}.cif",
    )
    bcif = verify_content(
        download(structure.bcif.url),
        structure.bcif,
        f"{pdb_id}.bcif",
    )
    atomic_write(destination / f"{pdb_id}.cif.gz", cif_gz)
    atomic_write(destination / f"{pdb_id}.cif", cif)
    atomic_write(destination / f"{pdb_id}.bcif", bcif)
    print(f"verified {pdb_id}: {len(cif):,} CIF bytes, {len(bcif):,} BinaryCIF bytes")


def main(argv: Sequence[str] | None = None) -> int:
    arguments = parse_args(argv)
    structures = select_structures(load_structures(arguments.manifest), arguments.ids)
    for structure in structures:
        fetch_structure(structure, arguments.destination)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
