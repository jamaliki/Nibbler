"""Exercise representative Nibbler paths in an instrumented PGO wheel."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any

if sys.version_info >= (3, 11):
    import tomllib
else:  # pragma: no cover - Python 3.10
    import tomli as tomllib

PROJECTED_COLUMNS = ("label_comp_id", "Cartn_x", "Cartn_y", "Cartn_z")
FULL_DOCUMENT_IDS = frozenset(("1crn", "6qnr"))


def load_manifest(file: Path) -> list[dict[str, Any]]:
    with file.open("rb") as stream:
        contents: dict[str, Any] = tomllib.load(stream)
    structures = contents.get("structures")
    if not isinstance(structures, list) or not structures:
        raise RuntimeError(f"{file}: missing [[structures]] entries")
    return structures


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    return parser.parse_args(argv)


def project(file: Path, expected_rows: int) -> None:
    import nibbler

    table = nibbler.chomp(
        file,
        category="atom_site",
        columns=PROJECTED_COLUMNS,
        schema="pdbx",
    )
    if not isinstance(table, nibbler.CifTable):
        raise RuntimeError(f"{file}: projection returned a document")
    if len(table) != expected_rows:
        raise RuntimeError(f"{file}: expected {expected_rows} rows, found {len(table)}")


def verify_file(file: Path, expected_bytes: object, expected_sha256: object) -> None:
    if not isinstance(expected_bytes, int) or not isinstance(expected_sha256, str):
        raise RuntimeError(f"{file}: manifest size and digest must be typed")
    content = file.read_bytes()
    if len(content) != expected_bytes:
        raise RuntimeError(
            f"{file}: expected {expected_bytes} bytes, found {len(content)}"
        )
    digest = hashlib.sha256(content).hexdigest()
    if digest != expected_sha256:
        raise RuntimeError(f"{file}: expected sha256 {expected_sha256}, found {digest}")


def main(argv: Sequence[str] | None = None) -> None:
    import nibbler
    import nibbler._core

    arguments = parse_args(argv)
    trained: list[dict[str, object]] = []
    for structure in load_manifest(arguments.manifest):
        pdb_id = str(structure["id"])
        expected_rows = int(structure["atom_rows"])
        for suffix in ("cif", "bcif"):
            file = arguments.corpus / f"{pdb_id}.{suffix}"
            verify_file(
                file, structure[f"{suffix}_bytes"], structure[f"{suffix}_sha256"]
            )
            project(file, expected_rows)
            trained.append({"file": file.name, "operation": "projection"})
        gzip_file = arguments.corpus / f"{pdb_id}.cif.gz"
        if gzip_file.is_file():
            verify_file(
                gzip_file,
                structure["cif_gz_bytes"],
                structure["cif_gz_sha256"],
            )
            project(gzip_file, expected_rows)
            trained.append({"file": gzip_file.name, "operation": "projection"})
        if pdb_id in FULL_DOCUMENT_IDS:
            for suffix in ("cif", "bcif"):
                file = arguments.corpus / f"{pdb_id}.{suffix}"
                document = nibbler.chomp(file)
                if not isinstance(document, nibbler.CifDocument):
                    raise RuntimeError(f"{file}: full-document construction failed")
                trained.append({"file": file.name, "operation": "full_document"})
    print(json.dumps({"extension": nibbler._core.__file__, "trained": trained}))


if __name__ == "__main__":
    main()
