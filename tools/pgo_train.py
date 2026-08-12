"""Train PGO with ``python -m tools.pgo_train`` from the repository root."""

from __future__ import annotations

import argparse
import json
from collections.abc import Sequence
from pathlib import Path

from benchmarks.corpus import load_structures, verify_file

PROJECTED_COLUMNS = ("label_comp_id", "Cartn_x", "Cartn_y", "Cartn_z")
FULL_DOCUMENT_IDS = frozenset(("1crn", "6qnr"))


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


def main(argv: Sequence[str] | None = None) -> None:
    arguments = parse_args(argv)

    import nibbler
    import nibbler._core

    trained: list[dict[str, object]] = []
    for structure in load_structures(arguments.manifest):
        pdb_id = structure.pdb_id
        for suffix in ("cif", "bcif"):
            file = arguments.corpus / f"{pdb_id}.{suffix}"
            verify_file(file, structure.file(suffix))
            project(file, structure.atom_rows)
            trained.append({"file": file.name, "operation": "projection"})
        gzip_file = arguments.corpus / f"{pdb_id}.cif.gz"
        if gzip_file.is_file():
            verify_file(gzip_file, structure.cif_gz)
            project(gzip_file, structure.atom_rows)
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
