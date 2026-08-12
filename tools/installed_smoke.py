"""Exercise an installed Nibbler wheel outside the source tree."""

from __future__ import annotations

import gzip
import io
from importlib.metadata import version

import pyarrow  # type: ignore[import-untyped]

import nibbler

SOURCE = b"""data_smoke
_entry.id smoke

loop_
_atom_site.id
_atom_site.label_comp_id
_atom_site.pdbx_PDB_model_num
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
1 ALA 1 1.25 2.5 3.75
2 GLY 2 ? . -4.0
"""


def project(source: bytes) -> nibbler.CifTable:
    """Project and type the smoke-test atom rows."""
    table = nibbler.chomp(
        source,
        category="atom_site",
        columns=[
            "id",
            "label_comp_id",
            "pdbx_PDB_model_num",
            "Cartn_x",
            "Cartn_y",
            "Cartn_z",
        ],
        schema="pdbx",
    )
    if not isinstance(table, nibbler.CifTable):
        raise RuntimeError("projection returned a document")
    return table


def arrow(source: bytes) -> pyarrow.Table:
    """Import one native table through the standard Arrow stream protocol."""
    return pyarrow.RecordBatchReader.from_stream(
        project(source).with_missing("columns")
    ).read_all()


def main() -> None:
    """Check import metadata and every shipped input/output encoding."""
    if nibbler.__version__ != version("nibbler-cif"):
        raise RuntimeError("native and package versions differ")

    document = nibbler.chomp(SOURCE)
    if not isinstance(document, nibbler.CifDocument) or document.block_count != 1:
        raise RuntimeError("full-document parse failed")
    canonical = document.to_canonical().encode()
    if not isinstance(nibbler.chomp(canonical), nibbler.CifDocument):
        raise RuntimeError("canonical round trip failed")
    if not isinstance(nibbler.chomp(gzip.compress(SOURCE)), nibbler.CifDocument):
        raise RuntimeError("gzip round trip failed")

    binary = io.BytesIO()
    nibbler.dump(document, binary, format="bcif")
    binary_table = arrow(binary.getvalue())
    text_table = arrow(SOURCE)
    if not binary_table.equals(text_table):
        raise RuntimeError("text and BinaryCIF projections differ")
    if text_table.column("pdbx_PDB_model_num").to_pylist() != [1, 2]:
        raise RuntimeError("integer schema typing failed")
    if text_table.column("Cartn_x").to_pylist() != [1.25, None]:
        raise RuntimeError("floating-point schema typing failed")
    if text_table.column("Cartn_x__missing_kind").to_pylist() != [0, 1]:
        raise RuntimeError("missing-state export failed")

    print(
        f"qualified nibbler-cif {nibbler.__version__} from {nibbler.__file__} "
        f"with pyarrow {pyarrow.__version__}"
    )


if __name__ == "__main__":
    main()
