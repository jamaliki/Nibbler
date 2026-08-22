"""Reference engines used for correctness-first benchmark comparisons."""

from __future__ import annotations

from dataclasses import dataclass
from importlib import import_module
from pathlib import Path
from typing import Protocol, cast

ProjectedRows = tuple[tuple[str, ...], ...]


class ProjectionAdapter(Protocol):
    """One parse-and-project adapter used by the isolated benchmark worker."""

    def __call__(self, file: Path, columns: tuple[str, ...]) -> ProjectedRows: ...


class ArrowColumn(Protocol):
    """The tiny PyArrow surface required by the benchmark."""

    def to_pylist(self) -> list[object]: ...


class ArrowTable(Protocol):
    """The tiny PyArrow table surface required by the benchmark."""

    def column(self, name: str) -> ArrowColumn: ...


def project_with_nibbler(file: Path, columns: tuple[str, ...]) -> ProjectedRows:
    """Project atom_site through Nibbler and consume its Arrow stream."""
    import nibbler

    table = nibbler.chomp(file, category="atom_site", columns=columns)
    if not isinstance(table, nibbler.CifTable):
        raise TypeError("projected Nibbler read did not return CifTable")
    pyarrow = import_module("pyarrow")
    arrow_table = cast(
        ArrowTable, pyarrow.RecordBatchReader.from_stream(table).read_all()
    )
    values = [arrow_table.column(column).to_pylist() for column in columns]
    return tuple(
        tuple(str(value) for value in row) for row in zip(*values, strict=True)
    )


def project_with_gemmi(file: Path, columns: tuple[str, ...]) -> ProjectedRows:
    """Parse and extract atom_site columns with Gemmi."""
    gemmi = import_module("gemmi")

    block = gemmi.cif.read_file(str(file)).sole_block()
    tags = [f"_atom_site.{column}" for column in columns]
    return tuple(tuple(str(value) for value in row) for row in block.find(tags))


def project_with_biotite(file: Path, columns: tuple[str, ...]) -> ProjectedRows:
    """Parse and extract atom_site columns with Biotite."""
    module = import_module("biotite.structure.io.pdbx")
    category = module.CIFFile.read(file).block["atom_site"]
    values = [category[column].as_array() for column in columns]
    return tuple(
        tuple(str(value) for value in row) for row in zip(*values, strict=True)
    )


def project_with_biopython(file: Path, columns: tuple[str, ...]) -> ProjectedRows:
    """Parse and extract atom_site columns with Bio.PDB."""
    module = import_module("Bio.PDB.MMCIF2Dict")
    values_by_tag = module.MMCIF2Dict(str(file))
    values = [values_by_tag[f"_atom_site.{column}"] for column in columns]
    return tuple(
        tuple(str(value) for value in row) for row in zip(*values, strict=True)
    )


@dataclass(frozen=True, slots=True)
class Engine:
    """One benchmark engine and everything needed to invoke and describe it."""

    distribution: str
    purpose: str
    project: ProjectionAdapter


ENGINES = {
    "nibbler": Engine("nibbler-cif", "Nibbler candidate", project_with_nibbler),
    "gemmi": Engine("gemmi", "fast general CIF parsing", project_with_gemmi),
    "biotite": Engine(
        "biotite", "columnar structural-biology workflows", project_with_biotite
    ),
    "biopython": Engine(
        "biopython", "Bio.PDB compatibility workloads", project_with_biopython
    ),
}
