"""Typed access to the pinned real-PDB benchmark corpus."""

from __future__ import annotations

import hashlib
import sys
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

if sys.version_info >= (3, 11):
    import tomllib
else:  # pragma: no cover - Python 3.10
    import tomli as tomllib

DEFAULT_MANIFEST = Path("benchmarks/pdb_corpus.toml")


@dataclass(frozen=True, slots=True)
class CorpusFile:
    """One hash-pinned corpus artifact."""

    bytes: int
    sha256: str


@dataclass(frozen=True, slots=True)
class RemoteCorpusFile(CorpusFile):
    """A corpus artifact fetched directly rather than derived locally."""

    url: str


@dataclass(frozen=True, slots=True)
class Structure:
    """One benchmark structure and its equivalent encodings."""

    pdb_id: str
    workload: str
    atom_rows: int
    cif: CorpusFile
    bcif: RemoteCorpusFile
    cif_gz: RemoteCorpusFile

    def file(self, file_format: str) -> CorpusFile:
        return {"cif": self.cif, "bcif": self.bcif, "cif.gz": self.cif_gz}[file_format]


def load_structures(manifest: Path = DEFAULT_MANIFEST) -> tuple[Structure, ...]:
    """Load and type-check every structure in a corpus manifest."""
    with manifest.open("rb") as stream:
        document: dict[str, Any] = tomllib.load(stream)
    raw_structures = document.get("structures")
    if not isinstance(raw_structures, list) or not raw_structures:
        raise RuntimeError(f"{manifest}: missing [[structures]] entries")
    return tuple(_structure(raw) for raw in raw_structures)


def select_structures(
    structures: Sequence[Structure], selected_ids: Sequence[str]
) -> tuple[Structure, ...]:
    """Select requested IDs in manifest order and reject unknown IDs."""
    selected = set(selected_ids)
    unknown = selected - {structure.pdb_id for structure in structures}
    if unknown:
        raise RuntimeError(f"unknown PDB IDs: {', '.join(sorted(unknown))}")
    return tuple(
        structure
        for structure in structures
        if not selected or structure.pdb_id in selected
    )


def verify_content(content: bytes, expected: CorpusFile, label: str) -> bytes:
    """Return content after checking its pinned size and digest."""
    if len(content) != expected.bytes:
        raise RuntimeError(
            f"{label}: expected {expected.bytes} bytes, received {len(content)}"
        )
    digest = hashlib.sha256(content).hexdigest()
    if digest != expected.sha256:
        raise RuntimeError(
            f"{label}: expected sha256 {expected.sha256}, received {digest}"
        )
    return content


def verify_file(file: Path, expected: CorpusFile) -> None:
    """Check one local corpus file against its manifest record."""
    verify_content(file.read_bytes(), expected, file.name)


def _structure(raw: object) -> Structure:
    if not isinstance(raw, dict):
        raise RuntimeError("each [[structures]] entry must be a table")

    def text(name: str) -> str:
        value = raw.get(name)
        if not isinstance(value, str):
            raise RuntimeError(f"corpus field {name!r} must be a string")
        return value

    def integer(name: str) -> int:
        value = raw.get(name)
        if isinstance(value, bool) or not isinstance(value, int):
            raise RuntimeError(f"corpus field {name!r} must be an integer")
        return value

    def file(prefix: str) -> CorpusFile:
        return CorpusFile(
            bytes=integer(f"{prefix}_bytes"),
            sha256=text(f"{prefix}_sha256"),
        )

    def remote_file(prefix: str) -> RemoteCorpusFile:
        return RemoteCorpusFile(
            bytes=integer(f"{prefix}_bytes"),
            sha256=text(f"{prefix}_sha256"),
            url=text(f"{prefix}_url"),
        )

    return Structure(
        pdb_id=text("id"),
        workload=text("workload"),
        atom_rows=integer("atom_rows"),
        cif=file("cif"),
        bcif=remote_file("bcif"),
        cif_gz=remote_file("cif_gz"),
    )
