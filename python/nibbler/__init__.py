"""Nibbler's small themed facade over conventional CIF APIs."""

from __future__ import annotations

from typing import Literal

from . import cif, components, mmcif, reduce
from ._core import MmcifModel, contract_version, native_version
from .cif import BatchDiagnostics, CifDocument, CifTable, Destination, ScanResult
from .contracts import (
    Diagnostic,
    MissingKind,
    Profile,
    Severity,
    SourceSpan,
    ValidationReport,
)
from .errors import (
    BatchError,
    ChemistryError,
    NibblerError,
    ParseError,
    ProjectionError,
    SchemaError,
    ValidationError,
    WriteError,
)

__all__ = [
    "BatchDiagnostics",
    "BatchError",
    "ChemistryError",
    "CifDocument",
    "CifTable",
    "Diagnostic",
    "MissingKind",
    "MmcifModel",
    "NibblerError",
    "ParseError",
    "Profile",
    "ProjectionError",
    "ScanResult",
    "SchemaError",
    "Severity",
    "SourceSpan",
    "ValidationError",
    "ValidationReport",
    "WriteError",
    "chomp",
    "cif",
    "components",
    "contract_version",
    "dump",
    "feast",
    "mmcif",
    "native_version",
    "reduce",
    "sniff",
]

__version__ = native_version()

# These two names are exact aliases, not alternate implementations.
chomp = cif.read
feast = cif.scan


def sniff(
    value: object,
    *,
    profile: Profile | str | None = None,
    schema: str | None = None,
) -> ValidationReport:
    """Validate generic CIF or an explicit macromolecular profile."""
    if profile is None:
        return cif.validate(value, schema=schema)
    if schema is not None:
        raise TypeError("schema and profile are mutually exclusive")
    return mmcif.validate(value, profile=profile)


def dump(
    value: object,
    destination: Destination,
    *,
    profile: Profile | str | None = None,
    mode: Literal["canonical", "preserve"] = "canonical",
    validate: Literal["none", "syntax", "dictionary", "profile"] | None = None,
    mirror_local_qa_metric: int | None = None,
    format: Literal["cif", "bcif"] | None = None,
) -> None:
    """Write generic CIF or an explicit macromolecular profile."""
    if profile is None:
        if mirror_local_qa_metric is not None:
            raise ValueError("mirror_local_qa_metric requires profile='modelcif'")
        generic_validation = "syntax" if validate is None else validate
        if generic_validation == "profile":
            raise ValueError("validate='profile' requires a profile")
        return cif.write(
            value,
            destination,
            mode=mode,
            validate=generic_validation,
            format=format,
        )
    return mmcif.write(
        value,
        destination,
        profile=profile,
        mode=mode,
        validate="profile" if validate is None else validate,
        mirror_local_qa_metric=mirror_local_qa_metric,
        format=format,
    )
