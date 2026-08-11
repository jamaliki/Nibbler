"""Nibbler's small themed facade over conventional CIF APIs."""

from __future__ import annotations

from typing import Literal, NoReturn

from . import cif, mmcif
from ._core import contract_version, native_version
from .cif import Destination
from .contracts import (
    Diagnostic,
    MissingKind,
    Profile,
    Severity,
    SourceSpan,
    ValidationReport,
)
from .errors import FeatureUnavailableError, NibblerError, ValidationError

__all__ = [
    "Diagnostic",
    "FeatureUnavailableError",
    "MissingKind",
    "NibblerError",
    "Profile",
    "Severity",
    "SourceSpan",
    "ValidationError",
    "ValidationReport",
    "chomp",
    "cif",
    "contract_version",
    "feast",
    "mmcif",
    "native_version",
    "sniff",
    "spit",
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
) -> NoReturn:
    """Validate generic CIF or an explicit macromolecular profile."""
    if profile is None:
        cif.validate(value, schema=schema)
    if schema is not None:
        raise TypeError("schema and profile are mutually exclusive")
    mmcif.validate(value, profile=profile)


def spit(
    value: object,
    destination: Destination,
    *,
    profile: Profile | str | None = None,
    mode: Literal["canonical", "preserve"] = "canonical",
    validate: Literal["none", "syntax", "dictionary", "profile"] = "syntax",
) -> NoReturn:
    """Write generic CIF or an explicit macromolecular profile."""
    if profile is None:
        if validate == "profile":
            raise ValueError("validate='profile' requires a profile")
        cif.write(
            value,
            destination,
            mode=mode,
            validate=validate,
        )
    mmcif.write(
        value,
        destination,
        profile=profile,
        mode=mode,
        validate=validate,
    )
