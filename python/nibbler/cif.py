"""Conventional generic CIF API."""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from os import PathLike
from typing import BinaryIO, Literal, NoReturn, TypeAlias

from .contracts import Profile
from .errors import FeatureUnavailableError

Source: TypeAlias = str | PathLike[str] | bytes | bytearray | memoryview | BinaryIO
Destination: TypeAlias = str | PathLike[str] | BinaryIO


def read(
    source: Source,
    *,
    category: str | None = None,
    columns: Sequence[str] | None = None,
    where: Mapping[str, object] | None = None,
    schema: str | None = None,
) -> NoReturn:
    """Parse one CIF source once the shared Rust parser is available."""
    del source, category, columns, where, schema
    raise FeatureUnavailableError("nibbler.cif.read", required_phase=2)


def scan(
    sources: Iterable[Source],
    *,
    category: str | None = None,
    columns: Sequence[str] | None = None,
    where: Mapping[str, object] | None = None,
    schema: str | None = None,
    workers: int | None = None,
    on_error: Literal["raise", "collect"] = "raise",
) -> NoReturn:
    """Parse many CIF sources once bounded native scanning is available."""
    del sources, category, columns, where, schema, workers, on_error
    raise FeatureUnavailableError("nibbler.cif.scan", required_phase=2)


def validate(value: object, *, schema: str | None = None) -> NoReturn:
    """Validate generic CIF syntax and dictionary constraints."""
    del value, schema
    raise FeatureUnavailableError("nibbler.cif.validate", required_phase=3)


def write(
    value: object,
    destination: Destination,
    *,
    mode: Literal["canonical", "preserve"] = "canonical",
    validate: Literal["none", "syntax", "dictionary"] = "syntax",
) -> NoReturn:
    """Write generic CIF once the shared serializer is available."""
    del value, destination, mode, validate
    raise FeatureUnavailableError("nibbler.cif.write", required_phase=3)


def validate_profile(profile: Profile | str) -> Profile:
    """Normalize a semantic profile without accepting invented values."""
    return profile if isinstance(profile, Profile) else Profile(profile)
