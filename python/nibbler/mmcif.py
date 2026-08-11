"""Conventional PDBx/mmCIF and ModelCIF API."""

from __future__ import annotations

from typing import Literal, NoReturn

from .cif import Destination, validate_profile
from .contracts import Profile
from .errors import FeatureUnavailableError


def validate(value: object, *, profile: Profile | str) -> NoReturn:
    """Validate a semantic macromolecular model against an explicit profile."""
    del value
    normalized = validate_profile(profile)
    raise FeatureUnavailableError(
        f"nibbler.mmcif.validate[{normalized.value}]", required_phase=4
    )


def write(
    value: object,
    destination: Destination,
    *,
    profile: Profile | str,
    mode: Literal["canonical", "preserve"] = "canonical",
    validate: Literal["none", "syntax", "dictionary", "profile"] = "profile",
) -> NoReturn:
    """Write a complete PDBx/mmCIF or ModelCIF model."""
    del value, destination, mode, validate
    normalized = validate_profile(profile)
    raise FeatureUnavailableError(
        f"nibbler.mmcif.write[{normalized.value}]", required_phase=4
    )
