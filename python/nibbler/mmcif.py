"""Conventional PDBx/mmCIF and ModelCIF semantic API."""

from __future__ import annotations

from typing import Literal, cast

from . import _core, cif
from ._input import Source
from ._native import raise_chemistry_error, raise_schema_error, validation_report
from ._objects import CifDocument, MmcifModel
from .cif import Destination
from .components import Registry
from .contracts import Profile, ValidationReport


def read(
    source: Source | CifDocument,
    *,
    profile: Profile | str = Profile.PDBX,
    registry: Registry | None = None,
) -> MmcifModel:
    """Build an immutable PDBx or ModelCIF semantic model."""
    selected = Profile(profile)
    document = (
        source
        if isinstance(source, CifDocument)
        else cast(CifDocument, cif.read(source, schema=selected.value))
    )
    try:
        native_registry = None if registry is None else registry._native
        native = _core.build_mmcif_model(
            document._native, selected.value, native_registry
        )
        return MmcifModel(native)
    except ValueError as error:
        raise_chemistry_error(error)


def validate(
    value: object,
    *,
    profile: Profile | str,
) -> ValidationReport:
    """Validate a document or semantic model against an explicit profile."""
    selected = Profile(profile)
    return _validate_value(value, selected)


def _validate_value(value: object, profile: Profile) -> ValidationReport:
    try:
        if isinstance(value, MmcifModel):
            if value.profile != profile.value:
                raise ValueError(
                    f"model profile is {value.profile!r}, not {profile.value!r}"
                )
            fields = _core.validate_mmcif_model(value._native)
        elif isinstance(value, CifDocument):
            fields = _core.validate_mmcif_document(value._native, profile.value)
        else:
            raise TypeError(
                "nibbler.mmcif.validate() requires CifDocument or MmcifModel"
            )
    except ValueError as error:
        raise_schema_error(error)
    return validation_report(fields)


def write(
    value: object,
    destination: Destination,
    *,
    profile: Profile | str,
    mode: Literal["canonical", "preserve"] = "canonical",
    validate: Literal["none", "syntax", "dictionary", "profile"] = "profile",
    mirror_local_qa_metric: int | None = None,
    format: Literal["cif", "bcif"] | None = None,
) -> None:
    """Write one complete semantic model with transactional path output."""
    selected = Profile(profile)
    if not isinstance(value, MmcifModel):
        raise TypeError("nibbler.mmcif.write() requires MmcifModel")
    if value.profile != selected.value:
        raise ValueError(f"model profile is {value.profile!r}, not {selected.value!r}")
    if mirror_local_qa_metric is not None:
        if selected is not Profile.MODELCIF:
            raise ValueError("mirror_local_qa_metric requires profile='modelcif'")
        if isinstance(mirror_local_qa_metric, bool) or not isinstance(
            mirror_local_qa_metric, int
        ):
            raise TypeError("mirror_local_qa_metric must be an integer metric ID")
    if mode not in {"canonical", "preserve"}:
        raise ValueError("mode must be 'canonical' or 'preserve'")
    if mode == "preserve" and mirror_local_qa_metric is not None:
        raise ValueError("mode='preserve' cannot mirror a QA metric into B factors")
    if validate not in {"none", "syntax", "dictionary", "profile"}:
        raise ValueError(
            "validate must be 'none', 'syntax', 'dictionary', or 'profile'"
        )
    if validate == "profile":
        _validate_value(value, selected).raise_for_errors()

    generic_validation: Literal["none", "syntax", "dictionary"]
    if validate == "profile":
        generic_validation = "dictionary"
    else:
        generic_validation = validate
    document = (
        CifDocument(value._native.source_document(), value.profile)
        if mode == "preserve"
        else value.to_document(mirror_local_qa_metric=mirror_local_qa_metric)
    )
    cif.write(
        document,
        destination,
        mode=mode,
        validate=generic_validation,
        format=format,
    )
