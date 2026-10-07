"""Conventional PDBx/mmCIF and ModelCIF semantic API."""

from __future__ import annotations

from typing import Literal

from . import _core, cif
from ._core import CifDocument, MmcifModel
from ._input import Source
from ._native import (
    raise_chemistry_error,
    raise_schema_error,
    raise_write_error,
    validation_report,
)
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
        else cif.read(source, schema=selected.value)
    )
    try:
        return _core.build_mmcif_model(document, selected.value, registry)
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
            fields = _core.validate_mmcif_model(value)
        elif isinstance(value, CifDocument):
            fields = _core.validate_mmcif_document(value, profile.value)
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
    try:
        document = (
            value.source_document()
            if mode == "preserve"
            else _core._model_document(value, mirror_local_qa_metric)
        )
    except ValueError as error:
        raise_write_error(error)
    cif.write(
        document,
        destination,
        mode=mode,
        validate=generic_validation,
        format=format,
    )


def assembly(
    source: Source | CifDocument | MmcifModel,
    assembly_id: str | None = None,
    *,
    registry: Registry | None = None,
) -> CifDocument:
    """Write out every copy of one biological assembly as explicit coordinates.

    ``assembly_id`` selects a ``_pdbx_struct_assembly`` (default: the first). Each
    distinct operator combination of its ``_pdbx_struct_assembly_gen`` rows is one copy;
    copy 1 keeps the chain identifiers and copy *n* appends ``-n`` to the
    ``label_asym_id`` and ``auth_asym_id`` of its chains. An atom on a symmetry axis of
    the assembly, whose copies coincide, is written once. The result holds the entry,
    entity, and component categories, and one row per copy in ``_struct_asym``, the
    sequence schemes, ``_atom_site`` (coordinates transformed, atom ids renumbered),
    ``_atom_site_anisotrop`` (tensors rotated), and ``_struct_conn`` (connections
    within one crystal copy). The crystal cell, symmetry, and assembly definitions are
    left out, so nothing downstream applies the crystal symmetry again.
    """
    if assembly_id is not None and not isinstance(assembly_id, str):
        raise TypeError("assembly_id must be a string")
    model = (
        source
        if isinstance(source, MmcifModel)
        else read(source, profile=Profile.PDBX, registry=registry)
    )
    try:
        return _core.build_assembly_document(model, assembly_id)
    except ValueError as error:
        raise_chemistry_error(error)
