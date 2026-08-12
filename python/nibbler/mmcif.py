"""Conventional PDBx/mmCIF and ModelCIF semantic API."""

from __future__ import annotations

from typing import Literal, NoReturn

from . import _core, cif
from ._input import Source
from ._objects import CifDocument, MmcifModel
from .cif import Destination, validate_profile
from .components import Registry
from .contracts import Diagnostic, Profile, Severity, ValidationReport
from .errors import ChemistryError, SchemaError

# The concise conventional spelling documented for model construction.
Model = MmcifModel


def read(
    source: Source | CifDocument,
    *,
    profile: Profile | str = Profile.PDBX,
    registry: Registry | None = None,
) -> MmcifModel:
    """Build an immutable PDBx or ModelCIF semantic model."""
    selected = _implemented_profile(profile)
    document = (
        source
        if isinstance(source, CifDocument)
        else cif.read(source, schema=selected.value)
    )
    if not isinstance(document, CifDocument):  # pragma: no cover - defensive narrowing
        raise TypeError("nibbler.mmcif.read() requires a full CIF document")
    try:
        native_registry = None if registry is None else registry._native
        native: _core._PdbxModel | _core._ModelCifModel
        if selected is Profile.MODELCIF:
            native = _core.build_modelcif_model(document._native, native_registry)
        else:
            native = _core.build_pdbx_model(document._native, native_registry)
        return MmcifModel(native, selected.value)
    except ValueError as error:
        _raise_model_error(error)


def validate(
    value: object,
    *,
    profile: Profile | str,
) -> ValidationReport:
    """Validate a document or semantic model against an explicit profile."""
    selected = _implemented_profile(profile)
    return _validate_value(value, selected)


def _validate_value(value: object, profile: Profile) -> ValidationReport:
    try:
        if isinstance(value, MmcifModel):
            if value.profile != profile.value:
                raise ValueError(
                    f"model profile is {value.profile!r}, not {profile.value!r}"
                )
            if isinstance(value._native, _core._ModelCifModel):
                fields = _core.validate_modelcif_model(value._native)
            else:
                fields = _core.validate_pdbx_model(value._native)
        elif isinstance(value, CifDocument):
            if profile is Profile.MODELCIF:
                fields = _core.validate_modelcif_document(value._native)
            else:
                fields = _core.validate_pdbx_document(value._native)
        else:
            raise TypeError(
                "nibbler.mmcif.validate() requires CifDocument or MmcifModel"
            )
    except ValueError as error:
        _raise_schema_error(error)
    return _validation_report(fields)


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
    selected = _implemented_profile(profile)
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
        value._source_document()
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


def _implemented_profile(profile: Profile | str) -> Profile:
    return validate_profile(profile)


def _validation_report(
    fields: tuple[str, str, list[str], list[tuple[str, str, str, list[str]]]],
) -> ValidationReport:
    schema_name, version, coverage, diagnostics = fields
    return ValidationReport(
        tuple(
            Diagnostic(
                code=code,
                severity=Severity(severity),
                message=message,
                context=tuple(context),
            )
            for code, severity, message, context in diagnostics
        ),
        schema=schema_name,
        dictionary_version=version,
        coverage=tuple(coverage),
    )


def _raise_model_error(error: ValueError) -> NoReturn:
    if len(error.args) != 3:
        raise error
    code, message, context = error.args
    raise ChemistryError(
        code=str(code),
        message=str(message),
        context=tuple(str(value) for value in context),
    ) from None


def _raise_schema_error(error: ValueError) -> NoReturn:
    if len(error.args) != 2:
        raise error
    code, message = error.args
    raise SchemaError(code=str(code), message=str(message)) from None
