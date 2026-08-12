"""Conventional generic CIF API."""

from __future__ import annotations

import gzip
import os
import tempfile
from collections.abc import Iterable, Mapping, Sequence
from os import PathLike
from pathlib import Path
from typing import BinaryIO, Literal, TypeAlias, overload

from . import _core
from ._core import CifDocument as CifDocument
from ._core import CifTable as CifTable
from ._input import (
    Source,
    _normalize_predicates,
    _normalize_source,
)
from ._native import (
    raise_read_error,
    raise_schema_error,
    raise_write_error,
    validation_report,
)
from ._scan import BatchDiagnostics as BatchDiagnostics
from ._scan import ScanResult as ScanResult
from .contracts import ValidationReport
from .errors import SchemaError, WriteError

Destination: TypeAlias = str | PathLike[str] | BinaryIO


@overload
def read(
    source: Source,
    *,
    category: None = None,
    columns: None = None,
    where: None = None,
    schema: str | None = None,
) -> CifDocument: ...


@overload
def read(
    source: Source,
    *,
    category: str,
    columns: Sequence[str] | None = None,
    where: Mapping[str, object] | None = None,
    schema: str | None = None,
) -> CifTable: ...


def read(
    source: Source,
    *,
    category: str | None = None,
    columns: Sequence[str] | None = None,
    where: Mapping[str, object] | None = None,
    schema: str | None = None,
) -> CifDocument | CifTable:
    """Parse one filesystem, bytes-like, or binary-file CIF source.

    With no category this returns :class:`CifDocument`; with one category it returns
    :class:`CifTable`. Gzip is detected by content. Parsing and decompression release
    the GIL.
    """
    if category is None and (columns is not None or where is not None):
        raise TypeError("columns and where require category")
    normalized_schema = None if schema is None else _normalize_schema(schema)

    normalized_columns = None if columns is None else list(columns)
    predicates = _normalize_predicates(where)
    source_name, content = _normalize_source(source)
    try:
        if content is None:
            native = _core.read_file(
                source_name,
                category,
                normalized_columns,
                predicates,
                normalized_schema,
            )
        else:
            native = _core.read_bytes(
                source_name,
                content,
                category,
                normalized_columns,
                predicates,
                normalized_schema,
            )
    except ValueError as error:
        raise_read_error(error)
    return native


def scan(
    sources: Iterable[Source],
    *,
    category: str | None = None,
    columns: Sequence[str] | None = None,
    where: Mapping[str, object] | None = None,
    schema: str | None = None,
    workers: int | None = None,
    on_error: Literal["raise", "collect"] = "raise",
) -> ScanResult:
    """Parse many sources through a bounded, deterministic native worker pool."""
    if category is None and (columns is not None or where is not None):
        raise TypeError("columns and where require category")
    normalized_schema = None if schema is None else _normalize_schema(schema)
    if on_error not in {"raise", "collect"}:
        raise ValueError("on_error must be 'raise' or 'collect'")
    worker_count = min(32, os.cpu_count() or 1) if workers is None else workers
    if not 1 <= worker_count <= 256:
        raise ValueError("workers must be between 1 and 256")
    return ScanResult(
        sources,
        category=category,
        columns=columns,
        predicates=_normalize_predicates(where),
        schema=normalized_schema,
        workers=worker_count,
        on_error=on_error,
    )


def validate(value: object, *, schema: str | None = None) -> ValidationReport:
    """Validate generic CIF syntax and optional dictionary constraints."""
    if not isinstance(value, CifDocument):
        raise TypeError("nibbler.cif.validate() requires CifDocument")
    selected = value.schema if schema is None else _normalize_schema(schema)
    if selected is None:
        return ValidationReport(coverage=("cif-1.1-syntax",))
    try:
        schema_name, version, coverage, fields = _core.validate_document(
            value, selected
        )
    except ValueError as error:
        raise_schema_error(error)
    return validation_report((schema_name, version, coverage, fields))


def write(
    value: object,
    destination: Destination,
    *,
    mode: Literal["canonical", "preserve"] = "canonical",
    validate: Literal["none", "syntax", "dictionary"] = "syntax",
    format: Literal["cif", "bcif"] | None = None,
) -> None:
    """Serialize CIF or BinaryCIF without leaving a partial path destination."""
    if not isinstance(value, CifDocument):
        raise TypeError("nibbler.cif.write() requires CifDocument")
    if mode not in {"canonical", "preserve"}:
        raise ValueError("mode must be 'canonical' or 'preserve'")
    if validate not in {"none", "syntax", "dictionary"}:
        raise ValueError("validate must be 'none', 'syntax', or 'dictionary'")
    if validate == "dictionary":
        if value.schema is None:
            raise SchemaError(
                code="CIF_SCHEMA_REQUIRED",
                message="validate='dictionary' requires a schema-attached document",
            )
        _validate_document(value).raise_for_errors()

    destination_file = _destination_file(destination)
    output_format = _output_format(destination_file, format)
    if mode == "preserve" and output_format == "bcif":
        raise ValueError("mode='preserve' requires format='cif'")
    payload = _document_payload(value, mode, output_format)
    if destination_file is None:
        _write_stream(destination, payload)
        return
    if destination_file.name.lower().endswith(".gz"):
        payload = gzip.compress(payload, mtime=0)
    _write_path_transactionally(destination_file, payload, value.schema, validate)


def _normalize_schema(schema: str) -> str:
    if not isinstance(schema, str):
        raise TypeError("schema must be a string")
    normalized = schema.lower()
    if normalized == "mmcif":
        return "pdbx"
    if normalized not in {"pdbx", "modelcif"}:
        raise SchemaError(
            code="CIF_SCHEMA_UNKNOWN",
            message=f"unknown schema {schema!r}; expected 'pdbx' or 'modelcif'",
        )
    return normalized


def _document_payload(
    document: CifDocument,
    mode: Literal["canonical", "preserve"],
    output_format: Literal["cif", "bcif"],
) -> bytes:
    try:
        if output_format == "bcif":
            return document.to_binary()
        text = (
            document.to_preserving() if mode == "preserve" else document.to_canonical()
        )
        return text.encode("utf-8")
    except ValueError as error:
        raise_write_error(error)


def _output_format(
    destination: Path | None,
    requested: Literal["cif", "bcif"] | None,
) -> Literal["cif", "bcif"]:
    if requested is not None:
        if requested not in {"cif", "bcif"}:
            raise ValueError("format must be 'cif' or 'bcif'")
        return requested
    if destination is not None:
        name = destination.name.lower()
        if name.endswith(".bcif") or name.endswith(".bcif.gz"):
            return "bcif"
    return "cif"


def _destination_file(destination: Destination) -> Path | None:
    if not isinstance(destination, (str, os.PathLike)):
        return None
    file = os.fspath(destination)
    if not isinstance(file, str):
        raise TypeError("filesystem destinations must resolve to a string path")
    return Path(file)


def _write_stream(destination: Destination, payload: bytes) -> None:
    writer = getattr(destination, "write", None)
    if writer is None:
        raise TypeError("destination must be a path or writable binary stream")
    try:
        written = writer(payload)
    except OSError as error:
        raise WriteError(code="CIF_WRITE_IO", message=str(error)) from error
    if written is not None and written != len(payload):
        raise WriteError(
            code="CIF_WRITE_SHORT",
            message=f"binary stream accepted {written} of {len(payload)} bytes",
        )


def _write_path_transactionally(
    destination: Path,
    payload: bytes,
    schema: str | None,
    validation: Literal["none", "syntax", "dictionary"],
) -> None:
    parent = destination.parent or Path(".")
    descriptor = -1
    temporary: Path | None = None
    try:
        descriptor, temporary_name = tempfile.mkstemp(
            dir=parent, prefix=f".{destination.name}.", suffix=".tmp"
        )
        temporary = Path(temporary_name)
        with os.fdopen(descriptor, "wb") as output:
            descriptor = -1
            output.write(payload)
            output.flush()
            os.fsync(output.fileno())
        if validation != "none":
            reparsed = read(temporary)
            if not isinstance(reparsed, CifDocument):
                raise WriteError(
                    code="CIF_WRITE_REPARSE",
                    message="transactional output did not reparse as a document",
                    destination=str(destination),
                )
            if validation == "dictionary":
                _validate_document(reparsed, schema=schema).raise_for_errors()
        os.replace(temporary, destination)
        temporary = None
    except OSError as error:
        raise WriteError(
            code="CIF_WRITE_IO",
            message=str(error),
            destination=str(destination),
        ) from error
    finally:
        if descriptor >= 0:
            os.close(descriptor)
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def _validate_document(
    document: CifDocument, *, schema: str | None = None
) -> ValidationReport:
    """Validate from write(), where the public `validate` parameter shadows the API."""
    return validate(document, schema=schema)
