"""Convert private native return values into public Python objects."""

from __future__ import annotations

from typing import NoReturn, TypeAlias, cast

from .contracts import Diagnostic, Severity, ValidationReport
from .errors import (
    BatchError,
    ChemistryError,
    ParseError,
    ProjectionError,
    SchemaError,
    WriteError,
)

ErrorFields: TypeAlias = tuple[
    str,
    str,
    str | None,
    int | None,
    int | None,
    int | None,
    int | None,
]
ValidationFields: TypeAlias = tuple[
    str,
    str,
    list[str],
    list[tuple[str, str, str, list[str]]],
]


def raise_read_error(error: ValueError) -> NoReturn:
    """Raise the public exception encoded by a native read failure."""
    if len(error.args) != 7:
        raise error
    code, message, source_name, line, column, byte_start, byte_end = error.args
    code = str(code)
    if code.startswith("CIF_SCHEMA_"):
        raise SchemaError(code=code, message=str(message)) from None
    exception_type = (
        ProjectionError if code.startswith("CIF_PROJECTION_") else ParseError
    )
    raise exception_type(
        code=code,
        message=str(message),
        source_name=None if source_name is None else str(source_name),
        line=cast(int | None, line),
        column=cast(int | None, column),
        byte_start=cast(int | None, byte_start),
        byte_end=cast(int | None, byte_end),
    ) from None


def batch_error(source_index: int, fields: ErrorFields) -> BatchError:
    """Build one public batch error from native fields."""
    code, message, source_name, line, column, byte_start, byte_end = fields
    return BatchError(
        source_index,
        code=code,
        message=message,
        source_name=source_name,
        line=line,
        column=column,
        byte_start=byte_start,
        byte_end=byte_end,
    )


def raise_schema_error(error: ValueError) -> NoReturn:
    """Raise a public schema error encoded by native code and message fields."""
    if len(error.args) != 2:
        raise error
    code, message = error.args
    raise SchemaError(code=str(code), message=str(message)) from None


def raise_chemistry_error(error: ValueError) -> NoReturn:
    """Raise a public chemistry error encoded by native fields."""
    if len(error.args) != 3:
        raise error
    code, message, context = error.args
    raise ChemistryError(
        code=str(code),
        message=str(message),
        context=tuple(str(value) for value in context),
    ) from None


def raise_write_error(error: ValueError) -> NoReturn:
    """Raise a public write error encoded by native code and message fields."""
    if len(error.args) != 2:
        raise error
    code, message = error.args
    raise WriteError(code=str(code), message=str(message)) from None


def validation_report(fields: ValidationFields) -> ValidationReport:
    """Convert native validation fields into an immutable public report."""
    schema, version, coverage, diagnostics = fields
    return ValidationReport(
        tuple(
            Diagnostic(code, Severity(severity), message, context=tuple(context))
            for code, severity, message, context in diagnostics
        ),
        schema=schema,
        dictionary_version=version,
        coverage=tuple(coverage),
    )
