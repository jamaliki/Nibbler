"""Public Nibbler exception hierarchy."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .contracts import Diagnostic


class NibblerError(Exception):
    """Base class for errors raised by Nibbler."""


class FeatureUnavailableError(NibblerError):
    """Raised when an API belongs to a later implementation phase."""

    operation: str
    required_phase: int

    def __init__(self, operation: str, required_phase: int) -> None:
        self.operation = operation
        self.required_phase = required_phase
        super().__init__(
            f"{operation} is not implemented in this build; "
            f"it is scheduled for Phase {required_phase}"
        )


class ParseError(NibblerError):
    """A structured source, decompression, or CIF syntax failure."""

    def __init__(
        self,
        *,
        code: str,
        message: str,
        source_name: str | None = None,
        line: int | None = None,
        column: int | None = None,
        byte_start: int | None = None,
        byte_end: int | None = None,
    ) -> None:
        self.code = code
        self.message = message
        self.source_name = source_name
        self.line = line
        self.column = column
        self.byte_start = byte_start
        self.byte_end = byte_end
        location = source_name or "<unknown>"
        if line is not None and column is not None:
            location = f"{location}:{line}:{column}"
        super().__init__(f"{location}: {code}: {message}")


class ProjectionError(ParseError):
    """A malformed projection plan or incompatible category occurrence."""


class SchemaError(NibblerError):
    """A schema selector, lock, artifact, or dictionary setup failure."""

    def __init__(self, *, code: str, message: str) -> None:
        self.code = code
        self.message = message
        super().__init__(f"{code}: {message}")


class ChemistryError(NibblerError):
    """A structured failure while constructing a semantic coordinate model."""

    def __init__(
        self, *, code: str, message: str, context: tuple[str, ...] = ()
    ) -> None:
        self.code = code
        self.message = message
        self.context = context
        super().__init__(f"{code}: {message}")


class WriteError(NibblerError):
    """A serialization or destination failure that leaves paths unmodified."""

    def __init__(
        self,
        *,
        code: str,
        message: str,
        destination: str | None = None,
    ) -> None:
        self.code = code
        self.message = message
        self.destination = destination
        location = "" if destination is None else f"{destination}: "
        super().__init__(f"{location}{code}: {message}")


class BatchError(ParseError):
    """One source failure collected by a deterministic batch scan."""

    def __init__(
        self,
        source_index: int,
        *,
        code: str,
        message: str,
        source_name: str | None = None,
        line: int | None = None,
        column: int | None = None,
        byte_start: int | None = None,
        byte_end: int | None = None,
    ) -> None:
        self.source_index = source_index
        super().__init__(
            code=code,
            message=message,
            source_name=source_name,
            line=line,
            column=column,
            byte_start=byte_start,
            byte_end=byte_end,
        )


class ValidationError(NibblerError):
    """Raised when a validation report contains error diagnostics."""

    diagnostics: tuple[Diagnostic, ...]

    def __init__(self, diagnostics: tuple[Diagnostic, ...]) -> None:
        self.diagnostics = diagnostics
        super().__init__(f"validation failed with {len(diagnostics)} error(s)")
