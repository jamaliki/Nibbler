"""Small immutable types shared by the public APIs."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from enum import Enum, unique

from .errors import ValidationError


@unique
class MissingKind(str, Enum):
    """The two missing-value states defined by CIF."""

    UNKNOWN = "?"
    NOT_APPLICABLE = "."


@unique
class Severity(str, Enum):
    """Severity of a structured diagnostic."""

    INFO = "info"
    WARNING = "warning"
    ERROR = "error"


@unique
class Profile(str, Enum):
    """Semantic validation and writing profiles."""

    PDBX = "pdbx"
    MODELCIF = "modelcif"


@dataclass(frozen=True, slots=True)
class SourceSpan:
    """Half-open byte span with an optional one-based display location."""

    byte_start: int
    byte_end: int
    line: int | None = None
    column: int | None = None

    def __post_init__(self) -> None:
        if self.byte_start < 0:
            raise ValueError("byte_start must be non-negative")
        if self.byte_end < self.byte_start:
            raise ValueError("byte_end must not precede byte_start")
        if self.line is not None and self.line < 1:
            raise ValueError("line must be one-based")
        if self.column is not None and self.column < 1:
            raise ValueError("column must be one-based")
        if self.column is not None and self.line is None:
            raise ValueError("column requires line")


@dataclass(frozen=True, slots=True)
class Diagnostic:
    """Stable machine-readable diagnostic with optional source context."""

    code: str
    severity: Severity
    message: str
    span: SourceSpan | None = None
    context: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not self.code:
            raise ValueError("diagnostic code must not be empty")
        if not self.message:
            raise ValueError("diagnostic message must not be empty")


@dataclass(frozen=True, slots=True)
class ValidationReport:
    """Immutable collection of validation diagnostics."""

    diagnostics: tuple[Diagnostic, ...] = ()

    @classmethod
    def from_iterable(cls, diagnostics: Iterable[Diagnostic]) -> ValidationReport:
        """Construct a report while taking an immutable snapshot."""
        return cls(tuple(diagnostics))

    @property
    def errors(self) -> tuple[Diagnostic, ...]:
        """Return only error diagnostics."""
        return tuple(
            diagnostic
            for diagnostic in self.diagnostics
            if diagnostic.severity is Severity.ERROR
        )

    @property
    def is_valid(self) -> bool:
        """Return whether the report contains no error diagnostics."""
        return not self.errors

    def raise_for_errors(self) -> None:
        """Raise :class:`ValidationError` if errors are present."""
        errors = self.errors
        if errors:
            raise ValidationError(errors)
