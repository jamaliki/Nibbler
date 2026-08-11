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
            f"{operation} is not implemented in Phase 0; "
            f"it is scheduled for Phase {required_phase}"
        )


class ValidationError(NibblerError):
    """Raised when a validation report contains error diagnostics."""

    diagnostics: tuple[Diagnostic, ...]

    def __init__(self, diagnostics: tuple[Diagnostic, ...]) -> None:
        self.diagnostics = diagnostics
        super().__init__(f"validation failed with {len(diagnostics)} error(s)")
