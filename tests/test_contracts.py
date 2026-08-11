"""Tests for immutable public contract types."""

from __future__ import annotations

import pytest

from nibbler import (
    Diagnostic,
    MissingKind,
    Severity,
    SourceSpan,
    ValidationError,
    ValidationReport,
)


def test_missing_states_remain_distinct() -> None:
    assert MissingKind.UNKNOWN.value == "?"
    assert MissingKind.NOT_APPLICABLE.value == "."
    assert MissingKind.UNKNOWN is not MissingKind.NOT_APPLICABLE


@pytest.mark.parametrize(
    ("arguments", "message"),
    [
        ({"byte_start": -1, "byte_end": 0}, "byte_start"),
        ({"byte_start": 2, "byte_end": 1}, "byte_end"),
        ({"byte_start": 0, "byte_end": 1, "line": 0}, "line"),
        ({"byte_start": 0, "byte_end": 1, "column": 1}, "column requires line"),
    ],
)
def test_source_span_rejects_invalid_locations(
    arguments: dict[str, int], message: str
) -> None:
    with pytest.raises(ValueError, match=message):
        SourceSpan(**arguments)


def test_validation_report_raises_only_for_errors() -> None:
    warning = Diagnostic("CIF0001", Severity.WARNING, "accepted deviation")
    error = Diagnostic(
        "CIF_LOOP_VALUE_COUNT",
        Severity.ERROR,
        "loop value count is not divisible by tag count",
        SourceSpan(10, 11, line=4, column=1),
    )

    warning_report = ValidationReport((warning,))
    warning_report.raise_for_errors()
    assert warning_report.is_valid

    error_report = ValidationReport.from_iterable([warning, error])
    assert not error_report.is_valid
    assert error_report.errors == (error,)
    with pytest.raises(ValidationError) as raised:
        error_report.raise_for_errors()
    assert raised.value.diagnostics == (error,)
