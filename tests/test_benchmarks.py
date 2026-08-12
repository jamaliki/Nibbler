"""Tests for safeguards in the executable benchmark harness."""

from __future__ import annotations

import pytest

from benchmarks.pdb_stress import selected_formats
from benchmarks.run import require_release_nibbler


def test_benchmark_accepts_release_extension(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr("nibbler._core.build_profile", lambda: "release")

    require_release_nibbler()


def test_benchmark_rejects_debug_extension(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def debug_profile() -> str:
        return "debug"

    monkeypatch.setattr("nibbler._core.build_profile", debug_profile)

    with pytest.raises(RuntimeError, match="make develop-release"):
        require_release_nibbler()


def test_pdb_stress_format_selection_includes_pinned_gzip() -> None:
    assert selected_formats("gzip") == ("cif.gz",)
    assert selected_formats("all") == ("cif", "bcif", "cif.gz")
