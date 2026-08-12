"""Tests for safeguards in the executable benchmark harness."""

from __future__ import annotations

import subprocess
import sys

import pytest

from benchmarks.check_pdb_report import Guardrails, check
from benchmarks.pdb_stress import selected_formats
from benchmarks.run import require_release_nibbler

GUARDRAILS = Guardrails(
    minimum_logical_input_bytes=8_000_000,
    maximum_peak_rss_bytes=4_000_000_000,
    peak_rss_fixed_allowance_bytes=160_000_000,
    maximum_incremental_rss_to_logical_input=8.0,
    minimum_projection={"cif": 200.0},
    minimum_full_document={"cif": 150.0},
)


def pdb_report(
    *,
    projections_equal: bool = True,
    logical_input_bytes: int = 10_000_000,
    projection_rate: float = 300.0,
    full_rate: float = 250.0,
    peak_rss: int = 200_000_000,
) -> dict[str, object]:
    """Return one minimal large-input PDB benchmark report."""
    return {
        "structures": [
            {
                "pdb_id": "test",
                "projections_equal": projections_equal,
                "formats": [
                    {
                        "format": "cif",
                        "logical_input_bytes": logical_input_bytes,
                        "native_project": {
                            "input_mb_per_second": projection_rate,
                        },
                        "full_document": {
                            "stage": {"input_mb_per_second": full_rate},
                            "peak_rss_bytes": peak_rss,
                        },
                    }
                ],
            }
        ]
    }


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


def test_pdb_guardrails_accept_healthy_large_input() -> None:
    assert check(pdb_report(), GUARDRAILS) == ()


def test_pdb_guardrails_report_independent_failures() -> None:
    failures = check(
        pdb_report(
            projections_equal=False,
            projection_rate=199.0,
            full_rate=149.0,
            peak_rss=241_000_000,
        ),
        GUARDRAILS,
    )

    assert len(failures) == 4
    assert any("projections differ" in failure for failure in failures)
    assert any("projection 199.0 MB/s" in failure for failure in failures)
    assert any("full document 149.0 MB/s" in failure for failure in failures)
    assert any("peak RSS" in failure for failure in failures)


def test_pdb_guardrails_ignore_noisy_small_inputs() -> None:
    assert (
        check(
            pdb_report(
                logical_input_bytes=1_000_000,
                projection_rate=0.0,
                full_rate=0.0,
                peak_rss=4_100_000_000,
            ),
            GUARDRAILS,
        )
        == ()
    )


@pytest.mark.parametrize("module", ("tools.fetch_pdb_corpus", "tools.pgo_train"))
def test_corpus_tools_support_module_entry_points(module: str) -> None:
    subprocess.run(
        [sys.executable, "-m", module, "--help"],
        check=True,
        capture_output=True,
        text=True,
    )
