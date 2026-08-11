"""Independent syntax checks against committed reference engines."""

from __future__ import annotations

from pathlib import Path

import gemmi
import pytest

FIXTURE_ROOT = Path(__file__).parent / "fixtures"


@pytest.mark.parametrize(
    "relative_file",
    [
        "syntax/minimal.cif",
        "syntax/missing_values.cif",
        "chemistry/ligand_ion_water.cif",
    ],
)
def test_gemmi_accepts_valid_syntax_fixtures(relative_file: str) -> None:
    document = gemmi.cif.read_file(str(FIXTURE_ROOT / relative_file))
    assert len(document) == 1


def test_gemmi_rejects_incomplete_loop_row() -> None:
    with pytest.raises(ValueError, match="Wrong number of values in loop"):
        gemmi.cif.read_file(str(FIXTURE_ROOT / "syntax/malformed_loop.cif"))
