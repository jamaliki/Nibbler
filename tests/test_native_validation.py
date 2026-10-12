"""End-to-end native validation checks on parsed documents."""

import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pytest

from nibbler import cif
from nibbler.geometry import validate

CASES = Path(__file__).with_name("native_validation_reference_cases")


@pytest.fixture
def chemistry():
    configured = os.environ.get("NIBBLER_VALIDATION_CHEM_DATA")
    if not configured:
        pytest.skip("set NIBBLER_VALIDATION_CHEM_DATA to the pinned chemistry tree")
    return Path(configured)


@pytest.mark.parametrize("name", ["valid_alanine", "valid_d_alanine"])
def test_native_graph_covers_valid_controls(chemistry, name):
    document = cif.read(CASES / f"{name}.cif")
    result = validate(document, chem_data=chemistry, use_cdl=False, include_graph=True)
    assert result["complete"]
    assert result["failures"] == []
    assert result["checks"]["bonds"]["count"] == 9
    assert result["checks"]["angles"]["count"] == 14
    assert result["checks"]["chiralities"]["count"] == 1
    assert result["contacts"]["clashscore"] == 0
    assert len(result["atoms"]) == 10


def test_native_hydrogen_outlier_keeps_input_coordinates(chemistry):
    result = validate(CASES / "displaced_hydrogen.cif", chem_data=chemistry,
                      use_cdl=False, include_graph=True)
    assert result["complete"]
    assert any(row["check"] == "bonds" for row in result["failures"])
    assert any(row["name"] == "HA" for row in result["atoms"])


def test_cached_chemistry_is_safe_across_threads(chemistry):
    document = cif.read(CASES / "displaced_hydrogen.cif")
    expected = validate(document, chem_data=chemistry, use_cdl=False)

    def check(_):
        return validate(document, chem_data=chemistry, use_cdl=False)

    with ThreadPoolExecutor(max_workers=2) as executor:
        assert list(executor.map(check, range(4))) == [expected] * 4


def test_default_findings_keep_source_atom_ids_without_a_full_graph(chemistry):
    result = validate(CASES / "displaced_hydrogen.cif", chem_data=chemistry,
                      use_cdl=False)
    assert "atoms" not in result
    assert "graph" not in result
    assert result["failures"]
    assert all(row["atom_ids"] and row["source_rows"] for row in result["failures"])


def test_multiple_coordinate_blocks_cannot_be_omitted(chemistry):
    source = b"data_a\nloop_\n_atom_site.id\n1\ndata_b\nloop_\n_atom_site.id\n2\n"
    with pytest.raises(ValueError, match="Multiple coordinate blocks"):
        validate(source, chem_data=chemistry)


def test_invalid_formal_charge_cannot_become_neutral(chemistry):
    source = (b"data_charge\nloop_\n_atom_site.id\n"
              b"_atom_site.pdbx_formal_charge\n1 bogus\n")
    with pytest.raises(ValueError, match="invalid formal charge"):
        validate(source, chem_data=chemistry)


@pytest.mark.parametrize("value", [0, -1, float("inf"), float("nan")])
def test_native_invalid_cutoff_fails(chemistry, value):
    with pytest.raises(ValueError, match="sigma_cutoff"):
        validate(CASES / "valid_alanine.cif", chem_data=chemistry, sigma_cutoff=value)


@pytest.mark.parametrize("flag", ["nuclear", "use_cdl", "include_graph"])
def test_native_flags_reject_integers(flag):
    with pytest.raises(TypeError, match=flag):
        validate(b"data_empty\n", chem_data="unused", **{flag: 1})
