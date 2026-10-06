"""Public Python coverage for running Reduce3 on documents in memory."""

from __future__ import annotations

import os
from pathlib import Path

import pytest

import nibbler
from nibbler import ChemistryError, CifDocument
from nibbler.reduce import Reduction

CRAMBIN = Path(__file__).parent / "reduce3" / "crambin_1-5.cif"

pytestmark = pytest.mark.skipif(
    not nibbler.reduce.available(), reason="built without the reduce3 feature"
)


def _chem_data_or_skip() -> None:
    try:
        nibbler.reduce.run(CRAMBIN, approach="remove")
    except ChemistryError as error:
        if error.code == "REDUCE_MONOMER_LIBRARY_NOT_FOUND":
            pytest.skip("chem_data not found (set REDUCE3_CHEM_DATA)")
        raise


def _atom_count(document: CifDocument) -> int:
    table = nibbler.cif.read(document.to_canonical().encode(), category="atom_site")
    return len(table)


def test_hydrogens_are_added_to_a_parsed_document() -> None:
    _chem_data_or_skip()
    document = nibbler.chomp(CRAMBIN)

    result = nibbler.reduce.run(document, add_flip_movers=True)

    assert isinstance(result, Reduction)
    assert result.document.block_count == 1
    assert result.document.to_canonical().startswith("data_crn5\n")
    assert _atom_count(document) == 33
    assert _atom_count(result.document) == 66
    assert "MoverSingleHydrogenRotator" in result.report


def test_paths_documents_and_compat_mode_agree() -> None:
    _chem_data_or_skip()
    from_path = nibbler.reduce.run(CRAMBIN, compat=True)
    from_document = nibbler.reduce.run(nibbler.chomp(CRAMBIN), compat=True)

    assert from_path.document.to_canonical() == from_document.document.to_canonical()


def test_remove_strips_the_added_hydrogens() -> None:
    _chem_data_or_skip()
    added = nibbler.reduce.run(CRAMBIN)

    removed = nibbler.reduce.run(added.document, approach="remove")

    assert _atom_count(removed.document) == 33


def test_named_library_directory_must_exist(tmp_path: Path) -> None:
    with pytest.raises(ChemistryError) as raised:
        nibbler.reduce.run(CRAMBIN, chem_data=tmp_path)
    assert raised.value.code == "REDUCE_MONOMER_LIBRARY_NOT_FOUND"


def test_documents_without_a_model_are_rejected() -> None:
    _chem_data_or_skip()
    with pytest.raises(ChemistryError) as raised:
        nibbler.reduce.run(b"data_empty\n_entry.id empty\n")
    assert raised.value.code == "REDUCE_INVALID_MODEL"


@pytest.mark.parametrize(
    ("options", "error"),
    [
        ({"approach": "replace"}, ValueError),
        ({"n_terminal_charge": "always"}, ValueError),
        ({"model_id": 0}, ValueError),
        ({"bonded_neighbor_depth": -1}, ValueError),
        ({"add_flip_movers": 1}, TypeError),
        ({"preference_magnitude": "1"}, TypeError),
        ({"probe": {"probe_radius": True}}, TypeError),
        ({"probe": {"set_polar_hydrogen_radius": 1}}, TypeError),
        ({"probe": {"radius": 0.3}}, ValueError),
    ],
)
def test_options_are_validated(
    options: dict[str, object], error: type[Exception]
) -> None:
    with pytest.raises(error):
        nibbler.reduce.run(CRAMBIN, **options)  # type: ignore[arg-type]


def test_explicit_chem_data_directory(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    chem_data = os.environ.get("REDUCE3_CHEM_DATA")
    if chem_data is None:
        pytest.skip("REDUCE3_CHEM_DATA is not set")
    monkeypatch.delenv("REDUCE3_CHEM_DATA")

    result = nibbler.reduce.run(CRAMBIN, chem_data=Path(chem_data))

    assert _atom_count(result.document) == 66
