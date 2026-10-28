"""Retained-target policy reuse keeps source checks and error order."""

import json

import gemmi
import nibbler
import pytest

from test_reduce_hydrogen_parents import document as parent_document, row
from test_reduce_source_model_ids import document, heavy_coordinates, target_keys
from test_reduce_target_masks import bridge_source


FLAGS = [(False, False), (False, True), (True, False), (True, True)]


def run(source, keys, flags, **options):
    return nibbler.reduce.run(
        source,
        model_id=options.pop("model_id", 1),
        required_hydrogen_target_keys=keys,
        stop_on_any_missing_hydrogen=True,
        allow_undefined_ligand_targets=flags[0],
        allow_unresolved_bond_targets=flags[1],
        **options,
    )


@pytest.mark.parametrize("flags", FLAGS)
@pytest.mark.parametrize(
    "keys,message",
    [
        ("[", "EOF"),
        (json.dumps([["A", "1", "?", "KHP", "C1"]] * 2), "duplicates"),
        ("[]", "Invalid coordinate number"),
    ],
)
def test_key_checks_precede_source_and_library_errors(flags, keys, message, tmp_path):
    source = parent_document([row("C1", "C", x="wrong")])
    with pytest.raises(nibbler.ChemistryError, match=message) as raised:
        run(source, keys, flags, chem_data=tmp_path)
    assert raised.value.code == "REDUCE_INVALID_MODEL"


@pytest.mark.parametrize("flags", FLAGS)
@pytest.mark.parametrize("models", [[2], [1, 2]])
def test_model_two_keeps_original_heavy_coordinates(flags, models):
    source = document(models)
    before = source.to_canonical()
    result = run(source, target_keys(source, 2), flags, model_id=2)
    assert heavy_coordinates(result.document, 2) == heavy_coordinates(source, 2)
    assert source.to_canonical() == before


@pytest.mark.parametrize("flags", FLAGS)
@pytest.mark.parametrize(
    "models,selected,message",
    [
        ([2], 1, "not present"),
        ([2, 3, 2], 2, "ambiguous"),
        ([2, "02"], 2, "ambiguous"),
    ],
)
def test_missing_or_ambiguous_models_keep_error_class(flags, models, selected, message):
    source = document(models)
    with pytest.raises(nibbler.ChemistryError, match=message) as raised:
        run(source, "[]", flags, model_id=selected)
    assert raised.value.code == "REDUCE_INVALID_MODEL"


def alternate_document(labels):
    source = gemmi.cif.read_string(document([2]).to_canonical())
    table = source.sole_block().find_mmcif_category("_atom_site.")
    column = list(table.tags).index("_atom_site.label_alt_id")
    for index, atom in enumerate(table):
        atom[column] = labels[index % len(labels)]
    return nibbler.cif.read(source.as_string().encode())


@pytest.mark.parametrize("flags", FLAGS)
def test_selected_alternate_keeps_source_coordinates(flags):
    source = alternate_document(["A"])
    result = run(source, target_keys(source, 2), flags, model_id=2, alt_id="A")
    assert heavy_coordinates(result.document, 2) == heavy_coordinates(source, 2)


@pytest.mark.parametrize("flags", FLAGS)
def test_unresolved_alternates_fail_before_library_loading(flags, tmp_path):
    source = alternate_document(["A", "B"])
    with pytest.raises(
        nibbler.ChemistryError, match="Resolve alternate conformers"
    ) as raised:
        run(source, target_keys(source, 2), flags, model_id=2, chem_data=tmp_path)
    assert raised.value.code == "REDUCE_INVALID_MODEL"


def test_late_bridge_library_error_keeps_plain_value_error(tmp_path):
    keys = json.dumps([["L", "601", "?", "EJ2", "FE"]])
    with pytest.raises(ValueError, match="REDUCE_MONOMER_LIBRARY_NOT_FOUND") as raised:
        run(bridge_source(), keys, (True, False), chem_data=tmp_path)
    assert not isinstance(raised.value, nibbler.ChemistryError)


def test_required_bridge_library_error_keeps_chemistry_error(tmp_path):
    keys = json.dumps([["L", "601", "?", "EJ2", "H"]])
    with pytest.raises(
        nibbler.ChemistryError, match="REDUCE_MONOMER_LIBRARY_NOT_FOUND"
    ) as raised:
        run(bridge_source(), keys, (True, False), chem_data=tmp_path)
    assert raised.value.code == "REDUCE_INVALID_MODEL"


@pytest.mark.parametrize("flags", FLAGS)
def test_ligand_library_error_keeps_reduce_code(flags, tmp_path):
    source = parent_document([row("C1", "C")])
    keys = json.dumps([["A", "1", "?", "KHP", "C1"]])
    with pytest.raises(nibbler.ChemistryError) as raised:
        run(source, keys, flags, chem_data=tmp_path)
    assert raised.value.code == "REDUCE_MONOMER_LIBRARY_NOT_FOUND"
