"""Source rows keep every retained target strict while preserving context."""

from pathlib import Path
import json

import pytest

import nibbler


ADDUCT = Path(__file__).parent / "reduce3" / "linked_9hj3.cif"


def source():
    gemmi = pytest.importorskip("gemmi")
    document = gemmi.cif.read_file(str(ADDUCT))
    document.sole_block().find_mmcif_category("_struct_conn.").erase()
    return nibbler.cif.read(document.as_string().encode())


def keys(document, component):
    pa = pytest.importorskip("pyarrow")
    table = pa.RecordBatchReader.from_stream(document.project("atom_site")).read_all().to_pylist()
    result = set()
    for atom in table:
        if atom["label_comp_id"] != component or atom["type_symbol"] in {"H", "D", "T"}:
            continue
        protein = atom["group_PDB"] == "ATOM"
        result.add((atom["label_asym_id"], str(atom["label_seq_id"] if protein else atom["auth_seq_id"]),
                    atom["pdbx_PDB_ins_code"] or "?", component, "protein" if protein else atom["label_atom_id"]))
    return json.dumps(sorted(result))


def run(document, **options):
    try:
        return nibbler.reduce.run(document, **options)
    except nibbler.ChemistryError as error:
        if "REDUCE_MONOMER_LIBRARY_NOT_FOUND" in str(error):
            pytest.skip("chem_data not found (set REDUCE3_CHEM_DATA)")
        raise


def test_context_parent_rows_keep_heavy_atoms_and_existing_output():
    document = source()
    actual = run(document, stop_on_any_missing_hydrogen=True,
                 model_id=1, required_hydrogen_target_keys=keys(document, "CYS"))
    reference = run(document)
    assert actual.document.to_canonical() == reference.document.to_canonical()
    assert "Removed unplaced context hydrogens" in actual.report
    assert "H351" in actual.report and "H352" in actual.report
    with pytest.raises(nibbler.ChemistryError, match="H351"):
        run(document, stop_on_any_missing_hydrogen=True,
            model_id=1, required_hydrogen_target_keys=keys(document, "Z41"))


@pytest.mark.parametrize("invalid", [b"[]", [], True, 0, 1.0, {}])
def test_target_keys_require_a_json_string(invalid):
    with pytest.raises(TypeError, match="JSON string"):
        nibbler.reduce.run(source(), required_hydrogen_target_keys=invalid)


@pytest.mark.parametrize("invalid", ["[]", "null", "{}", "not JSON", '[[]]', '[["X","9","?","CYS","protein"]]'])
def test_target_keys_reject_empty_invalid_or_absent_keys(invalid):
    with pytest.raises(nibbler.ChemistryError):
        run(source(), model_id=1, stop_on_any_missing_hydrogen=True, required_hydrogen_target_keys=invalid)


def test_duplicate_target_keys_fail():
    document = source()
    wanted = json.loads(keys(document, "CYS"))
    with pytest.raises(nibbler.ChemistryError, match="duplicates"):
        run(document, model_id=1, stop_on_any_missing_hydrogen=True,
            required_hydrogen_target_keys=json.dumps(wanted + wanted))


@pytest.mark.parametrize("options", [{"compat": True}, {"approach": "optimize"}, {"stop_on_any_missing_hydrogen": False}])
def test_target_keys_require_strict_fixed_mode_addition(options):
    document = source()
    with pytest.raises(nibbler.ChemistryError, match="strict fixed-mode"):
        run(document, **(dict(model_id=1, stop_on_any_missing_hydrogen=True,
            required_hydrogen_target_keys=keys(document, "CYS")) | options))


def test_target_keys_require_an_explicit_model():
    document = source()
    with pytest.raises(nibbler.ChemistryError, match="explicit model_id"):
        run(document, stop_on_any_missing_hydrogen=True, required_hydrogen_target_keys=keys(document, "CYS"))


def test_unpublished_source_row_option_is_not_a_compatibility_path():
    with pytest.raises(TypeError, match="unexpected keyword"):
        nibbler.reduce.run(source(), required_hydrogen_parent_sources=[0])
