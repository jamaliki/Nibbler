"""Keep excluded heavy contact context without weakening strict target checks."""

import os
from pathlib import Path

import numpy as np
import pytest

import nibbler
from token_tables import token_table

CRAMBIN = Path(__file__).parent / "reduce3/crambin_1-5.cif"
ADDUCT = Path(__file__).parent / "reduce3/linked_9hj3.cif"
EXCLUDED = ["D8U", "DOD", "HOH", "UNL"]


def document(components=("UNL",)):
    gemmi = pytest.importorskip("gemmi")
    structure = gemmi.read_structure(str(CRAMBIN))
    for index, component in enumerate(components):
        chain = gemmi.Chain(f"U{index}")
        residue = gemmi.Residue()
        residue.name, residue.seqid, residue.het_flag = (
            component,
            gemmi.SeqId(900 + index, " "),
            "H",
        )
        atom = gemmi.Atom()
        atom.name, atom.element, atom.pos = (
            "C1",
            gemmi.Element("C"),
            gemmi.Position(100 + index, 100, 100),
        )
        residue.add_atom(atom)
        chain.add_residue(residue)
        structure[0].add_chain(chain)
    return nibbler.cif.read(structure.make_mmcif_document().as_string().encode())


def run(source, **options):
    try:
        return nibbler.reduce.run(source, **options)
    except nibbler.ChemistryError as error:
        if error.code == "REDUCE_MONOMER_LIBRARY_NOT_FOUND":
            pytest.skip("chem_data not found (set REDUCE3_CHEM_DATA)")
        raise


def test_allowed_unl_keeps_heavy_atoms_and_matches_existing_unknown_mode():
    source = document()
    with pytest.raises(nibbler.ChemistryError, match="UNL"):
        run(source, stop_on_any_missing_hydrogen=True)
    allowed = run(
        source, stop_on_any_missing_hydrogen=True, allow_missing_restraints_for=EXCLUDED
    )
    reference = run(source)
    assert allowed.document.to_canonical() == reference.document.to_canonical()
    assert "Allowed missing restraints for components UNL" in allowed.report
    pa = pytest.importorskip("pyarrow")
    rows = (
        pa.RecordBatchReader.from_stream(allowed.document.project("atom_site"))
        .read_all()
        .to_pylist()
    )
    unl = [row for row in rows if row["label_comp_id"] == "UNL"]
    assert (
        len(unl) == 1
        and unl[0]["label_atom_id"] == "C1"
        and unl[0]["type_symbol"] == "C"
    )
    assert [float(unl[0][key]) for key in ("Cartn_x", "Cartn_y", "Cartn_z")] == [
        100.0,
        100.0,
        100.0,
    ]


def test_unlisted_missing_component_still_fails():
    source = document(("UNL", "NO_SUCH_COMPONENT"))
    with pytest.raises(nibbler.ChemistryError, match="NO_SUCH_COMPONENT") as caught:
        run(
            source,
            stop_on_any_missing_hydrogen=True,
            allow_missing_restraints_for=EXCLUDED,
        )
    assert "UNL" not in str(caught.value)


def test_known_component_outputs_are_unchanged_with_the_allowlist():
    strict = run(CRAMBIN, stop_on_any_missing_hydrogen=True)
    allowed = run(
        CRAMBIN,
        stop_on_any_missing_hydrogen=True,
        allow_missing_restraints_for=[*EXCLUDED, "THR", "CYS", "PRO"],
    )
    assert allowed.document.to_canonical() == strict.document.to_canonical()
    assert "Allowed missing restraints" not in allowed.report


def test_allowed_missing_dictionary_does_not_allow_unplaced_hydrogens():
    gemmi = pytest.importorskip("gemmi")
    source = gemmi.cif.read_file(str(ADDUCT))
    source.sole_block().find_mmcif_category("_struct_conn.").erase()
    with pytest.raises(nibbler.ChemistryError, match="H351") as caught:
        run(
            nibbler.cif.read(source.as_string().encode()),
            stop_on_any_missing_hydrogen=True,
            allow_missing_restraints_for=[*EXCLUDED, "Z41"],
        )
    assert "H352" in str(caught.value)


def test_known_native_labels_match_the_existing_lenient_control():
    schema = os.environ.get("KAVEH_ATOM32_SCHEMA")
    if not schema:
        pytest.skip("Set KAVEH_ATOM32_SCHEMA to the pinned model schema")
    from nibbler.observations import ObservationPacker

    source = document()
    allowed = run(
        source, stop_on_any_missing_hydrogen=True, allow_missing_restraints_for=EXCLUDED
    )
    control = run(source)
    packer = ObservationPacker(Path(schema).read_text(), token_table("C"), EXCLUDED)
    actual, keys = packer.pack(allowed.document, complete=True, validate=True)
    expected, expected_keys = packer.pack(
        control.document, complete=True, validate=True
    )
    assert keys == expected_keys and all(key[3] != "UNL" for key in keys)
    assert actual.keys() == expected.keys()
    for name in actual:
        dtype, shape, buffer = actual[name]
        old_dtype, old_shape, old_buffer = expected[name]
        assert dtype == old_dtype and shape == old_shape
        np.testing.assert_array_equal(
            np.frombuffer(buffer, dtype=dtype), np.frombuffer(old_buffer, dtype=dtype)
        )


@pytest.mark.parametrize("invalid", ["UNL", [1], {"UNL"}])
def test_python_rejects_non_sequence_component_ids(invalid):
    with pytest.raises(TypeError, match="sequence of component identifiers"):
        run(CRAMBIN, allow_missing_restraints_for=invalid)


@pytest.mark.parametrize("invalid", [[""], ["?"], ["."], ["UNL ZZZ"]])
def test_native_rejects_empty_or_joined_component_ids(invalid):
    with pytest.raises(nibbler.ChemistryError, match="component identifiers"):
        run(CRAMBIN, allow_missing_restraints_for=invalid)


def test_compat_mode_does_not_change_its_contract():
    with pytest.raises(nibbler.ChemistryError, match="require fixed mode"):
        run(CRAMBIN, compat=True, allow_missing_restraints_for=EXCLUDED)
