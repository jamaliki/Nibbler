"""Use the producer's dictionary and original source rows for H targets."""

import json
from pathlib import Path

import nibbler
import pytest


HEADER = ("data_parents\nloop_\n_atom_site.group_PDB\n_atom_site.label_entity_id\n"
          "_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.auth_seq_id\n"
          "_atom_site.label_comp_id\n_atom_site.label_atom_id\n_atom_site.type_symbol\n"
          "_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n"
          "_atom_site.Cartn_z\n_atom_site.occupancy\n_atom_site.pdbx_PDB_model_num\n")


def document(rows, protein=False):
    prefix = "data_parents\n_entity_poly.entity_id 1\n_entity_poly.type 'polypeptide(L)'\n" if protein else ""
    text = (HEADER if not prefix else prefix + HEADER.removeprefix("data_parents\n")) + "".join(rows)
    return nibbler.cif.read(text.encode())


def row(atom, element, *, component="KHP", chain="A", seq=".", auth="1", model=1, group="HETATM", entity=2, x=0, occupancy=1):
    return f"{group} {entity} {chain} {seq} {auth} {component} {atom} {element} ? {x} 0 0 {occupancy} {model}\n"


@pytest.fixture
def library(tmp_path):
    for folder in ("geostd/list", "mon_lib/list", "geostd/k", "chemical_components/k"):
        (tmp_path / folder).mkdir(parents=True)
    (tmp_path / "mon_lib/list/mon_lib_list.cif").write_text("data_empty\n")
    (tmp_path / "geostd/ener_lib.cif").write_text("data_empty\n")
    geo = ("data_comp_KHP\nloop_\n_chem_comp_atom.atom_id\n_chem_comp_atom.type_symbol\n"
           "_chem_comp_atom.type_energy\nC1 C CH1\nH1 H HCH1\nloop_\n"
           "_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n_chem_comp_bond.type\n"
           "_chem_comp_bond.value_dist\nC1 H1 single 1.09\n")
    (tmp_path / "geostd/k/data_KHP.cif").write_text(geo)
    ccd = ("data_KHP\n_chem_comp.id KHP\n_chem_comp.type NON-POLYMER\nloop_\n"
           "_chem_comp_atom.atom_id\n_chem_comp_atom.type_symbol\n_chem_comp_atom.charge\n"
           "_chem_comp_atom.pdbx_model_Cartn_x_ideal\n_chem_comp_atom.pdbx_model_Cartn_y_ideal\n"
           "_chem_comp_atom.pdbx_model_Cartn_z_ideal\nC2 C 0 0 0 0\nH1 H 0 1.09 0 0\nloop_\n"
           "_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n_chem_comp_bond.value_order\nC2 H1 SING\n")
    (tmp_path / "chemical_components/k/data_KHP.cif").write_text(ccd)
    return tmp_path


def test_dictionary_selection_matches_geostd_and_ccd_fallback(library):
    doc = document([row("C1", "C", model=1), row("C2", "C", model=2)])
    assert nibbler.reduce.hydrogen_parents(doc, chem_data=library, model_id=1) == {"KHP": {"H1": "C1"}}
    assert nibbler.reduce.hydrogen_parents(doc, chem_data=library, model_id=2) == {"KHP": {"H1": "C2"}}


def test_dictionary_cache_keeps_checked_parents_without_reopening_source(library):
    doc = document([row("C1", "C")])
    expected = nibbler.reduce.hydrogen_parents(doc, chem_data=library)
    (library / "geostd/k/data_KHP.cif").unlink()
    assert nibbler.reduce.hydrogen_parents(doc, chem_data=library) == expected


def test_same_component_with_conflicting_selected_dictionaries_fails(library):
    doc = document([row("C1", "C"), row("C2", "C", chain="B")])
    with pytest.raises(ValueError, match="conflicting hydrogen dictionaries"):
        nibbler.reduce.hydrogen_parents(doc, chem_data=library)


def test_unavailable_source_dictionary_stays_distinct_from_no_hydrogens(library):
    doc = document([row("C1", "C", component="UNL")])
    assert nibbler.reduce.hydrogen_parents(doc, chem_data=library) == {}
    (library / "geostd/k/data_KHP.cif").write_text("data_comp_KHP\nloop_\n"
        "_chem_comp_atom.atom_id\n_chem_comp_atom.type_symbol\nC1 C\n")
    assert nibbler.reduce.hydrogen_parents(document([row("C1", "C")]), chem_data=library) == {"KHP": {}}


@pytest.mark.parametrize("model", [0, -1, True, 1.5])
def test_wrapper_rejects_bad_models(model):
    with pytest.raises(ValueError, match="one-based"):
        nibbler.reduce.hydrogen_parents(document([row("C1", "C")]), model_id=model)


def test_wrapper_requires_parsed_document():
    with pytest.raises(TypeError, match="parsed CIF document"):
        nibbler.reduce.hydrogen_parents(Path("unused.cif"))


def sources(doc, keys, model=1):
    return nibbler.reduce.required_parent_sources(doc, json.dumps(keys), model_id=model)


@pytest.mark.parametrize("model", [0, -1, True, 1.5])
def test_parent_source_wrapper_rejects_bad_models(model):
    with pytest.raises(ValueError, match="one-based"):
        nibbler.reduce.required_parent_sources(document([row("C1", "C")]), "[]", model_id=model)


def test_parent_source_wrapper_requires_document_and_json_string():
    with pytest.raises(TypeError, match="parsed CIF document"):
        nibbler.reduce.required_parent_sources(Path("unused.cif"), "[]")
    with pytest.raises(TypeError, match="JSON string"):
        nibbler.reduce.required_parent_sources(document([row("C1", "C")]), [])


def test_required_parent_rows_keep_original_ordinals_and_selected_model():
    doc = document([row("C1", "C", model=1), row("H1", "H", model=2),
                    row("C1", "C", model=2), row("C2", "C", model=2)])
    key = ["A", "1", "?", "KHP", "C2"]
    assert sources(doc, [key], model=2) == [3]
    assert sources(doc, [], model=2) == []


def test_required_parent_rows_resolve_selected_duplicate_occupancy():
    doc = document([row("C1", "C", occupancy=.2), row("C1", "C", occupancy=.8)])
    assert sources(doc, [["A", "1", "?", "KHP", "C1"]]) == [1]


def test_required_parent_rows_exclude_context_but_keep_entire_retained_residue():
    doc = document([row("N", "N", group="ATOM", entity=1, component="ALA", seq="2"),
                    row("CA", "C", group="ATOM", entity=1, component="ALA", seq="2"),
                    row("HB1", "H", group="ATOM", entity=1, component="ALA", seq="2"),
                    row("N", "N", group="ATOM", entity=1, component="UNK", seq="3", auth="3"),
                    row("CB", "C", group="ATOM", entity=1, component="UNK", seq="3", auth="3")], protein=True)
    assert sources(doc, [["A", "2", "?", "ALA", "protein"]]) == [0, 1]


@pytest.mark.parametrize("keys", [[['A', '1', '?', 'KHP', 'C3']],
                                  [['A', '1', '?', 'KHP', 'C1']] * 2])
def test_required_parent_keys_reject_missing_or_duplicate_targets(keys):
    with pytest.raises(ValueError, match="Retained hydrogen-parent keys"):
        sources(document([row("C1", "C")]), keys)


def test_deposited_hydrogen_bond_precedes_dictionary_parent():
    import numpy as np
    from test_observation_protein_parents import packer

    native, _ = packer()
    native.add_hydrogen_parents("KHP", {"H1": "C1"})
    rows = [row("C1", "C"), row("C2", "C", x=2), row("H1", "H", x=2.97)]
    with pytest.raises(ValueError, match="Invalid ligand hydrogen parent"):
        native.pack(document(rows), complete=True, validate=True)
    connection = "#\nloop_\n" + "".join(f"_struct_conn.{name}\n" for name in (
        "conn_type_id", "ptnr1_label_asym_id", "ptnr1_label_seq_id", "ptnr1_label_comp_id",
        "ptnr1_label_atom_id", "ptnr2_label_asym_id", "ptnr2_label_seq_id", "ptnr2_label_comp_id",
        "ptnr2_label_atom_id", "ptnr1_symmetry", "ptnr2_symmetry"))
    connection += "covale A . KHP H1 A . KHP C2 1_555 1_555\n"
    buffers, _ = native.pack(document([*rows, connection]), complete=True, validate=True)
    masks = np.frombuffer(buffers["mask"][2], dtype=bool).reshape(buffers["mask"][1])
    assert not masks[0, 15] and masks[1, 15]
