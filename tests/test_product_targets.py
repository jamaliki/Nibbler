"""Occurrence-level product evidence must precede target selection."""

import json

import pytest

from test_observation_ligand_hydrogens import arrays, document, ligand_document, ligand_packer


def target(hydrogens, charge=-1, complete=True):
    return dict(parent_key=["G", "101", "?", "LIG", "C4"], sequence_id=None,
                formal_charge=charge, hybridization=4, hydrogens=hydrogens,
                complete=complete, modified=True)


@pytest.mark.parametrize("prepared", [False, True])
@pytest.mark.parametrize("count,known", [(3, True), (2, False), (None, False)])
def test_product_targets_control_hydrogen_and_charge_labels(prepared, count, known):
    native = ligand_packer()
    doc = ligand_document(3)
    result = (native.pack_prepared(native.prepare_source(doc), complete=True,
                                  product_targets=json.dumps([target(count)])) if prepared else
              native.pack(doc, complete=True, product_targets=json.dumps([target(count)])))
    values, _ = arrays(result)
    assert values["ligand_charge"].tolist() == [-1, 0]
    assert values["ligand_hydrogen_known"].tolist() == [known, False]
    assert values["ligand_charge_known"].tolist() == [known, False]
    assert values["mask"][0, 15:].sum() == 3


def test_empty_product_targets_never_attest_unchecked_ligands():
    values, _ = arrays(ligand_packer().pack(ligand_document(3), complete=True, product_targets="[]"))
    assert not values["ligand_charge_known"].any()
    assert not values["ligand_hydrogen_known"].any()


def endpoint_document(author="2"):
    rows = ("HETATM C1 C B LIG 2 . 1 . 0 0 0 1\n"
            "HETATM C1 C B LIG 2 . 2 . 5 0 0 1\n"
            "HETATM C1 C C LIG 3 . 3 . 6.4 0 0 1\n")
    links = """loop_
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr1_auth_seq_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.ptnr2_auth_seq_id
""" + f"covale B LIG . C1 {author} C LIG . C1 3\n"
    return document(rows, links)


def test_native_endpoint_uses_author_sequence():
    native = ligand_packer()
    native.prepare_source(endpoint_document())
    with pytest.raises(ValueError, match="ambiguous atom identity"):
        native.prepare_source(endpoint_document("?"))
