"""Hydrogen removal keeps every non-hydrogen row and category as authored."""

import pytest

import nibbler

SOURCE = (
    b"data_x\n_entry.id x\n"
    b"loop_\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n"
    b"1 N N\n2 H H\n3 C CA\n4 d D1\n5 ? X\n"
    b"loop_\n_atom_site_anisotrop.id\n_atom_site_anisotrop.U[1][1]\n"
    b"1 0.1\n2 0.2\n9 0.4\n"
    b"loop_\n_struct_conn.id\n_struct_conn.ptnr1_label_atom_id\n"
    b"_struct_conn.ptnr2_label_atom_id\nc1 H CA\n"
)


def test_removes_hydrogen_rows_and_their_tensors():
    document = nibbler.cif.read(SOURCE)
    result = nibbler.mmcif.without_hydrogens(document)
    expected = nibbler.cif.read(
        b"data_x\n_entry.id x\n"
        b"loop_\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n"
        b"1 N N\n3 C CA\n5 ? X\n"
        b"loop_\n_atom_site_anisotrop.id\n_atom_site_anisotrop.U[1][1]\n1 0.1\n"
        b"loop_\n_struct_conn.id\n_struct_conn.ptnr1_label_atom_id\n"
        b"_struct_conn.ptnr2_label_atom_id\nc1 H CA\n"
    )
    assert result.to_canonical() == expected.to_canonical()
    assert nibbler.mmcif.without_hydrogens(SOURCE).to_canonical() == (
        expected.to_canonical()
    )


def test_document_without_hydrogens_is_returned_itself():
    document = nibbler.cif.read(
        b"data_x\nloop_\n_atom_site.id\n_atom_site.type_symbol\n1 C\n2 O\n"
    )
    assert nibbler.mmcif.without_hydrogens(document) is document


def test_rejects_a_table_without_element_symbols():
    document = nibbler.cif.read(b"data_x\nloop_\n_atom_site.id\n1\n2\n")
    with pytest.raises(nibbler.ChemistryError) as raised:
        nibbler.mmcif.without_hydrogens(document)
    assert raised.value.code == "PDBX_HYDROGEN_REMOVAL"
