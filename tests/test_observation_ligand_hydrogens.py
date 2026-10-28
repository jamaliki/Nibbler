"""Keep ligand hydrogen capacity failures explicit and lossless."""

import json

import numpy as np
import pytest

import nibbler
from token_tables import token_table

observations = pytest.importorskip("nibbler.observations")
schema = pytest.importorskip("kaveh_data.atom32_schema")


def document(rows, connections=""):
    header = """data_ligand
loop_
_atom_site.group_PDB
_atom_site.label_atom_id
_atom_site.type_symbol
_atom_site.label_asym_id
_atom_site.label_comp_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
"""
    return nibbler.cif.read((header + rows + connections).encode())


def ligand_document(count):
    positions = [(0.97, 0, 0), (0, 0.97, 0), (0, 0, 0.97), (-0.97, 0, 0)]
    rows = "HETATM C4 C G LIG 1 . 101 . 0 0 0 1\nHETATM C3 C G LIG 1 . 101 . 5 0 0 1\n"
    rows += "".join(
        f"HETATM H{i + 1} H G LIG 1 . 101 . {x} {y} {z} 1\n"
        for i, (x, y, z) in enumerate(positions[:count])
    )
    return document(rows)


def ligand_packer():
    native = observations.ObservationPacker(schema.ENCODED, token_table("C"), [])
    native.add_component("LIG", {"C4": [0, 4, 3, 1], "C3": [0, 4, 0, 1]})
    return native


def arrays(result):
    buffers, keys = result
    return {
        name: np.frombuffer(buffer, dtype=dtype).reshape(shape)
        for name, (dtype, shape, buffer) in buffers.items()
    }, keys


# Exact ULW coordinates from the indexed 7NQC assembly after Reduce3.
ULW_ROWS = """HETATM C13 C G ULW 4 . 301 . 13.444 -11.316 2.962 1
HETATM C14 C G ULW 4 . 301 . 14.122 -11.071 4.273 1
HETATM H1 H G ULW 4 . 301 . 13.746 -10.904 3.786 1
HETATM H141 H G ULW 4 . 301 . 13.886 -11.756 4.918 1
HETATM H142 H G ULW 4 . 301 . 15.087 -11.076 4.170 1
HETATM H143 H G ULW 4 . 301 . 13.865 -10.211 4.640 1
"""
ULW_PARENTS = {"H1": "C13", "H141": "C14", "H142": "C14", "H143": "C14"}


def ulw_packer():
    native = observations.ObservationPacker(schema.ENCODED, token_table("C"), [])
    native.add_component("ULW", {"C13": [0, 4, 1, 1], "C14": [0, 4, 3, 1]})
    return native


def covalent_link(parent):
    return f"""loop_
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.ptnr2_label_atom_id
covale G LIG . H1 G LIG . {parent}
"""


def test_overflow_error_keeps_parent_and_all_assigned_hydrogens():
    with pytest.raises(
        ValueError, match="Ligand parent has more than three hydrogens"
    ) as caught:
        ligand_packer().pack(ligand_document(4), complete=True, validate=True)
    detail = json.loads(str(caught.value).split(": ", 1)[1])
    assert detail["parent_key"] == ["G", "101", ".", "LIG", "C4"]
    assert detail["parent_element"] == "C" and detail["parent_position"] == [0, 0, 0]
    assert detail["capacity"] == 3
    assert detail["ccd_features"]["total_hydrogen"] == 3
    assert [child["name"] for child in detail["hydrogens"]] == ["H1", "H2", "H3", "H4"]
    np.testing.assert_allclose(
        [child["distance"] for child in detail["hydrogens"]], 0.97
    )
    np.testing.assert_allclose(
        [child["position"] for child in detail["hydrogens"]],
        [[0.97, 0, 0], [0, 0.97, 0], [0, 0, 0.97], [-0.97, 0, 0]],
    )


def test_three_hydrogens_keep_the_existing_coordinates_and_capacity():
    buffers, keys = ligand_packer().pack(
        ligand_document(3), complete=True, validate=True
    )
    arrays = {
        name: np.frombuffer(buffer, dtype=dtype).reshape(shape)
        for name, (dtype, shape, buffer) in buffers.items()
    }
    assert keys[0] == ["G", "101", ".", "LIG", "C4"]
    assert arrays["mask"][0, 15:].sum() == 3
    np.testing.assert_allclose(
        arrays["pos"][0, 15:18], [[0.97, 0, 0], [0, 0.97, 0], [0, 0, 0.97]]
    )
    assert arrays["ligand_charge_known"][0]


def test_7nqc_ccd_bonds_fix_fourth_hydrogen_without_changing_coordinates():
    native = ulw_packer()
    with pytest.raises(ValueError, match="more than three hydrogens"):
        native.pack(document(ULW_ROWS), complete=True, validate=True)
    native.add_hydrogen_parents("ULW", ULW_PARENTS)
    actual, keys = arrays(native.pack(document(ULW_ROWS), complete=True, validate=True))
    assert keys == [["G", "301", ".", "ULW", "C13"], ["G", "301", ".", "ULW", "C14"]]
    np.testing.assert_array_equal(actual["mask"][:, 15:].sum(axis=1), [1, 3])
    np.testing.assert_array_equal(
        actual["pos"][:, 0],
        np.array([[13.444, -11.316, 2.962], [14.122, -11.071, 4.273]], np.float32),
    )
    np.testing.assert_array_equal(
        actual["pos"][0, 15], np.array([13.746, -10.904, 3.786], np.float32)
    )
    np.testing.assert_array_equal(
        actual["pos"][1, 15:18],
        np.array(
            [
                [13.886, -11.756, 4.918],
                [15.087, -11.076, 4.170],
                [13.865, -10.211, 4.640],
            ],
            np.float32,
        ),
    )
    assert actual["ligand_charge_known"].all()


@pytest.mark.parametrize(
    ("name", "element", "ccd_name"),
    [("D1", "D", "H1"), ("T1", "T", "H1"), ("1D1", "D", "1H1")],
)
def test_ccd_hydrogen_parents_accept_isotope_aliases(name, element, ccd_name):
    native = ulw_packer()
    native.add_hydrogen_parents("ULW", {**ULW_PARENTS, ccd_name: "C13"})
    actual, _ = arrays(
        native.pack(
            document(ULW_ROWS.replace("H1 H G", f"{name} {element} G")),
            complete=True,
            validate=True,
        )
    )
    np.testing.assert_array_equal(actual["mask"][:, 15:].sum(axis=1), [1, 3])
    assert actual["ligand_charge_known"].all()


def test_observed_covalent_parent_overrides_ccd_parent():
    rows = (
        "HETATM C4 C G LIG 1 . 101 . 0 0 0 1\nHETATM C3 C G LIG 1 . 101 . 1.5 0 0 1\n"
    )
    rows += "HETATM H1 H G LIG 1 . 101 . .6 0 0 1\n"
    native = ligand_packer()
    native.add_hydrogen_parents("LIG", {"H1": "C4"})
    actual, _ = arrays(
        native.pack(document(rows, covalent_link("C3")), complete=True, validate=True)
    )
    np.testing.assert_array_equal(actual["mask"][:, 15:].sum(axis=1), [0, 1])
    np.testing.assert_array_equal(
        actual["pos"][1, 15], np.array([0.6, 0, 0], np.float32)
    )


@pytest.mark.parametrize("source", ["CCD", "struct_conn"])
def test_absent_declared_parent_fails_instead_of_using_nearest_atom(source):
    native = ligand_packer()
    if source == "CCD":
        native.add_hydrogen_parents("LIG", {"H1": "ABSENT"})
        doc = ligand_document(1)
    else:
        rows = "HETATM C4 C G LIG 1 . 101 . 0 0 0 1\nHETATM H1 H G LIG 1 . 101 . .97 0 0 1\n"
        doc = document(rows, covalent_link("ABSENT"))
    with pytest.raises(ValueError, match="Invalid ligand hydrogen parent") as caught:
        native.pack(doc, complete=True, validate=True)
    detail = json.loads(str(caught.value).split(": ", 1)[1])
    assert detail["hydrogen_key"] == ["G", "101", ".", "LIG", "H1"]
    assert "absent" in detail["reason"]


def test_invalid_ccd_parent_distance_fails_without_reassignment():
    native = ligand_packer()
    native.add_hydrogen_parents("LIG", {"H1": "C3"})
    with pytest.raises(ValueError, match="bonded parent has invalid geometry"):
        native.pack(ligand_document(1), complete=True, validate=True)


def test_ccd_parent_resolves_ambiguous_distance_without_moving_hydrogen():
    rows = (
        "HETATM C4 C G LIG 1 . 101 . 0 0 0 1\nHETATM C3 C G LIG 1 . 101 . 1.8 0 0 1\n"
    )
    rows += "HETATM H1 H G LIG 1 . 101 . .9 0 0 1\n"
    native = ligand_packer()
    native.add_hydrogen_parents("LIG", {"H1": "C4"})
    actual, _ = arrays(native.pack(document(rows), complete=True, validate=True))
    np.testing.assert_array_equal(actual["mask"][:, 15:].sum(axis=1), [1, 0])
    np.testing.assert_array_equal(
        actual["pos"][0, 15], np.array([0.9, 0, 0], np.float32)
    )
    assert actual["ligand_hydrogen_known"].all()


def test_conflicting_observed_parents_fail_with_hydrogen_identity():
    links = covalent_link("C4") + "covale G LIG . H1 G LIG . C3\n"
    rows = (
        "HETATM C4 C G LIG 1 . 101 . 0 0 0 1\nHETATM C3 C G LIG 1 . 101 . 1.8 0 0 1\n"
    )
    rows += "HETATM H1 H G LIG 1 . 101 . .9 0 0 1\n"
    with pytest.raises(ValueError, match="conflicting covalent parents") as caught:
        ligand_packer().pack(document(rows, links), complete=True, validate=True)
    detail = json.loads(str(caught.value).split(": ", 1)[1])
    assert detail["hydrogen_key"] == ["G", "101", ".", "LIG", "H1"]


def test_observed_parent_can_belong_to_another_ligand_occurrence():
    rows = (
        "HETATM C4 C G LIG 1 . 101 . 0 0 0 1\nHETATM C3 C Q OTHER 2 . 202 . 1.5 0 0 1\n"
    )
    rows += "HETATM H1 H G LIG 1 . 101 . .6 0 0 1\n"
    native = ligand_packer()
    native.add_hydrogen_parents("LIG", {"H1": "C4"})
    links = covalent_link("C3").replace("G LIG . C3", "Q OTHER . C3")
    actual, keys = arrays(
        native.pack(document(rows, links), complete=True, validate=True)
    )
    assert keys[1] == ["Q", "202", ".", "OTHER", "C3"]
    np.testing.assert_array_equal(actual["mask"][:, 15:].sum(axis=1), [0, 1])
    np.testing.assert_array_equal(
        actual["pos"][1, 15], np.array([0.6, 0, 0], np.float32)
    )
