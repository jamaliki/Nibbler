"""Keep canonical hydrogen parents when another atom is closer."""

import json
import os
from pathlib import Path

import nibbler
import numpy as np
import pytest
from nibbler.observations import ObservationPacker


CYS_7F52 = [
    ("N", "N", (-4.965, 33.169, -45.688)),
    ("CA", "C", (-6.159, 33.997, -45.394)),
    ("C", "C", (-7.353, 34.952, -45.412)),
    ("O", "O", (-7.365, 35.961, -44.680)),
    ("CB", "C", (-5.768, 35.390, -44.918)),
    ("SG", "S", (-7.078, 36.230, -43.991)),
    ("H", "H", (-4.945, 32.454, -45.210)),
    ("HA", "H", (-6.439, 33.073, -45.487)),
    ("HB2", "H", (-4.994, 35.317, -44.338)),
    ("HB3", "H", (-5.554, 35.936, -45.691)),
    ("HG", "H", (-7.817, 36.773, -44.765)),
]


ASP_1R24 = [
    ("N", "N", (75.097, -21.262, 8.356)),
    ("CA", "C", (76.052, -20.122, 8.423)),
    ("C", "C", (76.202, -21.086, 9.488)),
    ("O", "O", (75.551, -22.048, 9.016)),
    ("CB", "C", (75.224, -18.914, 8.769)),
    ("CG", "C", (74.850, -18.615, 7.355)),
    ("OD1", "O", (75.372, -19.394, 6.549)),
    ("OD2", "O", (74.117, -17.664, 7.060)),
    ("H1", "H", (75.534, -22.024, 8.501)),
    ("H2", "H", (74.467, -21.158, 8.976)),
    ("H3", "H", (74.719, -21.286, 7.550)),
    ("HA", "H", (76.728, -19.722, 7.854)),
    ("HB2", "H", (74.454, -19.123, 9.321)),
    ("HB3", "H", (75.740, -18.201, 9.176)),
]


def packer():
    file = os.environ.get("KAVEH_ATOM32_SCHEMA")
    if file:
        encoded = Path(file).read_text()
    else:
        encoded = pytest.importorskip("kaveh_data.atom32_schema").ENCODED
    return ObservationPacker(encoded, ["S"], []), json.loads(encoded)


def document(child="HG", element="H", linked=False):
    header = ("data_parent\nloop_\n_entity_poly.entity_id\n_entity_poly.type\n1 'polypeptide(L)'\n"
              "loop_\n_entity_poly_seq.entity_id\n_entity_poly_seq.num\n1 570\nloop_\n"
              "_atom_site.group_PDB\n_atom_site.label_entity_id\n_atom_site.label_asym_id\n"
              "_atom_site.label_seq_id\n_atom_site.auth_seq_id\n_atom_site.label_comp_id\n"
              "_atom_site.label_atom_id\n_atom_site.type_symbol\n_atom_site.Cartn_x\n"
              "_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.pdbx_PDB_model_num\n")
    rows = "".join(f"ATOM 1 A 569 569 CYS {child if name == 'HG' else name} "
                   f"{element if name == 'HG' else kind} {x} {y} {z} 1\n"
                   for name, kind, (x, y, z) in CYS_7F52)
    if linked:
        rows += "HETATM 2 B . 1 LIG S S -5.048 36.230 -43.991 1\n#\nloop_\n"
        rows += "".join(f"_struct_conn.{name}\n" for name in (
            "conn_type_id", "ptnr1_label_asym_id", "ptnr1_label_seq_id", "ptnr1_label_comp_id",
            "ptnr1_label_atom_id", "ptnr2_label_asym_id", "ptnr2_label_seq_id", "ptnr2_label_comp_id",
            "ptnr2_label_atom_id", "ptnr1_symmetry", "ptnr2_symmetry"))
        rows += "disulf A 569 CYS SG B . LIG S 1_555 1_555\n"
    return nibbler.cif.read((header + rows).encode())


@pytest.mark.parametrize("name,element", [("HG", "H"), ("HG", "D"), ("DG", "D"), ("TG", "T")])
def test_7f52_named_hydrogen_stays_on_sulfur(name, element):
    native, schema = packer()
    buffers, _ = native.pack(document(name, element), complete=True, protein_only=True,
                            topology_policy="entity_sequence", validate=True)
    result = {key: np.frombuffer(data, dtype=dtype).reshape(shape)
              for key, (dtype, shape, data) in buffers.items()}
    cys = schema["residues"][4]
    assert result["mask"][0, 21] and not result["mask"][0, 22]
    assert result["state_known"][0] and result["chemical_state"][0] == 0
    assert result["atom_presence_known"][0].all()
    np.testing.assert_array_equal(result["mask"][0], cys["states"][0]["present"])
    for atom, _, xyz in CYS_7F52:
        slot = cys["names"].index(atom)
        np.testing.assert_array_equal(result["pos"][0, slot], np.asarray(xyz, dtype=np.float32))


def test_named_parent_does_not_relax_disulfide_validation():
    native, _ = packer()
    with pytest.raises(ValueError, match="Atom-presence labels contradict every allowed chemical state"):
        native.pack(document(linked=True), complete=True, validate=True)


def terminal_document(child="H1", element="H", distant_n=False):
    header = ("data_terminal_parent\nloop_\n_entity_poly.entity_id\n_entity_poly.type\n"
              "1 'polypeptide(L)'\nloop_\n_entity_poly_seq.entity_id\n_entity_poly_seq.num\n"
              "1 1\n1 2\nloop_\n_atom_site.group_PDB\n_atom_site.label_entity_id\n"
              "_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.auth_seq_id\n"
              "_atom_site.label_comp_id\n_atom_site.label_atom_id\n_atom_site.type_symbol\n"
              "_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n"
              "_atom_site.pdbx_PDB_model_num\n")
    rows = "".join(f"ATOM 1 C 1 1 ASP {child if name == 'H1' else name} "
                   f"{element if name == 'H1' else kind} "
                   f"{x + 5 if distant_n and name == 'N' else x} {y} {z} 1\n"
                   for name, kind, (x, y, z) in ASP_1R24)
    return nibbler.cif.read((header + rows).encode())


@pytest.mark.parametrize("name,element", [("H1", "H"), ("D1", "D"), ("T1", "T")])
def test_1r24_terminal_hydrogen_stays_on_nitrogen(name, element):
    native, schema = packer()
    source = {atom: np.asarray(xyz) for atom, _, xyz in ASP_1R24}
    assert np.linalg.norm(source["H1"] - source["O"]) < np.linalg.norm(source["H1"] - source["N"])
    buffers, _ = native.pack(terminal_document(name, element), complete=True, protein_only=True,
                            topology_policy="entity_sequence", validate=True)
    result = {key: np.frombuffer(data, dtype=dtype).reshape(shape)
              for key, (dtype, shape, data) in buffers.items()}
    asp = schema["residues"][3]
    assert result["mask"][0, 15:18].all() and not result["mask"][0, 23]
    assert result["state_known"][0] and result["atom_presence_known"][0].all()
    state = asp["states"][result["chemical_state"][0]]
    assert state["n_terminal"] and not state["c_terminal"]
    np.testing.assert_array_equal(result["mask"][0], state["present"])
    terminal = sorted((name if atom == "H1" else atom, xyz)
                      for atom, _, xyz in ASP_1R24 if atom in {"H1", "H2", "H3"})
    np.testing.assert_array_equal(result["pos"][0, 15:18],
                                  np.asarray([xyz for _, xyz in terminal], dtype=np.float32))
    for atom, _, xyz in ASP_1R24:
        if atom not in {"H1", "H2", "H3"}:
            slot = asp["names"].index(atom)
            np.testing.assert_array_equal(result["pos"][0, slot], np.asarray(xyz, dtype=np.float32))


@pytest.mark.parametrize("name,element", [("H1", "H"), ("D1", "D"), ("T1", "T")])
def test_terminal_hydrogen_rejects_distant_nitrogen_without_oxygen_fallback(name, element):
    native, _ = packer()
    doc = terminal_document(name, element, distant_n=True)
    buffers, _ = native.pack(doc, complete=False, protein_only=True,
                            topology_policy="entity_sequence", validate=False)
    dtype, shape, data = buffers["mask"]
    mask = np.frombuffer(data, dtype=dtype).reshape(shape)
    assert not mask[0, 15:18].any() and not mask[0, 23]
    buffers, _ = native.pack(doc, complete=True, protein_only=True,
                             topology_policy="entity_sequence", validate=True)
    arrays = {key: np.frombuffer(data, dtype=dtype).reshape(shape)
              for key, (dtype, shape, data) in buffers.items()}
    assert not arrays["mask"][0, 15:18].any() and not arrays["mask"][0, 23]
    assert not arrays["state_known"][0] and not arrays["atom_presence_known"][0, 15:18].any()
