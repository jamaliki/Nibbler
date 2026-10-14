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
