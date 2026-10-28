"""Ordinary molecular views retain spatial crops without metal-hydride bridges."""

import numpy as np
import pytest

import nibbler
from nibbler.observations import ObservationPacker
from token_tables import frozen, token_table


def complex_document():
    raw = """data_crop
    _entry.id CROP
    loop_
    _entity_poly.entity_id
    _entity_poly.type
    1 'polypeptide(L)'
    loop_
    _atom_site.group_PDB
    _atom_site.id
    _atom_site.label_entity_id
    _atom_site.label_asym_id
    _atom_site.label_seq_id
    _atom_site.auth_seq_id
    _atom_site.label_comp_id
    _atom_site.label_atom_id
    _atom_site.type_symbol
    _atom_site.Cartn_x
    _atom_site.Cartn_y
    _atom_site.Cartn_z
    _atom_site.pdbx_PDB_model_num
    _atom_site.occupancy
    _atom_site.auth_comp_id
    _atom_site.auth_asym_id
    _atom_site.auth_atom_id
    _atom_site.label_alt_id
    _atom_site.pdbx_PDB_ins_code
    _atom_site.B_iso_or_equiv
    """
    for residue in range(60):
        names = [("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O"), ("CB", "C")]
        for slot, (name, element) in enumerate(names):
            raw += (
                f"ATOM {5 * residue + slot + 1} 1 A {residue + 1} {residue + 1} "
                f"ALA {name} {element} {3 * residue + slot * 0.2 + 1} 2 3 "
                f"1 1 ALA A {name} . ? 10\n"
            )
    raw += (
        "HETATM 301 2 B . 1 LIG C1 C 178 3 3 1 1 LIG B C1 . ? 10\n"
        "HETATM 302 2 B . 1 LIG C2 C 179 3 3 1 1 LIG B C2 . ? 10\n"
    )
    return nibbler.cif.read(raw.encode())


def arrays(buffers):
    return {
        name: np.frombuffer(data, dtype=dtype).reshape(shape)
        for name, (dtype, shape, data) in buffers.items()
    }


@pytest.mark.parametrize("index", [False, True])
def test_molecular_crop_keeps_the_ligand_center(tmp_path, index):
    # A minimal library keeps this regression independent of external chemistry data.
    (tmp_path / "mon_lib/list").mkdir(parents=True)
    (tmp_path / "geostd").mkdir()
    (tmp_path / "mon_lib/list/mon_lib_list.cif").write_text("data_empty\n")
    (tmp_path / "geostd/ener_lib.cif").write_text("data_empty\n")
    packer = ObservationPacker(
        frozen("atom32_schema.json.gz"), token_table("C", "O", "N"), []
    )
    packer.set_molecular_schema(frozen("molecular_schema.json.gz"))
    doc = complex_document()
    legacy, _ = packer.pack(doc)
    expected = arrays(legacy)["pos"][:, 0]
    assert expected.shape == (62, 3) and np.all(expected != 0)
    assert "ligand_bridge_index" not in legacy
    view = packer.molecular_view(doc, chem_data=tmp_path)
    perceived = [bytes(32 * len(graph["atoms"])) for graph in view.perception_graphs()]
    buffers, _, _, _ = packer.finish_molecular(view, perceived, sasa=False, index=index)
    result = arrays(buffers)
    np.testing.assert_array_equal(result["crop_positions"], expected)
    assert result["crop_center_eligible"].all()
    plans, _ = packer.crop_plans(
        result["crop_positions"].tobytes(),
        result["chain_idx"].tobytes(),
        result["is_ligand"].tobytes(),
        (result["is_ligand"] | result["mask"][:, 1]).tobytes(),
        crop_size=52,
    )
    plans = arrays(plans)
    ligand_plan = np.flatnonzero(plans["crop_center"] == 60).item()
    assert plans["crop_ligand_tokens"][ligand_plan] == 2
    assert plans["crop_radius"][ligand_plan] > 0
