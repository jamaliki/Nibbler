"""Retain full source atoms before the caller selects a dense token layout."""

import numpy as np
import pytest
from test_observation_ligand_hydrogens import covalent_link, document, ligand_packer


def unpack(source):
    buffers, labels = source.atom_table()
    return {
        name: np.frombuffer(data, dtype=dtype).reshape(shape)
        for name, (dtype, shape, data) in buffers.items()
    }, labels


def test_atom_table_retains_unknown_protein_sidechains_and_hydrogens():
    doc = document(
        "ATOM N N A MSE 1 1 1 . 0 0 0 1\n"
        "ATOM SE SE A MSE 1 1 1 . 1 0 0 1\n"
        "ATOM HE H A MSE 1 1 1 . 2 0 0 1\n"
        "HETATM P P B DA 2 1 8 . 3 0 0 1\n"
        "HETATM H1 H B DA 2 1 8 . 4 0 0 1\n"
        "loop_\n_entity_poly.entity_id\n_entity_poly.type\n"
        "1 'polypeptide(L)'\n2 polydeoxyribonucleotide\n"
    )
    source = ligand_packer().prepare_source(doc)
    arrays, labels = unpack(source)
    np.testing.assert_array_equal(arrays["source_atom_row"], np.arange(5))
    np.testing.assert_array_equal(arrays["source_atom_group"], [0, 0, 0, 1, 1])
    np.testing.assert_array_equal(arrays["source_atom_is_protein"], [1, 1, 1, 0, 0])
    assert labels["source_atom_element"] == ["N", "SE", "H", "P", "H"]
    assert (
        labels["source_atom_polymer_type"]
        == ["polypeptide(L)"] * 3 + ["polydeoxyribonucleotide"] * 2
    )
    assert arrays["source_atom_mask"].all()
    assert arrays["source_atom_occupancy"].tolist() == [1] * 5
    arrays["source_atom_pos"][0] = 100
    assert unpack(source)[0]["source_atom_pos"][0, 0] == 0


def test_atom_table_preserves_atom_endpoints_type_and_order():
    doc = document(
        "HETATM C4 C G LIG 1 . 101 . 0 0 0 1\n"
        "HETATM C3 C G LIG 1 . 101 . 1 0 0 1\n"
        "HETATM H1 H G LIG 1 . 101 . 2 0 0 1\n",
        covalent_link("C3"),
    )
    arrays, labels = unpack(ligand_packer().prepare_source(doc))
    np.testing.assert_array_equal(arrays["source_connection_index"], [[2], [1]])
    np.testing.assert_array_equal(arrays["source_hydrogen_parent"], [-1, -1, 1])
    assert labels["source_connection_type"] == ["covale"]
    assert labels["source_connection_order"] == ["."]


def crop_policy(minimum=None, protein=None):
    count = 80
    pos = np.zeros((count, 3), dtype=np.float32)
    pos[:, 0] = np.arange(count)
    chain = np.repeat([0, 1], [60, 20]).astype(np.int64)
    mask = np.zeros(count, dtype=bool)
    kwargs = (
        {}
        if minimum is None
        else dict(
            chain_minimum=np.array(minimum, np.int64).tobytes(),
            chain_protein=np.array(protein, bool).tobytes(),
        )
    )
    values, _ = ligand_packer().crop_plans(
        pos.tobytes(),
        chain.tobytes(),
        mask.tobytes(),
        (~mask).tobytes(),
        crop_size=80,
        **kwargs,
    )
    return {
        name: np.frombuffer(data, dtype=dtype).reshape(shape)
        for name, (dtype, shape, data) in values.items()
    }


def test_typed_crop_policy_retains_short_nucleotide_chain():
    assert set(crop_policy()["crop_chain_bits"]) == {1}
    assert set(crop_policy([50, 1], [1, 0])["crop_chain_bits"]) == {3}
    assert len(crop_policy([1, 1], [0, 0])["crop_center"]) == 0


@pytest.mark.parametrize(
    "minimum,protein", [([1, 1], [1, 0]), ([50], [1]), ([50, 0], [1, 0])]
)
def test_typed_crop_policy_rejects_invalid_chain_contract(minimum, protein):
    with pytest.raises(ValueError):
        crop_policy(minimum, protein)
