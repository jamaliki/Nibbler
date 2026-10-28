"""Match nucleotide dictionary names without moving atoms or changing labels."""

import os
from pathlib import Path

import numpy as np
import pytest
from token_tables import token_table

import nibbler
from nibbler.observations import ObservationPacker

HEADER = (
    "data_nucleotide\nloop_\n_atom_site.group_PDB\n_atom_site.label_atom_id\n"
    "_atom_site.type_symbol\n_atom_site.label_asym_id\n_atom_site.label_comp_id\n"
    "_atom_site.label_entity_id\n_atom_site.label_seq_id\n_atom_site.auth_seq_id\n"
    "_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n"
    "_atom_site.Cartn_z\n_atom_site.pdbx_PDB_model_num\n"
)


def packer(parent="C5M", child="H5M1", component="DT"):
    file = os.environ.get("KAVEH_ATOM32_SCHEMA")
    encoded = (
        Path(file).read_text()
        if file
        else pytest.importorskip("kaveh_data.atom32_schema").ENCODED
    )
    native = ObservationPacker(encoded, token_table("C", "O"), [])
    native.add_component(
        component,
        {
            name: [0, 4, int(name != "C2"), 1]
            for name in ("C7", "C5M", "C1*", "C1'", "C2")
        },
    )
    native.add_hydrogen_parents(component, {child: parent})
    return native


def row(name, element, x, *, chain="A", component="DT"):
    return f"HETATM {name} {element} {chain} {component} 1 . 1 ? {x} 0 0 1\n"


def pack(native, rows):
    buffers, keys = native.pack(
        nibbler.cif.read((HEADER + "".join(rows)).encode()),
        complete=True,
        validate=True,
    )
    return {
        name: np.frombuffer(data, dtype=dtype).reshape(shape)
        for name, (dtype, shape, data) in buffers.items()
    }, keys


@pytest.mark.parametrize("parent", ["C7", "C5M"])
@pytest.mark.parametrize("child,element", [("H5M1", "H"), ("H71", "H"), ("D71", "D")])
def test_thymine_methyl_alias_keeps_parent_coordinates_and_labels(
    parent, child, element
):
    actual, keys = pack(
        packer(), [row(parent, "C", 0), row("C2", "C", 0.8), row(child, element, 0.97)]
    )
    index = next(i for i, key in enumerate(keys) if key[-1] == parent)
    other = next(i for i, key in enumerate(keys) if key[-1] == "C2")
    assert (
        actual["mask"][index, 15:].sum() == 1 and not actual["mask"][other, 15:].any()
    )
    np.testing.assert_array_equal(
        actual["pos"][index, 15], np.asarray([0.97, 0, 0], np.float32)
    )
    np.testing.assert_array_equal(actual["pos"][index, 0], np.zeros(3, np.float32))
    assert actual["ligand_charge_known"].all() and actual["ligand_hydrogen_known"].all()


def test_one_component_accepts_mixed_spellings_across_residues():
    actual, keys = pack(
        packer(),
        [
            row("C7", "C", 0),
            row("H5M1", "H", 0.97),
            row("C5M", "C", 3, chain="B"),
            row("H5M1", "H", 3.97, chain="B"),
        ],
    )
    assert [key[-1] for key in keys] == ["C7", "C5M"]
    np.testing.assert_array_equal(actual["mask"][:, 15:].sum(axis=1), [1, 1])
    np.testing.assert_array_equal(
        actual["pos"][:, 15, 0], np.asarray([0.97, 3.97], np.float32)
    )


def test_prime_alias_keeps_the_dictionary_parent_when_another_atom_is_closer():
    actual, keys = pack(
        packer(parent="C1'", child="H1'", component="DA"),
        [
            row("C1*", "C", 0, component="DA"),
            row("C2", "C", 0.8, component="DA"),
            row("H1'", "H", 0.97, component="DA"),
        ],
    )
    index = next(i for i, key in enumerate(keys) if key[-1] == "C1*")
    assert actual["mask"][index, 15:].sum() == 1


def test_alias_does_not_relax_geometry_or_missing_parent_checks():
    with pytest.raises(ValueError, match="invalid geometry"):
        pack(packer(), [row("C7", "C", 3), row("C2", "C", 0), row("H5M1", "H", 0.97)])
    with pytest.raises(ValueError, match="Dictionary parent is absent"):
        pack(packer(), [row("C2", "C", 0), row("H5M1", "H", 0.97)])
    with pytest.raises(ValueError, match="Dictionary parent is absent"):
        pack(
            packer(component="LIG"),
            [
                row("C7", "C", 0, component="LIG"),
                row("H5M1", "H", 0.97, component="LIG"),
            ],
        )


def test_ambiguous_alias_and_conflicting_hydrogen_aliases_fail():
    with pytest.raises(ValueError, match="parent alias is ambiguous"):
        pack(packer(), [row("C7", "C", 0), row("c7", "C", 0.1), row("H5M1", "H", 0.97)])
    native = packer(child="H71")
    native.add_hydrogen_parents("DT", {"H71": "C5M", "H5M1": "C2"})
    with pytest.raises(ValueError, match="conflicting parents"):
        pack(native, [row("C7", "C", 0), row("C2", "C", 0.8), row("h5m1", "H", 0.97)])
