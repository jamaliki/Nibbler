"""Seeded random selection keeps source identity, models and qualified links."""

import pytest

import nibbler


def document(rows, links="", *, insertion=False):
    names = (
        "group_PDB id type_symbol label_atom_id label_alt_id label_comp_id "
        "label_asym_id label_entity_id label_seq_id Cartn_x Cartn_y Cartn_z occupancy "
        "auth_seq_id auth_comp_id auth_asym_id auth_atom_id pdbx_PDB_model_num"
    ).split()
    header = "data_choices\nloop_\n" + "\n".join("_atom_site." + name for name in names)
    if insertion:
        header += "\n_atom_site.pdbx_PDB_ins_code"
    return nibbler.cif.read((header + "\n" + "\n".join(rows) + "\n" + links).encode())


def test_random_policy_varies_unequal_occupancies_without_changing_source():
    source = document(
        [
            "ATOM 1 C C1 A LIG A 1 . 0 0 0 .95 1 LIG A C1 1",
            "ATOM 2 C C1 B LIG A 1 . 10 0 0 .05 1 LIG A C1 1",
            "ATOM 3 C C1 . LIG A 1 . 20 0 0 1 1 LIG A C1 3",
        ]
    )
    original, choices = source.to_canonical(), set()
    for seed in range(32):
        first = nibbler.mmcif.select_conformer(
            source, seed=seed, policy="random_clash_aware"
        )
        replay = nibbler.mmcif.select_conformer(
            source, seed=seed, policy="random_clash_aware"
        )
        assert first.atom_rows == replay.atom_rows
        assert first.schema_version == "random-clash-aware-altloc.v1"
        choices.add(tuple(first.atom_rows))
        assert nibbler.mmcif.select_conformer(source, seed=seed).atom_rows == [0]
    assert choices == {(0,), (1,)}
    assert nibbler.mmcif.select_conformer(
        source, model_number=3, policy="random_clash_aware"
    ).atom_rows == [2]
    assert source.to_canonical() == original


def test_random_policy_keeps_explicit_link_alternatives_together():
    names = (
        "id conn_type_id ptnr1_label_asym_id ptnr1_label_atom_id ptnr2_label_asym_id "
        "ptnr2_label_atom_id pdbx_ptnr1_label_alt_id pdbx_ptnr2_label_alt_id"
    ).split()
    links = "loop_\n" + "\n".join("_struct_conn." + name for name in names)
    source = document(
        [
            "ATOM 1 C C1 A LIG A 1 . 0 0 0 .95 1 LIG A C1 1",
            "ATOM 2 C C1 B LIG A 1 . 10 0 0 .05 1 LIG A C1 1",
            "ATOM 3 C C1 A LIG B 2 . 3 0 0 .05 1 LIG B C1 1",
            "ATOM 4 C C1 B LIG B 2 . 13 0 0 .95 1 LIG B C1 1",
        ],
        links + "\na covale A C1 B C1 A A\nb covale A C1 B C1 B B\n",
    )
    choices = {
        tuple(
            nibbler.mmcif.select_conformer(
                source, seed=seed, policy="random_clash_aware"
            ).atom_rows
        )
        for seed in range(32)
    }
    assert choices == {(0, 2), (1, 3)}


def test_unknown_policy_fails():
    source = document(["ATOM 1 C C1 . LIG A 1 . 0 0 0 1 1 LIG A C1 1"])
    with pytest.raises(ValueError, match="policy"):
        nibbler.mmcif.select_conformer(source, policy="guessed")


def test_link_endpoint_insertion_codes_do_not_correlate_another_residue():
    names = (
        "id conn_type_id ptnr1_label_asym_id ptnr1_label_atom_id ptnr2_label_asym_id "
        "ptnr2_label_atom_id pdbx_ptnr1_label_alt_id pdbx_ptnr2_label_alt_id "
        "pdbx_ptnr1_PDB_ins_code"
    ).split()
    links = "loop_\n" + "\n".join("_struct_conn." + name for name in names)
    source = document(
        [
            "ATOM 1 C C1 A LIG A 1 . 0 0 0 .5 1 LIG A C1 1 A",
            "ATOM 2 C C1 B LIG A 1 . 10 0 0 .5 1 LIG A C1 1 A",
            "ATOM 3 C C1 A LIG A 1 . 100 0 0 .5 1 LIG A C1 1 B",
            "ATOM 4 C C1 B LIG A 1 . 110 0 0 .5 1 LIG A C1 1 B",
            "ATOM 5 C C1 A LIG L 2 . 2 0 0 .5 1 LIG L C1 1 .",
            "ATOM 6 C C1 B LIG L 2 . 12 0 0 .5 1 LIG L C1 1 .",
        ],
        links + "\na covale A C1 L C1 A A A\nb covale A C1 L C1 B B A\n",
        insertion=True,
    )
    choices = {
        tuple(
            nibbler.mmcif.select_conformer(
                source, seed=seed, policy="random_clash_aware"
            ).atom_rows
        )
        for seed in range(64)
    }
    assert choices == {(0, 2, 4), (0, 3, 4), (1, 2, 5), (1, 3, 5)}
