"""Clash-aware views and special-position ligands use the native path."""

import pyarrow as pa
import pytest

import nibbler


def rows(document, category="atom_site"):
    return (
        pa.RecordBatchReader.from_stream(
            nibbler.chomp(document.to_canonical().encode(), category=category)
        )
        .read_all()
        .to_pylist()
    )


@pytest.mark.parametrize("case", ["1cbn", "6lzb", "1rcq"])
def test_real_conformer_selection_preserves_chemical_identity(case, preparation_data):
    document = nibbler.chomp(preparation_data.entry(case))
    before = document.to_canonical()
    selection = nibbler.mmcif.select_conformer(document, seed=42)
    repeat = nibbler.mmcif.select_conformer(document, seed=42)
    assert selection.atom_rows == repeat.atom_rows
    original = rows(document)
    selected = rows(selection.document)
    assert selected == [original[i] for i in selection.atom_rows]
    assert document.to_canonical() == before
    grouped = {}
    for a in selected:
        key = (a["label_asym_id"], a["auth_seq_id"], a.get("pdbx_PDB_ins_code"))
        if a["label_alt_id"] not in {".", "?", None}:
            grouped.setdefault(key, set()).add(a["label_alt_id"])
    assert all(len(alts) == 1 for alts in grouped.values())


def test_modified_residue_requires_explicit_fallback_policy(preparation_data):
    document = nibbler.chomp(preparation_data.entry("6qhp"))
    before = document.to_canonical()
    with pytest.raises(nibbler.ChemistryError, match="clash-free conformer"):
        nibbler.mmcif.select_conformer(document, seed=42)
    assert document.to_canonical() == before


@pytest.mark.parametrize("case,ligand_chains", [("7mub", {"D", "E"}), ("1xan", {"C"})])
def test_special_position_ligands_keep_one_complete_source_copy(
    case, ligand_chains, preparation_data
):
    document = nibbler.chomp(preparation_data.entry(case))
    expanded, aliases = nibbler.mmcif.assembly_with_mapping(document)
    assert aliases and all(
        removed != retained and retained in ligand_chains
        for _, removed, retained in aliases
    )
    assert expanded.to_canonical() == nibbler.mmcif.assembly(document).to_canonical()
    chains = {r["label_asym_id"] for r in rows(expanded)}
    assert ligand_chains <= chains
    assert not any(
        chain.startswith(tuple(name + "-" for name in ligand_chains))
        for chain in chains
    )
    source_codes = {
        error.code
        for error in nibbler.sniff(nibbler.mmcif.read(document), profile="pdbx").errors
    }
    result_codes = {
        error.code
        for error in nibbler.sniff(nibbler.mmcif.read(expanded), profile="pdbx").errors
    }
    # These deposited snapshots use an older dictionary and omit water schemes.
    assert result_codes <= source_codes


def test_array_selection_matches_document_selection(preparation_data):
    source = nibbler.chomp(preparation_data.entry("1cbn"))
    atoms = rows(source)
    group_ids = {}
    residue_ids = []
    for a in atoms:
        key = (a["label_asym_id"], a["auth_seq_id"], a.get("pdbx_PDB_ins_code"))
        residue_ids.append(group_ids.setdefault(key, len(group_ids)))
    selected = nibbler.mmcif.select_atom_conformers(
        [[float(a[k]) for k in ("Cartn_x", "Cartn_y", "Cartn_z")] for a in atoms],
        [float(a["occupancy"]) for a in atoms],
        residue_ids,
        [a["label_asym_id"] for a in atoms],
        [
            int(a["label_seq_id"])
            if a["label_seq_id"] not in {None, ".", "?"}
            else None
            for a in atoms
        ],
        [a["label_atom_id"] for a in atoms],
        [
            a["label_alt_id"] if a["label_alt_id"] not in {None, ".", "?"} else None
            for a in atoms
        ],
        bonds=[],
        clash_mask=[a["type_symbol"].upper() not in {"H", "D", "T"} for a in atoms],
    )
    assert selected.atom_rows == nibbler.mmcif.select_conformer(source).atom_rows
