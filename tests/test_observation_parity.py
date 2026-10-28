"""Compare native observations with the established training target contract."""

import json
from unittest.mock import patch

import numpy as np
import pytest
from token_tables import token_table

import nibbler

schema = pytest.importorskip("kaveh_data.atom32_schema")
atom32 = pytest.importorskip("kaveh_data.atom32")
fixtures = pytest.importorskip("test_pdb_ligands")
observations = pytest.importorskip("nibbler.observations")


def packer(catalog=None):
    from kaveh_data.ap_novo_pdb import EXCLUDED_COMPONENTS
    from kaveh_data.chemistry import AtomFeature
    from kaveh_data.pdb_tokens import LIGAND_ELEMENTS

    result = observations.ObservationPacker(
        schema.ENCODED, token_table(*LIGAND_ELEMENTS), sorted(EXCLUDED_COMPONENTS)
    )
    if catalog is not None:
        indices = [
            AtomFeature.FORMAL_CHARGE,
            AtomFeature.HYBRIDIZATION,
            AtomFeature.TOTAL_HYDROGEN,
            AtomFeature.FORMAL_CHARGE_KNOWN,
        ]
        for component in catalog.components:
            result.add_component(
                component.component_id,
                {
                    str(name): row[indices].tolist()
                    for name, row in zip(
                        (atom.atom_id for atom in component.atoms),
                        component.graph.atom_features,
                        strict=True,
                    )
                },
            )
    return result


def arrays(result):
    buffers, keys = result
    return {
        name: np.frombuffer(buffer, dtype=dtype).reshape(shape)
        for name, (dtype, shape, buffer) in buffers.items()
    }, keys


def product_options(document, catalog):
    from kaveh_data.product import prepare_product

    product = prepare_product(document, catalog)
    return dict(
        product_targets=json.dumps([] if product is None else product.targets())
    )


def legacy_tokens(raw, catalog, complete):
    from kaveh_data import ap_novo_pdb
    from kaveh_data.pdb_dataset import StoredPdbView
    from kaveh_data.pdb_structure import parse_pdb_structure
    from kaveh_data.pdb_token_view import prepare_token_view
    from kaveh_data.pdb_tokens import tokenize

    structure = parse_pdb_structure(raw)
    prepared = prepare_token_view(structure, catalog, 1)
    rows = prepared.source_atom_site_index
    components = np.full(prepared.data.n_components, "?", dtype="<U16")
    for atom, component in enumerate(prepared.data.atom_component):
        if components[component] == "?":
            components[component] = str(structure.atoms.label_comp_id[rows[atom]])
    provenance = dict(
        component_ids=components,
        chain_instance_id=structure.atoms.label_asym_id[rows],
        **ap_novo_pdb.token_fields(structure, prepared),
    )
    view = StoredPdbView(
        str(structure.entry_id), 1, "explicit", prepared.data, provenance
    )
    mapping = provenance["ap_token_index"]
    mapped = np.flatnonzero(mapping >= 0)
    _, offsets = np.unique(mapping[mapped], return_index=True)
    keys = []
    for atom in mapped[offsets]:
        row = rows[atom]
        protein = prepared.data.protein_atom_mask[atom]
        sequence = (
            structure.atoms.label_seq_id[row]
            if protein
            else structure.atoms.auth_seq_id[row]
        )
        keys.append(
            [
                str(structure.atoms.label_asym_id[row]),
                str(sequence),
                str(structure.atoms.pdbx_PDB_ins_code[row]),
                str(structure.atoms.label_comp_id[row]),
                "protein" if protein else str(structure.atoms.label_atom_id[row]),
            ]
        )
    return tokenize(view, hydrogen_complete=complete), keys


@pytest.mark.parametrize("complete", [False, True])
@pytest.mark.parametrize("kind", ["ligand-altloc", "metal"])
def test_selected_pdb_targets_match_python_tokens(kind, complete):
    from kaveh_data.native_tokens import document_tokens

    if kind == "metal":
        raw, catalog = (
            fixtures._coordination_cif().encode(),
            fixtures._coordination_catalog(),
        )
    else:
        raw, catalog = fixtures._cif().encode(), fixtures._catalog()
    document = nibbler.mmcif.select_conformer(nibbler.cif.read(raw)).document
    expected, expected_keys = legacy_tokens(
        document.to_canonical().encode(), catalog, complete
    )
    actual, keys = arrays(
        packer(catalog).pack(
            document,
            complete=complete,
            validate=True,
            **product_options(document, catalog),
        )
    )
    assert keys == expected_keys
    assert set(actual) == set(expected)
    for name in actual:
        np.testing.assert_array_equal(actual[name], expected[name], err_msg=name)
    wrapped, wrapped_keys = document_tokens(
        document, catalog, hydrogen_complete=complete, with_keys=True
    )
    assert wrapped_keys == expected_keys
    for name in wrapped:
        np.testing.assert_array_equal(wrapped[name], expected[name], err_msg=name)


def state_observation(aa, state, *, drop=(), overrides=None):
    residue = schema.SPEC["residues"][aa]
    xyz = np.zeros((32, 3), dtype=np.float32)
    xyz[:15, 0] = np.arange(15) * 4
    counts = {}
    for slot, parent in enumerate(residue["parents"]):
        if parent >= 0:
            index = counts.get(parent, 0)
            counts[parent] = index + 1
            xyz[slot] = xyz[parent] + np.eye(3)[index]
    mask = np.asarray(state["present"])
    names = np.asarray(residue["names"])[mask]
    xyz = xyz[mask]
    for name, position in (overrides or {}).items():
        xyz[names == name] = position
    keep = ~np.isin(names, drop)
    return names[keep], np.asarray([name[0] for name in names[keep]]), xyz[keep]


def state_document(aa, state, *, declared_termini=None, **changes):
    residue = schema.SPEC["residues"][aa]["residue"]
    n_terminal, c_terminal = (
        declared_termini
        if declared_termini is not None
        else (state["n_terminal"], state["c_terminal"])
    )
    index = 1 if n_terminal else 2
    length = index if c_terminal else index + 1
    names, elements, xyz = state_observation(aa, state, **changes)
    header = (
        "data_state\nloop_\n_entity_poly.entity_id\n_entity_poly.type\n"
        "1 'polypeptide(L)'\n"
        f"loop_\n_entity_poly_seq.entity_id\n_entity_poly_seq.num\n1 {length}\n"
        "loop_\n_atom_site.label_entity_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n"
        "_atom_site.auth_seq_id\n_atom_site.label_comp_id\n_atom_site.label_atom_id\n"
        "_atom_site.type_symbol\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n"
        "_atom_site.pdbx_PDB_model_num\n"
    )
    rows = "".join(
        f"1 A {index} {index} {residue} {name} {element} {x} {y} {z} 1\n"
        for name, element, (x, y, z) in zip(names, elements, xyz, strict=True)
    )
    return nibbler.cif.read((header + rows).encode())


@pytest.mark.parametrize("complete", [False, True])
def test_every_catalogue_state_matches_protein_packing(complete):
    native = packer()
    for aa, residue in enumerate(schema.SPEC["residues"]):
        for state in residue["states"]:
            termini = (state["n_terminal"], state["c_terminal"])
            expected = atom32.pack_residue(
                aa, *state_observation(aa, state), complete=complete, termini=termini
            )
            actual, _ = arrays(
                native.pack(
                    state_document(aa, state),
                    complete=complete,
                    topology_policy="entity_sequence",
                    protein_only=True,
                    validate=True,
                )
            )
            for name, value in expected.items():
                np.testing.assert_array_equal(
                    actual[name][0],
                    value,
                    err_msg=f"{residue['residue']} {state['name']} {name}",
                )


@pytest.mark.parametrize(
    "kind", ["missing-terminal", "incomplete-sidechain", "hydrogen-boundary"]
)
def test_partial_protein_evidence_matches_python(kind):
    aa = 8 if kind != "hydrogen-boundary" else 0
    wanted = "delta:N0:C1" if aa == 8 else "standard:N0:C0"
    state = next(
        state
        for state in schema.SPEC["residues"][aa]["states"]
        if state["name"] == wanted
    )
    changes = (
        dict(drop=["OXT"])
        if kind == "missing-terminal"
        else dict(drop=["NE2"])
        if kind == "incomplete-sidechain"
        else dict(overrides={"H": [0.5, 0.0, 0.0]})
    )
    expected = atom32.pack_residue(
        aa,
        *state_observation(aa, state, **changes),
        complete=True,
        termini=(state["n_terminal"], state["c_terminal"]),
    )
    actual, _ = arrays(
        packer().pack(
            state_document(aa, state, **changes),
            complete=True,
            topology_policy="entity_sequence",
            protein_only=True,
            validate=True,
        )
    )
    for name in expected:
        np.testing.assert_array_equal(actual[name][0], expected[name], err_msg=name)
    assert not actual["state_known"][0]


def add_atom_rows(raw, rows):
    marker = "#\nloop_\n_struct_conn.id"
    assert marker in raw
    return raw.replace(marker, rows + marker)


@pytest.mark.parametrize("ambiguous", [False, True])
def test_ligand_hydrogen_coverage_matches_python(ambiguous):
    xyz = "0.7 1.4 0" if ambiguous else "1.4 0 0.97"
    raw = add_atom_rows(
        fixtures._cif(), f"HETATM 15 H H1 . LIG B 2 1 H1 LIG B 1 ? {xyz} 1 10 ? 1\n"
    )
    document = nibbler.mmcif.select_conformer(nibbler.cif.read(raw.encode())).document
    expected, keys = legacy_tokens(
        document.to_canonical().encode(), fixtures._catalog(), True
    )
    actual, actual_keys = arrays(
        packer(fixtures._catalog()).pack(
            document,
            complete=True,
            validate=True,
            **product_options(document, fixtures._catalog()),
        )
    )
    assert actual_keys == keys
    for name in actual:
        np.testing.assert_array_equal(actual[name], expected[name], err_msg=name)
    ligand = actual["is_ligand"] & (actual["chain_idx"] == 1)
    assert not actual["ligand_hydrogen_known"][ligand].any()
    assert actual["mask"][ligand, 15:].sum() == (not ambiguous)


@pytest.mark.parametrize("complete_protein", [False, True])
def test_complete_filter_and_identity_alignment_remap_attachments(complete_protein):
    raw = fixtures._cif()
    if complete_protein:
        raw = add_atom_rows(
            raw,
            "ATOM 15 C CA . ALA A 1 1 CA ALA A 1 ? 1.45 0 0 1 10 ? 1\n"
            "ATOM 16 O O . ALA A 1 1 O ALA A 1 ? 10 1.2 0 1 10 ? 1\n"
            "ATOM 17 C CB . ALA A 1 1 CB ALA A 1 ? 1.45 1.5 0 1 10 ? 1\n",
        )
    catalog = fixtures._catalog()
    document = nibbler.mmcif.select_conformer(nibbler.cif.read(raw.encode())).document
    expected, keys = legacy_tokens(document.to_canonical().encode(), catalog, True)
    selected = np.flatnonzero(expected["is_ligand"] | complete_protein)[::-1]
    wanted_keys = [keys[row] for row in selected]
    inverse = np.full(len(keys), -1, dtype=np.int64)
    inverse[selected] = np.arange(len(selected))
    edge = expected["protein_attachment_index"]
    edge = inverse[edge[:, (inverse[edge] >= 0).all(0)]]
    expected = {
        name: value[selected]
        for name, value in expected.items()
        if name != "protein_attachment_index"
    }
    expected["protein_attachment_index"] = np.unique(edge, axis=1)
    expected["covalently_bound_to_protein"][:] = False
    expected["covalently_bound_to_protein"][edge[0]] = True
    for chain in np.unique(expected["chain_idx"][expected["is_ligand"]]):
        rows = np.flatnonzero(expected["is_ligand"] & (expected["chain_idx"] == chain))
        expected["res_idx"][rows] = np.arange(len(rows))
    actual, actual_keys = arrays(
        packer(catalog).pack(
            document,
            complete=True,
            complete_tokens=True,
            expected_keys=json.dumps(wanted_keys),
            validate=True,
            **product_options(document, catalog),
        )
    )
    assert actual_keys == wanted_keys
    for name in actual:
        np.testing.assert_array_equal(actual[name], expected[name], err_msg=name)


def test_modified_lysine_clears_exact_state_and_absent_hydrogen_labels():
    aa = next(
        i
        for i, residue in enumerate(schema.SPEC["residues"])
        if residue["residue"] == "LYS"
    )
    state = schema.SPEC["residues"][aa]["states"][0]
    document = state_document(aa, state)
    native = packer()
    before, keys = arrays(
        native.pack(
            document,
            complete=True,
            topology_policy="entity_sequence",
            protein_only=True,
        )
    )
    assert before["state_known"].tolist() == [True]
    target = dict(
        parent_key=[*keys[0][:4], "NZ"],
        sequence_id=int(keys[0][1]),
        formal_charge=0,
        hybridization=4,
        hydrogens=1,
        complete=False,
        modified=True,
    )
    after, _ = arrays(
        native.pack(
            document,
            complete=True,
            topology_policy="entity_sequence",
            protein_only=True,
            product_targets=json.dumps([target]),
        )
    )
    assert not after["state_known"].any()
    np.testing.assert_array_equal(
        after["atom_presence_known"][:, 15:], after["mask"][:, 15:]
    )


@pytest.mark.parametrize("with_keys", [False, True])
def test_incompatible_hydrogens_report_row_and_token_context(with_keys):
    state = next(
        state
        for state in schema.SPEC["residues"][0]["states"]
        if state["name"] == "standard:N2:C0"
    )
    document = state_document(0, state, declared_termini=(False, False))
    native = packer()
    unchecked, _ = arrays(
        native.pack(
            document,
            complete=True,
            protein_only=True,
            topology_policy="entity_sequence",
        )
    )
    assert unchecked["topology"][0].tolist() == [0, 0, -1]
    with pytest.raises(
        ValueError, match="Atom-presence labels contradict every allowed chemical state"
    ) as error:
        native.pack(
            document,
            complete=True,
            protein_only=True,
            topology_policy="entity_sequence",
            validate=True,
            with_keys=with_keys,
        )
    context = json.loads(str(error.value).partition(": ")[2])
    assert context["row"] == 0 and context["aatype"] == 0
    assert context["token_key"] == ["A", "2", ".", "ALA", "protein"]
    assert context["topology"] == [0, 0, -1]
    assert context["mask"] == unchecked["mask"][0].astype(int).tolist()
    assert (
        context["atom_presence_known"]
        == unchecked["atom_presence_known"][0].astype(int).tolist()
    )


def test_suppressing_returned_keys_keeps_arrays_identical():
    document = nibbler.mmcif.select_conformer(
        nibbler.cif.read(fixtures._cif().encode())
    ).document
    native = packer(fixtures._catalog())
    expected, keys = arrays(native.pack(document, complete=True, validate=True))
    actual, no_keys = arrays(
        native.pack(
            document,
            complete=True,
            validate=True,
            with_keys=False,
            expected_keys=json.dumps(keys),
        )
    )
    assert keys and no_keys == []
    assert set(actual) == set(expected)
    for name in actual:
        assert actual[name].dtype == expected[name].dtype
        np.testing.assert_array_equal(actual[name], expected[name], err_msg=name)


def reference_afdb_fixture(pos, mask):
    """Match the pre-native AFDB constructor for this three-residue fixture."""
    fields = (
        "group_PDB id type_symbol label_atom_id label_alt_id label_comp_id "
        "label_asym_id label_entity_id label_seq_id auth_atom_id auth_comp_id "
        "auth_asym_id auth_seq_id pdbx_PDB_ins_code Cartn_x Cartn_y Cartn_z occupancy "
        "B_iso_or_equiv pdbx_PDB_model_num"
    ).split()
    header = (
        "data_AF-test\n_entry.id AF-test\nloop_\n_entity.id\n_entity.type\n"
        "1 polymer\n2 polymer\n"
        "loop_\n_entity_poly.entity_id\n_entity_poly.type\n"
        "_entity_poly.pdbx_seq_one_letter_code_can\n"
        "1 'polypeptide(L)' AA\n2 'polypeptide(L)' A\n"
        "loop_\n_struct_asym.id\n_struct_asym.entity_id\nA1 1\nA2 2\nloop_\n"
        + "".join(f"_atom_site.{field}\n" for field in fields)
    )
    rows = []
    for row, (chain, sequence) in enumerate(((1, 1), (1, 2), (2, 1))):
        for slot in np.flatnonzero(mask[row]):
            name = schema.SPEC["residues"][0]["names"][slot]
            x, y, z = map(float, pos[row, slot])
            rows.append(
                f"ATOM {len(rows) + 1} {name[0]} {name} . "
                f"ALA A{chain} {chain} {sequence} "
                f"{name} ALA A{chain} {sequence} . {x} {y} {z} 1 0 1\n"
            )
    return nibbler.cif.read((header + "".join(rows)).encode())


def test_afdb_document_matches_python_without_text_construction_or_reparse():
    import gemmi
    import pyarrow as pa
    from kaveh_data.atom32_export import protein_arrays_from_bytes
    from kaveh_data.native_tokens import protein_arrays

    pos = np.zeros((3, 14, 3), dtype=np.float32)
    pos[:, :5, 0] = np.arange(3)[:, None] * 3.8 + np.array([0, 1.45, 2.6, 2.6, 1.45])
    pos[2] += 40
    mask = np.zeros((3, 14), dtype=bool)
    mask[:, :5] = True
    aa = np.zeros(3, dtype=np.int64)
    expected = reference_afdb_fixture(pos, mask)
    native = packer()
    with (
        patch.object(
            nibbler.cif,
            "read",
            side_effect=AssertionError("native construction reparsed text"),
        ),
        patch.object(
            gemmi.cif,
            "Document",
            side_effect=AssertionError("native construction used Gemmi"),
        ),
    ):
        actual, actual_chains, actual_sequence = native.protein_document(
            pos.tobytes(), mask.tobytes(), aa.tobytes(), "AF-test"
        )
        targets, _ = arrays(native.pack(actual, protein_only=True, validate=True))
    assert (actual_chains, actual_sequence) == (2, "AAA")
    for category in ("entry", "entity", "entity_poly", "struct_asym", "atom_site"):
        projected = [
            pa.RecordBatchReader.from_stream(
                nibbler.cif.project(document, category, schema="pdbx")
            ).read_all()
            for document in (expected, actual)
        ]
        assert projected[0].equals(projected[1]), category
    expected_targets, _ = protein_arrays_from_bytes(
        expected.to_canonical().encode(), topology_policy="observed_connectivity"
    )
    for name in targets:
        np.testing.assert_array_equal(
            targets[name], expected_targets[name], err_msg=name
        )
    wrapped = protein_arrays(actual, topology_policy="observed_connectivity")
    for name in wrapped:
        np.testing.assert_array_equal(
            wrapped[name], expected_targets[name], err_msg=name
        )
    np.testing.assert_array_equal(
        targets["pos"][:, :14], np.where(mask[..., None], pos, 0)
    )
    np.testing.assert_array_equal(targets["mask"][:, :14], mask)
    np.testing.assert_array_equal(
        targets["topology"][:, :2], [[-1, 0], [0, -1], [-1, -1]]
    )
    np.testing.assert_array_equal(targets["res_idx"], [1, 2, 1])
    np.testing.assert_array_equal(targets["chain_idx"], [0, 0, 1])
