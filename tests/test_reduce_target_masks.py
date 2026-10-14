"""Typed target masking on saved views; set KAVEH_TARGET_MASK_RUN for this gate."""

import ast
import json
import os
from hashlib import sha256
from pathlib import Path

import gemmi
import numpy as np
import pyarrow as pa
import pytest

import nibbler

CASES = {
    "5LBM": (
        "owner-07/train-000039-0053",
        "2387796b0aaf75813b3dfa456597dd6e0ab97a919c3afe08f59f160afd8e3729",
    ),
    "9L0K": (
        "owner-08/train-000072-0054",
        "6e9d2d9219069fef5b39c58c578426f69f6992089dedb023f4f8d6a7625677f8",
    ),
    "2H59": (
        "owner-04/train-000676-0148",
        "6c06a5eff5b24a0a6369ddfbcef688ec6de937cd073a37354e91c5017bdf79ae",
    ),
    "1NLM": (
        "owner-03/train-000963-0188",
        "69c356026a4f625d17223b8025f570697274d4430f2f7f5b657618642325fe99",
    ),
}
ADDUCT = Path(__file__).parent / "reduce3/linked_9hj3.cif"


def rows(document):
    return (
        pa.RecordBatchReader.from_stream(document.project("atom_site"))
        .read_all()
        .to_pylist()
    )


@pytest.fixture(scope="module")
def gate():
    configured = os.environ.get("KAVEH_TARGET_MASK_RUN")
    if configured is None:
        pytest.skip("Set KAVEH_TARGET_MASK_RUN for the saved target-mask gate")
    from kaveh_data.ap_novo_pdb import EXCLUDED_COMPONENTS
    from kaveh_data.atom32_schema import ENCODED
    from kaveh_data.ccd_catalog import CCDCatalog
    from kaveh_data.native_observations import NativeObservations

    run = Path(configured)
    assets = run / "cluster-preparation"
    source = run.parents[1] / "kaveh_novo/data/records.py"
    function = next(
        node
        for node in ast.parse(source.read_text()).body
        if isinstance(node, ast.FunctionDef) and node.name == "decode_record"
    )
    namespace = {"np": np}
    exec(
        compile(ast.Module(body=[function], type_ignores=[]), str(source), "exec"),
        namespace,
    )
    native = NativeObservations(
        CCDCatalog.open(assets / "ccd_catalog_v2"), chem_data=assets / "chem_data"
    )
    return (
        run,
        assets / "chem_data",
        namespace["decode_record"],
        native,
        json.loads(ENCODED),
        sorted(EXCLUDED_COMPONENTS),
    )


@pytest.fixture(scope="module")
def saved(gate):
    run, _, decode, _, _, _ = gate
    result = {}
    for entry, (relative, digest) in CASES.items():
        descriptor = run / "full-pdb-audit-v1-failures" / f"{relative}.json"
        metadata = json.loads(descriptor.read_text())
        frame = descriptor.with_suffix(".zstd").read_bytes()
        assert (
            sha256(frame).hexdigest()
            == digest
            == metadata["frame"]["sha256"]
            == metadata["record_index"]["sha256"]
        )
        assert (
            len(frame) == metadata["frame"]["bytes"] == metadata["record_index"]["size"]
        )
        raw = pa.decompress(
            frame, decompressed_size=metadata["record_index"]["raw_size"], codec="zstd"
        )
        assert len(raw) == metadata["record_index"]["raw_size"]
        record = decode(raw, metadata["fields"])
        model = int(record["model"][0])
        assert (
            entry == metadata["provenance"]["entry_id"]
            and model == metadata["provenance"]["model_id"]
        )
        result[entry] = (
            nibbler.cif.read(record["cif"].tobytes()),
            model,
            record["token_keys"].tobytes().decode(),
        )
    return result


def reduce(gate, source, model, keys, **options):
    return nibbler.reduce.run(
        source,
        chem_data=gate[1],
        model_id=model,
        required_hydrogen_target_keys=keys,
        stop_on_any_missing_hydrogen=True,
        add_flip_movers=True,
        allow_missing_restraints_for=gate[5],
        **options,
    )


def packed(gate, output, model, keys, statuses=(), validate=False):
    return gate[3].pack(
        output,
        model=model,
        complete=True,
        complete_tokens=True,
        expected_keys=keys,
        validate=validate,
        target_uncertainties=statuses,
    )


def expected_protein_masks(expected, before, row, parent_name, status, schema):
    residue = schema["residues"][int(before["aatype"][row])]
    parent = residue["names"].index(parent_name)
    if status["chemical_unknown"]:
        expected["state_known"][row] = False
        flag = {"N": 0, "C": 1, "OXT": 1, "SG": 2}.get(parent_name)
        if flag is not None:
            expected["topology"][row, flag] = -1
    for slot in range(15, 32):
        if residue["parents"][slot] == parent and not before["mask"][row, slot]:
            expected["atom_presence_known"][row, slot] = False


def check_selective_masks(before, after, keys, statuses, schema):
    expected = {name: value.copy() for name, value in before.items()}
    lookup = {tuple(key): index for index, key in enumerate(keys)}
    for status in statuses:
        key = tuple(status["parent_key"])
        row = lookup.get(key, lookup.get((*key[:4], "protein")))
        assert row is not None
        if before["is_ligand"][row]:
            if status["chemical_unknown"]:
                expected["ligand_charge_known"][row] = False
            if status["hydrogens_unknown"]:
                expected["ligand_hydrogen_known"][row] = False
        else:
            expected_protein_masks(expected, before, row, key[4], status, schema)
    assert after.keys() == expected.keys()
    for name in after:
        np.testing.assert_array_equal(after[name], expected[name], err_msg=name)


def source_keys(document, retained):
    wanted = {tuple(key) for key in json.loads(retained)}
    result = {}
    for row, atom in enumerate(rows(document)):
        if atom["type_symbol"] in {"H", "D", "T"}:
            continue
        prefix = (
            atom["label_asym_id"],
            str(atom["label_seq_id"]),
            atom["pdbx_PDB_ins_code"] or "?",
            atom["label_comp_id"],
        )
        sequence = (
            prefix[1] if (*prefix, "protein") in wanted else str(atom["auth_seq_id"])
        )
        key = (prefix[0], sequence, prefix[2], prefix[3], atom["label_atom_id"])
        assert key not in result
        result[key] = row
    return result


def check_source_positions(gate, source, after, keys, affected, expected):
    atoms = rows(source)
    for row in expected:
        atom = atoms[row]
        key = affected[row]["parent_key"]
        token = keys.index(key if key in keys else [*key[:4], "protein"])
        slot = (
            0
            if after["is_ligand"][token]
            else gate[4]["residues"][int(after["aatype"][token])]["names"].index(key[4])
        )
        np.testing.assert_array_equal(
            after["pos"][token, slot],
            np.asarray([atom[f"Cartn_{axis}"] for axis in "xyz"], dtype=np.float32),
        )
        assert after["mask"][token, slot]


@pytest.mark.parametrize(
    "flag", ["allow_undefined_ligand_targets", "allow_unresolved_bond_targets"]
)
@pytest.mark.parametrize("value", [1, "yes", None])
def test_target_flags_require_booleans(flag, value):
    with pytest.raises(TypeError):
        nibbler.reduce.run(ADDUCT, **{flag: value})


@pytest.mark.parametrize(
    "flag", ["allow_undefined_ligand_targets", "allow_unresolved_bond_targets"]
)
@pytest.mark.parametrize(
    "bad",
    [
        {"required_hydrogen_target_keys": None},
        {"model_id": None},
        {"stop_on_any_missing_hydrogen": False},
        {"compat": True},
        {"approach": "remove"},
    ],
)
def test_target_flags_require_strict_fixed_addition_and_retained_keys(flag, bad):
    options = dict(
        model_id=1,
        required_hydrogen_target_keys="[]",
        stop_on_any_missing_hydrogen=True,
    )
    options.update(bad)
    exception = (
        ValueError
        if bad.get("required_hydrogen_target_keys", "[]") is None
        else nibbler.ChemistryError
    )
    with pytest.raises(exception) as error:
        nibbler.reduce.run(ADDUCT, **options, **{flag: True})
    if exception is nibbler.ChemistryError:
        assert error.value.code == "REDUCE_INVALID_MODEL"


@pytest.mark.parametrize(
    "entry,expected,reason",
    [
        ("5LBM", {2744, 2745}, "insufficient_restraints"),
        ("9L0K", {45389}, "insufficient_restraints"),
        ("2H59", {983, 1175, 1197}, "multiple_sulfur_partners"),
        ("1NLM", {4063, 4064, 4157, 4158}, "degenerate_heavy_bond"),
    ],
)
def test_exact_source_rows_reach_only_their_stable_native_target_keys(
    gate, saved, entry, expected, reason
):
    source, model, retained = saved[entry]
    original = source.to_canonical()
    output = reduce(
        gate,
        source,
        model,
        retained,
        allow_undefined_ligand_targets=True,
        allow_unresolved_bond_targets=True,
    )
    assert source.to_canonical() == original
    mapping = source_keys(source, retained)
    required = set(
        nibbler.reduce.required_parent_sources(
            source, retained, model_id=model, chem_data=gate[1]
        )
    )
    statuses = output.target_uncertainties
    affected = {mapping[tuple(status["parent_key"])]: status for status in statuses}
    assert expected <= affected.keys() <= required
    if entry in {"2H59", "1NLM"}:
        assert affected.keys() == expected
    for row in expected:
        assert affected[row]["chemical_unknown"] and affected[row]["hydrogens_unknown"]
        assert reason in affected[row]["reasons"]
    before, keys = packed(gate, output.document, model, retained)
    after, same_keys = packed(
        gate, output.document, model, retained, statuses, validate=True
    )
    assert keys == same_keys == json.loads(retained)
    check_selective_masks(before, after, keys, statuses, gate[4])
    check_source_positions(gate, source, after, keys, affected, expected)


@pytest.mark.parametrize("entry", ["5LBM", "9L0K", "1NLM"])
def test_exact_failures_remain_strict_without_target_allowances(gate, saved, entry):
    with pytest.raises(nibbler.ChemistryError) as error:
        reduce(gate, *saved[entry])
    assert error.value.code == "REDUCE_FAILED"


def observed_hg(tags, parent):
    h = list(parent)
    values = {
        "id": "999001",
        "label_atom_id": "HG",
        "auth_atom_id": "HG",
        "type_symbol": "H",
        "occupancy": "0.43",
        "B_iso_or_equiv": "27.1",
        "pdbx_formal_charge": "-1",
    }
    for axis, offset in zip("xyz", (0.27, 0.31, 0.87), strict=True):
        tag = f"Cartn_{axis}"
        values[tag] = str(float(h[tags.index(f"_atom_site.{tag}")]) + offset)
    for tag, value in values.items():
        h[tags.index(f"_atom_site.{tag}")] = value
    return h


def changed_adduct(observed=False, missing_cb=False, appended=False):
    doc = gemmi.cif.read_file(str(ADDUCT))
    table = doc.sole_block().find_mmcif_category("_atom_site.")
    tags, atoms = list(table.tags), [list(row) for row in table]
    name, comp = (
        tags.index("_atom_site.label_atom_id"),
        tags.index("_atom_site.label_comp_id"),
    )
    if observed:
        parent = next(
            index
            for index, row in enumerate(atoms)
            if row[comp] == "CYS" and row[name] == "SG"
        )
        atoms.insert(
            len(atoms) if appended else parent + 1, observed_hg(tags, atoms[parent])
        )
    if missing_cb:
        atoms = [row for row in atoms if row[comp] == "CYS" and row[name] != "CB"]
        doc.sole_block().find_mmcif_category("_struct_conn.").erase()
    table.erase()
    loop = doc.sole_block().init_loop("", tags)
    for atom in atoms:
        loop.add_row(atom)
    return nibbler.cif.read(doc.as_string().encode())


@pytest.mark.parametrize("keep", [False, True])
@pytest.mark.parametrize("appended", [False, True], ids=["contiguous", "unsorted"])
def test_observed_hg_keeps_coordinates_presence_and_metadata(gate, keep, appended):
    source = changed_adduct(observed=True, appended=appended)
    retained = json.dumps([["A", "1", "?", "CYS", "protein"]])
    original = next(atom for atom in rows(source) if atom["label_atom_id"] == "HG")
    output = reduce(
        gate,
        source,
        1,
        retained,
        keep_existing_h=keep,
        allow_unresolved_bond_targets=True,
    )
    assert output.target_uncertainties == (
        {
            "parent_key": ["A", "1", "?", "CYS", "SG"],
            "hydrogens_unknown": False,
            "chemical_unknown": True,
            "reasons": ["unsupported_residue_environment"],
        },
    )
    actual = next(
        atom for atom in rows(output.document) if atom["label_atom_id"] == "HG"
    )
    assert {name: value for name, value in actual.items() if name != "id"} == {
        name: value for name, value in original.items() if name != "id"
    }
    before, keys = packed(gate, output.document, 1, None)
    after, same_keys = packed(
        gate, output.document, 1, None, output.target_uncertainties, validate=True
    )
    assert same_keys == keys
    check_selective_masks(before, after, keys, output.target_uncertainties, gate[4])
    slot = gate[4]["residues"][4]["names"].index("HG")
    row = keys.index(json.loads(retained)[0])
    assert after["mask"][row, slot] and after["atom_presence_known"][row, slot]


def test_ligand_and_bond_allowances_do_not_permit_ordinary_protein_missing_h(gate):
    with pytest.raises(nibbler.ChemistryError) as error:
        reduce(
            gate,
            changed_adduct(missing_cb=True),
            1,
            json.dumps([["A", "1", "?", "CYS", "protein"]]),
            allow_undefined_ligand_targets=True,
            allow_unresolved_bond_targets=True,
        )
    assert error.value.code == "REDUCE_FAILED"


def bridge_source():
    return nibbler.cif.read(b"""data_bridge
loop_
_atom_site.group_PDB
_atom_site.label_atom_id
_atom_site.type_symbol
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_PDB_model_num
HETATM FE FE . EJ2 L L 601 ? -4.600 -5.078 26.969 1 20 1
HETATM NI NI . EJ2 L L 601 ? -5.053 -2.663 26.213 1 20 1
HETATM C1 C . EJ2 L L 601 ? -3.495 -5.027 28.482 1 20 1
HETATM N1 N . EJ2 L L 601 ? -2.796 -4.998 29.405 1 20 1
HETATM C2 C . EJ2 L L 601 ? -3.227 -5.893 25.921 1 20 1
HETATM N2 N . EJ2 L L 601 ? -2.353 -6.357 25.319 1 20 1
HETATM C3 C . EJ2 L L 601 ? -5.225 -6.686 27.290 1 20 1
HETATM O3 O . EJ2 L L 601 ? -5.613 -7.733 27.462 1 20 1
""")


def test_missing_bridge_coordinate_keeps_h65_chemistry_and_both_bonds(gate):
    source = bridge_source()
    original, keys = gate[3].pack(source, complete=True, validate=True)
    retained = json.dumps(keys)
    output = reduce(gate, source, 1, retained, allow_undefined_ligand_targets=True)
    assert output.target_uncertainties == ()
    after, same_keys = packed(
        gate, output.document, 1, retained, output.target_uncertainties, validate=True
    )
    assert same_keys == keys
    h, fe, ni = (
        keys.index(["L", "601", "?", "EJ2", atom]) for atom in ("H", "FE", "NI")
    )
    assert after["aatype"][h] == 65 and after["ligand_charge"][h] == -1
    assert after["ligand_charge_known"][h]
    np.testing.assert_array_equal(
        after["ligand_charge_known"][[fe, ni]],
        original["ligand_charge_known"][[fe, ni]],
    )
    assert after["ligand_hydrogen_known"][[h, fe, ni]].all()
    assert not after["pos"][h].any() and not after["mask"][h].any()
    assert not after["mask"][[h, fe, ni], 15:].any()
    np.testing.assert_array_equal(after["ligand_bridge_index"], [[h, h], [fe, ni]])
    np.testing.assert_array_equal(after["ligand_bridge_type"], [1, 1])
    np.testing.assert_array_equal(
        after["crop_positions"][h], after["pos"][[fe, ni], 0].mean(axis=0)
    )
    assert not after["crop_center_eligible"][h]
