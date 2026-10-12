"""Check reference serialization and real CCTBX capture without skipped checks."""

import gzip
import importlib.util
import json
import os
from pathlib import Path

import pytest

MODULE = Path(__file__).parents[1] / "benchmarks/native_validation_reference_capture.py"
SPEC = importlib.util.spec_from_file_location("validation_reference", MODULE)
REFERENCE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(REFERENCE)


@pytest.fixture(autouse=True)
def standalone_benchmark_modules(monkeypatch):
    monkeypatch.syspath_prepend(str(MODULE.parent))


def test_proxy_rows_preserve_weights_and_origins():
    class Proxy:
        i_seqs = (3, 1)
        weight = 125.5
        origin_id = 7
        slack = 0.25

    assert REFERENCE.proxy_rows(
        [Proxy()], ("i_seqs", "weight", "origin_id", "slack")
    ) == [dict(i_seqs=[3, 1], weight=125.5, origin_id=7, slack=0.25)]


def test_proxy_rows_preserve_plane_weight_vector():
    class Proxy:
        i_seqs = (4, 8, 1, 0)
        weights = (1.0, 2.0, 3.0, 4.0)

    assert REFERENCE.proxy_rows([Proxy()], ("i_seqs", "weights")) == [
        dict(i_seqs=[4, 8, 1, 0], weights=[1.0, 2.0, 3.0, 4.0])
    ]


@pytest.fixture
def reference_chem_data():
    configured = os.environ.get("NIBBLER_REFERENCE_CHEM_DATA")
    if not configured:
        pytest.skip("Set NIBBLER_REFERENCE_CHEM_DATA for live reference capture")
    return Path(configured)


def test_capture_declares_unknown_components(tmp_path, reference_chem_data):
    pytest.importorskip("mmtbx.model", reason="Live reference capture needs CCTBX")
    fixture = tmp_path / "unknown.pdb"
    fixture.write_text(
        "HETATM    1  C1  ZZZ A   1       0.000   0.000   0.000"
        "  1.00 10.00           C\nEND\n"
    )
    report = REFERENCE.capture(fixture, reference_chem_data)
    assert report["complete"] is False
    assert report["error"]
    assert report["chemistry_sources"]["list/mon_lib_list.cif"]["sha256"]


@pytest.mark.parametrize(
    "name",
    [
        "valid_alanine",
        "valid_d_alanine",
        "metal_oxygen_contact",
        "overlapping_residues",
        "displaced_hydrogen",
        "microheterogeneity",
        "crn5-electron",
    ],
)
def test_frozen_probe_reference_covers_processed_graph(name):
    folder = Path(__file__).with_name("native_validation_reference_cases")
    with gzip.open(folder / f"{name}-reference.json.gz", "rt") as stream:
        report = json.load(stream)
    assert report["complete"]
    assert report["graph"]["bonds"]
    assert report["chemistry_sources"]["reduce_legacy_list"]["sha256"]
    initial = report["probe2"]["initial_atom_info_calls"][0]
    assert len(initial["atoms"]) == len(initial["extra_atom_info"])
    assert [row["index"] for row in initial["extra_atom_info"]] == list(
        range(len(initial["atoms"]))
    )
    assert all(row["vdw_radius"] > 0 for row in initial["extra_atom_info"])
    assert all(
        0 <= i < len(initial["atoms"]) for pair in initial["bonds"] for i in pair
    )
    assert report["probe2"]["program_runs"][0]["effective_phil"]


def test_independent_graph_comparison_maps_atom_order_and_detects_target_error():
    import copy

    spec = importlib.util.spec_from_file_location(
        "compare_reference", MODULE.with_name("native_validation_reference_compare.py")
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    atoms = [
        dict(
            model="",
            chain="A",
            resseq="1",
            icode="",
            component="ALA",
            altloc="",
            name=name,
            element="C",
            xyz=[0, 0, i],
        )
        for i, name in enumerate(("CA", "CB"))
    ]
    reference = dict(
        complete=True,
        atoms=atoms,
        graph=dict(
            bonds=[
                dict(
                    i_seqs=[0, 1], distance_ideal=1.5, weight=100, slack=0, origin_id=0
                )
            ]
        ),
    )
    native = dict(
        complete=True,
        atoms=list(reversed(atoms)),
        graph=dict(bonds=[[1, 0, 1.5, 100, 0, True]]),
    )
    assert module.compare(reference, native)["graph_parity"]
    wrong = copy.deepcopy(native)
    wrong["graph"]["bonds"][0][2] = 1.6
    assert not module.compare(reference, wrong)["graph_parity"]


def test_incomplete_native_contacts_do_not_claim_required_parity():
    spec = importlib.util.spec_from_file_location(
        "compare_reference", MODULE.with_name("native_validation_reference_compare.py")
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    folder = Path(__file__).with_name("native_validation_reference_cases")
    with gzip.open(folder / "hydrogen_bond_contact-reference.json.gz", "rt") as stream:
        reference = json.load(stream)
    native = dict(
        complete=False,
        atoms=reference["atoms"],
        graph=dict(bonds=[], angles=[], chiralities=[], planes=[]),
        contacts=None,
        coverage=[dict(check="hydrogens", unresolved=2, unsupported=["fixture"])],
    )
    result = module.compare(reference, native)
    assert not result["all_required_parity"]
    assert not result["initial_probe_typing_parity"]
    assert result["native_coverage"][0]["unresolved"] == 2


def test_packed_reference_preserves_scalar_results_and_rejects_unknown_atoms(
    monkeypatch,
):
    import copy

    monkeypatch.syspath_prepend(str(MODULE.parent))
    spec = importlib.util.spec_from_file_location(
        "pack_reference", MODULE.with_name("native_validation_reference_pack.py")
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    atom = dict(
        index=0,
        model="",
        chain="A",
        resseq="1",
        icode="",
        component="ALA",
        altloc="",
        name="CA",
        element="C",
    )
    info = dict(
        model_id="",
        chain_id="A",
        resseq="1",
        icode="",
        resname="ALA",
        altloc="",
        name="CA",
        element="C",
        symop=None,
    )
    source = dict(
        atoms=[atom],
        restraint_results=dict(
            bonds=dict(
                rows=[
                    dict(target=1.5, model=1.6, sigma=0.02, score=5, atoms_info=[info])
                ]
            )
        ),
    )
    result = module.pack(source)
    row = result["restraint_results"]["bonds"]["rows"][0]
    assert row == dict(target=1.5, model=1.6, sigma=0.02, score=5, atom_indices=[0])
    assert "atoms_info" in source["restraint_results"]["bonds"]["rows"][0]
    wrong = copy.deepcopy(source)
    wrong["restraint_results"]["bonds"]["rows"][0]["atoms_info"][0]["name"] = "CB"
    with pytest.raises(ValueError, match="lacks a processed atom-table match"):
        module.pack(wrong)


def test_missing_chirality_policy_flag_is_not_full_graph_parity():
    spec = importlib.util.spec_from_file_location(
        "compare_reference", MODULE.with_name("native_validation_reference_compare.py")
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    atoms = [
        dict(
            model="",
            chain="A",
            resseq="1",
            icode="",
            component="VAL",
            altloc="",
            name=name,
            element="C",
            xyz=[0, 0, i],
        )
        for i, name in enumerate(("CB", "CA", "CG1", "CG2"))
    ]
    proxy = dict(i_seqs=[0, 1, 2, 3], volume_ideal=2.5, weight=25, both_signs=False)
    reference = dict(complete=True, atoms=atoms, graph=dict(chiralities=[proxy]))
    native = dict(
        complete=True,
        atoms=atoms,
        graph=dict(chiralities=[[0, 1, 2, 3, 2.5, 25, False, True]]),
    )
    result = module.compare(reference, native)
    assert result["graph_target_weight_parity"]
    assert not result["graph_parity"]
    proxy["pseudochiral"] = True
    assert module.compare(reference, native)["graph_parity"]


@pytest.mark.parametrize(
    "name",
    [
        "valid_alanine",
        "valid_d_alanine",
        "metal_oxygen_contact",
        "hydrogen_bond_contact",
        "wrong_handedness",
        "displaced_hydrogen",
        "overlapping_residues",
        "collapsed_stereocentre",
        "1crn-prepared",
        "protein-2160",
        "crn5-electron",
        "crn5-electron-no-cdl",
        "water_phantom",
    ],
)
def test_live_native_matches_independent_frozen_graph_and_probe(
    name, reference_chem_data
):
    import hashlib

    pytest.importorskip("nibbler.geometry", reason="Live native parity needs Nibbler")
    import nibbler._core
    import nibbler.geometry

    if not hasattr(nibbler._core, "validate_geometry_document"):
        pytest.skip("Live native parity needs a native-validation build")
    folder = Path(__file__).with_name("native_validation_reference_cases")
    with gzip.open(folder / f"{name}-reference.json.gz", "rt") as stream:
        reference = json.load(stream)
    input_name = {
        "protein-2160": "protein-2160-cctbx-sites.cif",
        "crn5-electron-no-cdl": "crn5-electron.cif",
    }.get(name, f"{name}.cif")
    source = folder / input_name
    payload = source.read_bytes()
    assert hashlib.sha256(payload).hexdigest() == reference["input"]["sha256"]
    native = nibbler.geometry.validate(
        payload,
        chem_data=reference_chem_data,
        nuclear=reference["nuclear"],
        use_cdl=reference["use_cdl"],
        include_graph=True,
    )
    native["input"] = dict(sha256=hashlib.sha256(payload).hexdigest())
    spec = importlib.util.spec_from_file_location(
        "compare_reference", MODULE.with_name("native_validation_reference_compare.py")
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    comparison = module.compare(reference, native)
    assert comparison["all_required_parity"], comparison


def test_scalar_comparison_uses_plane_max_score_and_strict_cutoff():
    from native_validation_reference_decisions import compare_decisions

    atoms = [dict(index=i, name=f"C{i}", element="C") for i in range(3)]
    reference = dict(
        atoms=atoms,
        graph={},
        restraint_results={
            "planarities": {
                "rows": [
                    dict(
                        atom_indices=[0, 1, 2],
                        score=6,
                        residual=3,
                        rms_deltas=0.2,
                        delta_max=0.4,
                    )
                ]
            },
            "bonds": {"rows": [dict(atom_indices=[0, 1], score=5)]},
        },
    )
    native = dict(
        atoms=atoms,
        failures=[
            dict(
                check="planarities",
                indices=[2, 0, 1],
                score=6,
                residual=3,
                model=0.2,
                delta=0.4,
            )
        ],
    )
    assert compare_decisions(reference, native)["decision_scalar_parity"]
    native["failures"][0]["model"] = 0.4
    assert not compare_decisions(reference, native)["decision_scalar_parity"]
    del native["failures"][0]["model"]
    assert not compare_decisions(reference, native)["decision_scalar_parity"]


def test_scalar_comparison_excludes_pseudochirality():
    from native_validation_reference_decisions import compare_decisions

    atoms = [dict(index=i, name=f"C{i}", element="C") for i in range(4)]
    reference = dict(
        atoms=atoms,
        graph={"chiralities": [dict(i_seqs=[0, 1, 2, 3], pseudochiral=True)]},
        restraint_results={
            "chiralities": {"rows": [dict(atom_indices=[0, 1, 2, 3], score=10)]}
        },
    )
    assert compare_decisions(reference, dict(atoms=atoms, failures=[]))[
        "decision_scalar_parity"
    ]


def test_scalar_comparison_rejects_duplicate_findings():
    from native_validation_reference_decisions import compare_decisions

    atoms = [dict(index=i, name=f"C{i}", element="C") for i in range(2)]
    row = dict(atom_indices=[0, 1], score=6)
    reference = dict(
        atoms=atoms, graph={}, restraint_results={"bonds": {"rows": [row, row]}}
    )
    with pytest.raises(ValueError, match="duplicate"):
        compare_decisions(reference, dict(atoms=atoms, failures=[]))


def test_live_missing_hydrogens_report_incomplete_coverage(reference_chem_data):
    pytest.importorskip("nibbler.geometry", reason="Live coverage needs Nibbler")
    import nibbler._core
    import nibbler.geometry

    if not hasattr(nibbler._core, "validate_geometry_document"):
        pytest.skip("Live coverage needs a native-validation build")
    folder = Path(__file__).with_name("native_validation_reference_cases")
    report = nibbler.geometry.validate(
        (folder / "microheterogeneity.cif").read_bytes(),
        chem_data=reference_chem_data,
        use_cdl=False,
        include_graph=True,
    )
    assert not report["complete"]
    assert report["checks"]["hydrogens"]["status"] == "incomplete"
    assert report["checks"]["hydrogens"]["unresolved"] > 0
    assert report["checks"]["clashes"]["status"] == "incomplete"
    assert report["contacts"] is None
