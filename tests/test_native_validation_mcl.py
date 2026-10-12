"""Keep independent metal-selection and annotated-image captures."""

import builtins
import gzip
import hashlib
import json
import os
from pathlib import Path

import pytest

from nibbler.geometry import validate

CASES = Path(__file__).with_name("native_validation_mcl")


@pytest.fixture
def chemistry():
    root = os.environ.get("NIBBLER_VALIDATION_CHEM_DATA")
    if not root:
        pytest.skip("Set NIBBLER_VALIDATION_CHEM_DATA to the pinned chemistry tree")
    return root


@pytest.mark.parametrize(
    "case",
    json.loads((CASES / "manifest.json").read_text())["cases"],
    ids=lambda row: row["name"],
)
def test_native_mcl_matches_unmodified_reference(case, chemistry, monkeypatch):
    monkeypatch.syspath_prepend(str(CASES.parents[1] / "benchmarks"))
    from native_validation_reference_compare import compare

    source, reference_file = CASES / case["input"], CASES / case["reference"]
    assert hashlib.sha256(source.read_bytes()).hexdigest() == case["input_sha256"]
    assert (
        hashlib.sha256(reference_file.read_bytes()).hexdigest()
        == case["reference_sha256"]
    )
    with gzip.open(reference_file, "rt") as stream:
        reference = json.load(stream)
    dictionaries = [CASES / row["file"] for row in case["dictionaries"]]
    native = validate(
        source,
        chem_data=chemistry,
        nuclear=case["nuclear"],
        use_cdl=case["use_cdl"],
        include_graph=True,
        restraint_dictionaries=dictionaries,
        restraint_links=case["links"],
    )
    native["input"] = dict(sha256=case["input_sha256"])
    native["restraint_dictionaries"] = [
        dict(sha256=hashlib.sha256(file.read_bytes()).hexdigest())
        for file in dictionaries
    ]
    comparison = compare(reference, native)
    assert comparison["all_required_parity"], comparison


def test_full_native_path_does_not_load_reference_packages(chemistry, monkeypatch):
    original = builtins.__import__

    def guarded(name, *args, **kwargs):
        if name.split(".")[0] in {"cctbx", "mmtbx", "iotbx", "libtbx", "scitbx"}:
            raise AssertionError(f"Native validation imported {name}")
        return original(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", guarded)
    source = CASES / "zinc-cys2-his2-ND1-NE2.cif"
    before = source.read_bytes()
    report = validate(source, chem_data=chemistry)
    assert report["complete"]
    assert source.read_bytes() == before


@pytest.mark.parametrize("code", ["1_755", "2_555", "?", "bogus"])
def test_connection_images_cannot_match_the_wrong_bond(code, chemistry):
    parser = pytest.importorskip("iotbx.cif")
    document = parser.reader(file_path=str(CASES / "zinc-declared-image.cif")).model()
    block = next(iter(document.values()))
    block["_struct_conn.ptnr2_symmetry"][0] = code
    result = validate(str(document).encode(), chem_data=chemistry)
    assert not result["complete"]
    assert result["checks"]["bonds"]["status"] == "incomplete"
    assert result["checks"]["clashes"]["status"] == "incomplete"


@pytest.mark.parametrize("kind", ["covale", "covale_sugar", "covale_base", "modres"])
def test_covalent_annotation_subtypes_keep_coverage(kind, chemistry):
    parser = pytest.importorskip("iotbx.cif")
    document = parser.reader(file_path=str(CASES / "zinc-declared-image.cif")).model()
    block = next(iter(document.values()))
    block["_struct_conn.conn_type_id"][0] = kind
    block["_struct_conn.ptnr2_symmetry"][0] = "1_755"
    result = validate(str(document).encode(), chem_data=chemistry)
    assert not result["complete"]
    assert result["checks"]["bonds"]["status"] == "incomplete"
