"""Keep weighted metal, glycosidic, RNA and symmetry reference controls."""

import gzip
import hashlib
import json
import os
from pathlib import Path

import pytest

from nibbler.geometry import validate

CASES = Path(__file__).with_name("native_validation_remaining")
MANIFEST = json.loads((CASES / "manifest.json").read_text())


@pytest.mark.parametrize("case", MANIFEST["cases"], ids=lambda case: case["name"])
def test_native_remaining_chemistry_matches_reference(case, monkeypatch):
    chemistry = os.environ.get("NIBBLER_VALIDATION_CHEM_DATA")
    if not chemistry:
        pytest.skip("Set NIBBLER_VALIDATION_CHEM_DATA to the chemistry fixture tree")
    monkeypatch.syspath_prepend(str(CASES.parents[1] / "benchmarks"))
    from native_validation_reference_compare import compare

    source = CASES / case["input"]
    reference_file = CASES / case["reference"]
    assert hashlib.sha256(source.read_bytes()).hexdigest() == case["input_sha256"]
    assert (
        hashlib.sha256(reference_file.read_bytes()).hexdigest()
        == case["reference_sha256"]
    )
    with gzip.open(reference_file, "rt") as stream:
        reference = json.load(stream)
    dictionaries = [CASES / row["file"] for row in case["dictionaries"]]
    result = validate(
        source,
        chem_data=chemistry,
        nuclear=case["nuclear"],
        use_cdl=case["use_cdl"],
        restraint_dictionaries=dictionaries,
        restraint_links=case["links"],
        include_graph=True,
    )
    result["input"] = dict(sha256=case["input_sha256"])
    result["restraint_dictionaries"] = [
        dict(sha256=hashlib.sha256(file.read_bytes()).hexdigest())
        for file in dictionaries
    ]
    comparison = compare(reference, result)
    assert comparison["all_required_parity"], comparison
