"""Retain reference parity and private dictionary coverage across extensions."""

import gzip
import hashlib
import io
import json
import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pytest

from nibbler import cif
from nibbler.geometry import validate

CASES = Path(__file__).with_name("native_validation_extensions")
MANIFEST = json.loads((CASES / "manifest.json").read_text())


@pytest.fixture
def chemistry():
    configured = os.environ.get("NIBBLER_VALIDATION_CHEM_DATA")
    if not configured:
        pytest.skip(
            "Set NIBBLER_VALIDATION_CHEM_DATA to the declared extension chemistry tree"
        )
    return Path(configured)


@pytest.mark.parametrize("case", MANIFEST["cases"], ids=lambda case: case["name"])
def test_extension_matches_independent_reference(case, chemistry, monkeypatch):
    monkeypatch.syspath_prepend(str(CASES.parents[1] / "benchmarks"))
    from native_validation_reference_compare import compare

    source = CASES / case["input"]
    assert hashlib.sha256(source.read_bytes()).hexdigest() == case["input_sha256"]
    reference_file = CASES / case["reference"]
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
    native["input"] = dict(sha256=hashlib.sha256(source.read_bytes()).hexdigest())
    native["restraint_dictionaries"] = [
        dict(sha256=hashlib.sha256(file.read_bytes()).hexdigest())
        for file in dictionaries
    ]
    result = compare(reference, native)
    assert result["all_required_parity"], result


def test_unknown_component_requires_its_declared_dictionary(chemistry):
    source = CASES / "custom-ZZZ-model.cif"
    assert not validate(source, chem_data=chemistry)["complete"]
    assert validate(
        source, chem_data=chemistry, restraint_dictionaries=[CASES / "custom-ZZZ.cif"]
    )["complete"]
    assert not validate(source, chem_data=chemistry)["complete"]


def test_inline_dictionary_keeps_the_source_and_matches_the_private_overlay(chemistry):
    source = CASES / "custom-ZZZ-model.cif"
    dictionary = CASES / "custom-ZZZ.cif"
    document = cif.read(source.read_bytes() + b"\n" + dictionary.read_bytes())
    before = io.BytesIO()
    cif.write(document, before, mode="preserve")
    result = validate(document, chem_data=chemistry, include_graph=True)
    expected = validate(
        source,
        chem_data=chemistry,
        include_graph=True,
        restraint_dictionaries=[dictionary],
    )
    assert result["complete"]
    assert result["graph"] == expected["graph"]
    assert result["contacts"] == expected["contacts"]
    after = io.BytesIO()
    cif.write(document, after, mode="preserve")
    assert after.getvalue() == before.getvalue()


def test_dictionary_cache_is_safe_across_threads(chemistry):
    source, dictionary = CASES / "custom-ZZZ-model.cif", CASES / "custom-ZZZ.cif"
    expected = validate(
        source, chem_data=chemistry, restraint_dictionaries=[dictionary]
    )
    with ThreadPoolExecutor(max_workers=2) as executor:
        results = list(
            executor.map(
                lambda _: validate(
                    source, chem_data=chemistry, restraint_dictionaries=[dictionary]
                ),
                range(8),
            )
        )
    assert results == [expected] * 8


@pytest.mark.parametrize(
    "dictionary,reason",
    [
        (
            b"data_bad\nloop_\n_lib_bond.atom_type_1\n_lib_bond.atom_type_2\nC N\n",
            "Unsupported energy",
        ),
        (b"data_bad\nsave_targets\n_chem_comp.id ZZZ\nsave_\n", "Save-frame"),
        (
            b"data_bad\nloop_\n_chem_comp_atom.atom_id\n_atom_site.id\nC 1\n",
            "cannot share",
        ),
    ],
)
def test_supplied_inputs_cannot_omit_unsupported_definitions(
    chemistry, dictionary, reason
):
    with pytest.raises(ValueError, match=reason):
        validate(
            CASES / "custom-ZZZ-model.cif",
            chem_data=chemistry,
            restraint_dictionaries=[dictionary],
        )


@pytest.mark.parametrize(
    "links", [[("XLINK", True, 3)], [("XLINK", -1, 3)], ["XLINK"], "XLINK"]
)
def test_link_rows_require_explicit_nonnegative_integers(links):
    with pytest.raises(TypeError, match="restraint_links"):
        validate(b"data_empty\n", chem_data="unused", restraint_links=links)


@pytest.mark.parametrize(
    "links,reason",
    [
        ([("XLINK", 0, 999)], "outside"),
        ([("MISSING", 0, 3)], "Unknown dictionary link"),
        ([("XLINK", 3, 0)], "does not match"),
        ([("XLINK", 0, 3), ("XLINK", 0, 3)], "Duplicate"),
    ],
)
def test_declared_link_errors_do_not_change_the_shared_library(
    chemistry, links, reason
):
    with pytest.raises(ValueError, match=reason):
        validate(
            CASES / "cross-link.cif",
            chem_data=chemistry,
            restraint_dictionaries=[CASES / "cross-link-dictionary.cif"],
            restraint_links=links,
        )
    assert validate(
        CASES / "cross-link.cif",
        chem_data=chemistry,
        restraint_dictionaries=[CASES / "cross-link-dictionary.cif"],
        restraint_links=[("XLINK", 0, 3)],
    )["complete"]
