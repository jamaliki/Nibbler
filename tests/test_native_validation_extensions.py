"""Check native chemistry against pinned, independently captured controls."""

import builtins
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
from tools.fetch_validation_data import DEFAULT_DESTINATION

ROOT = Path(os.environ.get("NIBBLER_VALIDATION_DATA", DEFAULT_DESTINATION))
MANIFEST = (
    json.loads((ROOT / "manifest.json").read_text())
    if (ROOT / "manifest.json").is_file()
    else None
)


@pytest.fixture
def chemistry():
    configured = os.environ.get("NIBBLER_VALIDATION_CHEM_DATA")
    if not configured:
        pytest.skip(
            "Set NIBBLER_VALIDATION_CHEM_DATA to the declared extension chemistry tree"
        )
    return Path(configured)


@pytest.fixture
def assets(chemistry):
    if MANIFEST is None:
        pytest.fail("Run python -m tools.fetch_validation_data before this suite")
    assert MANIFEST["schema"] == "nibbler.validation-data.v1"
    assert len(MANIFEST["cases"]) == 84
    return lambda name: ROOT / MANIFEST["files"]["tests/" + name]


@pytest.fixture
def cases(assets):
    return lambda name: assets("native_validation_extensions/" + name)


@pytest.mark.parametrize(
    "case",
    MANIFEST["cases"] if MANIFEST else [None],
    ids=lambda row: f"{row['family']}/{row['name']}" if row else "data-unavailable",
)
def test_extension_matches_independent_reference(case, chemistry, assets, monkeypatch):
    monkeypatch.syspath_prepend(str(Path(__file__).parents[1] / "benchmarks"))
    from native_validation_reference_compare import compare

    source = assets(case["family"] + "/" + case["input"])
    assert hashlib.sha256(source.read_bytes()).hexdigest() == case["input_sha256"]
    reference_file = assets(case["family"] + "/" + case["reference"])
    assert (
        hashlib.sha256(reference_file.read_bytes()).hexdigest()
        == case["reference_sha256"]
    )
    with gzip.open(reference_file, "rt") as stream:
        reference = json.load(stream)
    dictionaries = [
        assets(case["family"] + "/" + row["file"]) for row in case["dictionaries"]
    ]
    for file, row in zip(dictionaries, case["dictionaries"], strict=True):
        assert hashlib.sha256(file.read_bytes()).hexdigest() == row["sha256"]
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


def test_unknown_component_requires_its_declared_dictionary(chemistry, cases):
    source = cases("custom-ZZZ-model.cif")
    assert not validate(source, chem_data=chemistry)["complete"]
    assert validate(
        source, chem_data=chemistry, restraint_dictionaries=[cases("custom-ZZZ.cif")]
    )["complete"]
    assert not validate(source, chem_data=chemistry)["complete"]


def test_inline_dictionary_keeps_the_source_and_matches_the_private_overlay(
    chemistry, cases
):
    source = cases("custom-ZZZ-model.cif")
    dictionary = cases("custom-ZZZ.cif")
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


def test_dictionary_cache_is_safe_across_threads(chemistry, cases):
    source, dictionary = cases("custom-ZZZ-model.cif"), cases("custom-ZZZ.cif")
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
    chemistry, cases, dictionary, reason
):
    with pytest.raises(ValueError, match=reason):
        validate(
            cases("custom-ZZZ-model.cif"),
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
    chemistry, cases, links, reason
):
    with pytest.raises(ValueError, match=reason):
        validate(
            cases("cross-link.cif"),
            chem_data=chemistry,
            restraint_dictionaries=[cases("cross-link-dictionary.cif")],
            restraint_links=links,
        )
    assert validate(
        cases("cross-link.cif"),
        chem_data=chemistry,
        restraint_dictionaries=[cases("cross-link-dictionary.cif")],
        restraint_links=[("XLINK", 0, 3)],
    )["complete"]


def test_full_native_path_does_not_load_reference_packages(
    chemistry, assets, monkeypatch
):
    original = builtins.__import__

    def guarded(name, *args, **kwargs):
        if name.split(".")[0] in {"cctbx", "mmtbx", "iotbx", "libtbx", "scitbx"}:
            raise AssertionError(f"Native validation imported {name}")
        return original(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", guarded)
    source = assets("native_validation_mcl/zinc-cys2-his2-ND1-NE2.cif")
    before = source.read_bytes()
    assert validate(source, chem_data=chemistry)["complete"]
    assert source.read_bytes() == before


@pytest.mark.parametrize("code", ["1_755", "2_555", "?", "bogus"])
def test_connection_images_cannot_match_the_wrong_bond(code, chemistry, assets):
    parser = pytest.importorskip("iotbx.cif")
    document = parser.reader(
        file_path=str(assets("native_validation_mcl/zinc-declared-image.cif"))
    ).model()
    block = next(iter(document.values()))
    block["_struct_conn.ptnr2_symmetry"][0] = code
    result = validate(str(document).encode(), chem_data=chemistry)
    assert not result["complete"]
    assert result["checks"]["bonds"]["status"] == "incomplete"
    assert result["checks"]["clashes"]["status"] == "incomplete"


@pytest.mark.parametrize("kind", ["covale", "covale_sugar", "covale_base", "modres"])
def test_covalent_annotation_subtypes_keep_coverage(kind, chemistry, assets):
    parser = pytest.importorskip("iotbx.cif")
    document = parser.reader(
        file_path=str(assets("native_validation_mcl/zinc-declared-image.cif"))
    ).model()
    block = next(iter(document.values()))
    block["_struct_conn.conn_type_id"][0] = kind
    block["_struct_conn.ptnr2_symmetry"][0] = "1_755"
    result = validate(str(document).encode(), chem_data=chemistry)
    assert not result["complete"]
    assert result["checks"]["bonds"]["status"] == "incomplete"
