"""Reuse one output source without changing targets or query order."""

import gc
import json
import os
from pathlib import Path

import numpy as np
import pytest

import nibbler
from nibbler.observations import PreparedObservationSource
from test_observation_ligand_hydrogens import (
    ULW_PARENTS,
    ULW_ROWS,
    arrays,
    document,
    ligand_document,
    ligand_packer,
    ulw_packer,
)
from test_observation_protein_parents import document as protein_document


def assert_equal(first, second):
    left, left_keys = arrays(first)
    right, right_keys = arrays(second)
    assert left_keys == right_keys and left.keys() == right.keys()
    for name in left:
        assert left[name].dtype == right[name].dtype
        np.testing.assert_array_equal(left[name], right[name], err_msg=name)
        assert right[name].flags.writeable


@pytest.mark.parametrize("complete", [False, True])
@pytest.mark.parametrize("complete_tokens", [False, True])
@pytest.mark.parametrize("with_keys", [False, True])
@pytest.mark.parametrize(
    "hydrogens_unknown,chemical_unknown",
    [(False, False), (False, True), (True, False), (True, True)],
)
def test_prepared_pack_matches_all_target_options(
    complete, complete_tokens, with_keys, hydrogens_unknown, chemical_unknown
):
    native = ligand_packer()
    doc = ligand_document(3)
    source = native.prepare_source(doc)
    _, keys = native.pack(doc)
    options = dict(
        complete=complete,
        complete_tokens=complete_tokens,
        with_keys=with_keys,
        expected_keys=json.dumps(keys[::-1]),
        validate=True,
        target_uncertainties=json.dumps(
            [
                dict(
                    parent_key=keys[0],
                    hydrogens_unknown=hydrogens_unknown,
                    chemical_unknown=chemical_unknown,
                    reasons=["test"],
                )
            ]
        ),
    )
    assert_equal(native.pack(doc, **options), native.pack_prepared(source, **options))


def test_source_owns_its_selected_model_and_cannot_change_document():
    native = ligand_packer()
    doc = document(
        "HETATM C4 C G LIG 1 . 101 . 0 0 0 1\nHETATM C4 C G LIG 1 . 101 . 20 0 0 2\n"
    )
    first = native.prepare_source(doc, model=1)
    second = native.prepare_source(doc, model=2)
    assert_equal(native.pack(doc, model=2), native.pack_prepared(second))
    assert_equal(native.pack(doc, model=1), native.pack_prepared(first))
    del doc
    gc.collect()
    positions, _ = arrays(native.pack_prepared(second))
    assert positions["pos"][0, 0, 0] == 20
    with pytest.raises(TypeError):
        native.pack_prepared(second, model=1)
    with pytest.raises(TypeError):
        native.pack_prepared(second, document=ligand_document(1))
    with pytest.raises((AttributeError, TypeError)):
        first.model = 2
    with pytest.raises(TypeError):
        PreparedObservationSource()


def test_prepared_source_reads_current_parent_tables_and_keeps_owned_buffers():
    native = ulw_packer()
    doc = document(ULW_ROWS)
    source = native.prepare_source(doc)
    with pytest.raises(ValueError, match="more than three hydrogens"):
        native.pack_prepared(source, complete=True, validate=True)
    native.add_hydrogen_parents("ULW", ULW_PARENTS)
    value = native.pack_prepared(source, complete=True, validate=True)
    assert_equal(native.pack(doc, complete=True, validate=True), value)
    del source, doc
    gc.collect()
    tokens, _ = arrays(value)
    assert tokens["mask"][:, 15:].sum() == 4
    tokens["pos"][0, 0, 0] = 5
    assert tokens["pos"][0, 0, 0] == 5


@pytest.mark.parametrize(
    "options",
    [
        dict(expected_keys="not JSON"),
        dict(target_uncertainties="not JSON"),
        dict(expected_keys=json.dumps([["missing"] * 5])),
    ],
)
def test_prepared_target_failures_keep_document_api_errors(options):
    native = ligand_packer()
    doc = ligand_document(3)
    source = native.prepare_source(doc)
    with pytest.raises(ValueError) as expected:
        native.pack(doc, **options)
    with pytest.raises(ValueError) as actual:
        native.pack_prepared(source, **options)
    assert str(actual.value) == str(expected.value)


def test_protein_topology_and_complete_targets_match():
    native = ligand_packer()
    doc = protein_document()
    options = dict(
        complete=True,
        protein_only=True,
        complete_tokens=True,
        topology_policy="entity_sequence",
        validate=True,
        with_keys=False,
    )
    assert_equal(
        native.pack(doc, **options),
        native.pack_prepared(native.prepare_source(doc), **options),
    )
    with pytest.raises(ValueError) as expected:
        native.pack(doc, protein_only=True, topology_policy="bad")
    with pytest.raises(ValueError) as actual:
        native.pack_prepared(
            native.prepare_source(doc), protein_only=True, topology_policy="bad"
        )
    assert str(actual.value) == str(expected.value)


def test_source_and_model_errors_precede_dictionary_loading():
    native = ligand_packer()
    invalid = document("HETATM C4 C G LIG 1 . 101 . bad 0 0 1\n")
    with pytest.raises(ValueError) as expected:
        native.missing_hydrogen_dictionaries(invalid, chem_data="/absent/chem-data")
    with pytest.raises(ValueError) as actual:
        native.prepare_source(invalid)
    assert str(actual.value) == str(expected.value)
    doc = ligand_document(3)
    with pytest.raises(ValueError) as expected:
        native.ligand_components(doc, model=7)
    with pytest.raises(ValueError) as actual:
        native.prepare_source(doc, model=7)
    assert str(actual.value) == str(expected.value)


@pytest.mark.skipif(
    not os.environ.get("REDUCE3_CHEM_DATA"), reason="Set REDUCE3_CHEM_DATA"
)
def test_dictionary_queries_match_each_heavy_atom_snapshot():
    native = ulw_packer()
    chem = os.environ["REDUCE3_CHEM_DATA"]
    docs = [document(ULW_ROWS), document(ULW_ROWS.replace("C14 C G", "UNKNOWN C G"))]
    missing = []
    for doc in docs:
        source = native.prepare_source(doc)
        assert source.ligand_components() == native.ligand_components(doc)
        value = source.missing_hydrogen_dictionaries(chem_data=chem)
        assert value == native.missing_hydrogen_dictionaries(doc, chem_data=chem)
        missing.append(value)
    # An unknown heavy name does not remove an available dictionary.
    assert missing == [[], []]


class PathLike:
    def __init__(self, value):
        self.value = value

    def __fspath__(self):
        return self.value


@pytest.mark.skipif(
    not os.environ.get("REDUCE3_CHEM_DATA"), reason="Set REDUCE3_CHEM_DATA"
)
@pytest.mark.parametrize("convert", [str, Path, PathLike])
def test_supported_chem_paths_keep_queries_and_invalid_model_errors(convert):
    native = ulw_packer()
    chem = convert(os.environ["REDUCE3_CHEM_DATA"])
    doc = document(ULW_ROWS)
    source = native.prepare_source(doc)
    assert source.missing_hydrogen_dictionaries(
        chem_data=chem
    ) == native.missing_hydrogen_dictionaries(doc, chem_data=chem)
    with pytest.raises(ValueError) as expected:
        native.missing_hydrogen_dictionaries(doc, chem_data=chem, model=7)
    with pytest.raises(ValueError) as actual:
        native.prepare_source(doc, model=7)
    assert str(actual.value) == str(expected.value)
