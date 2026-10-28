"""Fixed-mode selection uses deposited model IDs and keeps source identity."""

import json
import os
from hashlib import sha256
from pathlib import Path

import pytest

import nibbler

FIXTURE = Path(__file__).parent / "reduce3/crambin_1-5.cif"


def rows(document):
    pa = pytest.importorskip("pyarrow")
    return (
        pa.RecordBatchReader.from_stream(document.project("atom_site"))
        .read_all()
        .to_pylist()
    )


def document(model_ids):
    gemmi = pytest.importorskip("gemmi")
    source = gemmi.cif.read_file(str(FIXTURE))
    table = source.sole_block().find_mmcif_category("_atom_site.")
    names, original = list(table.tags), [list(row) for row in table]
    model_column = names.index("_atom_site.pdbx_PDB_model_num")
    identifier_column = names.index("_atom_site.id")
    table.erase()
    loop = source.sole_block().init_loop("", names)
    for model in model_ids:
        for atom in original:
            copied = list(atom)
            copied[model_column] = str(model)
            copied[identifier_column] = str(loop.length() + 1)
            loop.add_row(copied)
    return nibbler.cif.read(source.as_string().encode())


def target_keys(document, model_id):
    wanted = {
        (
            row["label_asym_id"],
            str(row["label_seq_id"]),
            row["pdbx_PDB_ins_code"] or "?",
            row["label_comp_id"],
            "protein",
        )
        for row in rows(document)
        if int(row["pdbx_PDB_model_num"]) == model_id
        and row["type_symbol"] not in {"H", "D", "T"}
    }
    return json.dumps(sorted(wanted))


def reduce(document, **options):
    try:
        return nibbler.reduce.run(document, **options)
    except nibbler.ChemistryError as error:
        if error.code == "REDUCE_MONOMER_LIBRARY_NOT_FOUND":
            pytest.skip("Set REDUCE3_CHEM_DATA to test native model selection")
        raise


def heavy_coordinates(document, model_id):
    return {
        (
            row["label_asym_id"],
            row["label_seq_id"],
            row["label_comp_id"],
            row["label_atom_id"],
        ): tuple(float(row[key]) for key in ("Cartn_x", "Cartn_y", "Cartn_z"))
        for row in rows(document)
        if int(row["pdbx_PDB_model_num"]) == model_id
        and row["type_symbol"] not in {"H", "D", "T"}
    }


@pytest.mark.parametrize(
    "model_ids,selected", [([18], 18), ([35], 35), ([35, 18], 18), ([18, 35], 35)]
)
def test_source_model_selection_preserves_ids_coordinates_and_strict_keys(
    model_ids, selected
):
    source = document(model_ids)
    before = source.to_canonical()
    result = reduce(
        source,
        model_id=selected,
        stop_on_any_missing_hydrogen=True,
        required_hydrogen_target_keys=target_keys(source, selected),
    )
    assert {int(row["pdbx_PDB_model_num"]) for row in rows(result.document)} == {
        selected
    }
    assert heavy_coordinates(result.document, selected) == heavy_coordinates(
        source, selected
    )
    assert source.to_canonical() == before
    assert len(rows(result.document)) > len(heavy_coordinates(source, selected))


@pytest.mark.parametrize(
    "model_ids,selected,message",
    [
        ([18], 1, "not present"),
        ([18, 35], 2, "not present"),
        ([18, 35, 18], 18, "ambiguous"),
        ([18, "018"], 18, "ambiguous"),
    ],
)
def test_absent_and_ambiguous_source_ids_fail_without_ordinal_fallback(
    model_ids, selected, message
):
    with pytest.raises(nibbler.ChemistryError, match=message) as caught:
        reduce(document(model_ids), model_id=selected)
    assert caught.value.code == "REDUCE_INVALID_MODEL"


def test_nonsequential_model_keeps_expected_target_keys_strict():
    source = document([18])
    keys = json.loads(target_keys(source, 18))
    keys.append(["absent", "1", "?", "ALA", "protein"])
    with pytest.raises(nibbler.ChemistryError):
        reduce(
            source,
            model_id=18,
            stop_on_any_missing_hydrogen=True,
            required_hydrogen_target_keys=json.dumps(keys),
        )


@pytest.mark.parametrize(
    "relative,model_id,cif_sha256",
    [
        (
            "owner-01/train-000001-0137",
            18,
            "623e09163398b7037c1a2d9ed03bc8a0fd36e168e2dceec36de239eae37dab08",
        ),
        (
            "owner-11/train-000075-0163",
            35,
            "41822b6e18697dd84086d27303b7f5346b2e78ee543a7afbdf1afcd9c9c80278",
        ),
    ],
)
def test_exact_audit_failure_keeps_source_model_and_required_keys(
    relative, model_id, cif_sha256
):
    folder = os.environ.get("NIBBLER_SOURCE_MODEL_FAILURES")
    if not folder:
        pytest.skip(
            "Set NIBBLER_SOURCE_MODEL_FAILURES to the checked census failure directory"
        )
    pa = pytest.importorskip("pyarrow")
    np = pytest.importorskip("numpy")
    base = Path(folder) / relative
    details = json.loads(base.with_suffix(".json").read_text())
    frame = base.with_suffix(".zstd").read_bytes()
    assert (
        sha256(frame).hexdigest()
        == details["frame"]["sha256"]
        == details["record_index"]["sha256"]
    )
    payload = pa.decompress(frame, details["record_index"]["raw_size"], codec="zstd")
    fields = details["fields"]
    cursor = 8 * len(fields)
    lengths = np.frombuffer(payload, dtype="<u8", count=len(fields))
    assert int(lengths.sum()) == len(payload) - cursor
    record = {}
    for field, size in zip(fields, map(int, lengths), strict=True):
        dtype = np.dtype(field["dtype"])
        record[field["name"]] = np.frombuffer(
            payload, dtype=dtype, count=size // dtype.itemsize, offset=cursor
        ).reshape(field["shape"])
        cursor += size
    raw = record["cif"].tobytes()
    assert sha256(raw).hexdigest() == cif_sha256 and record["model"].tolist() == [
        model_id
    ]
    source = nibbler.cif.read(raw)
    before = source.to_canonical()
    result = reduce(
        source,
        model_id=model_id,
        stop_on_any_missing_hydrogen=True,
        required_hydrogen_target_keys=record["token_keys"].tobytes().decode(),
    )
    assert {int(row["pdbx_PDB_model_num"]) for row in rows(result.document)} == {
        model_id
    }
    gemmi = pytest.importorskip("gemmi")
    reference = gemmi.cif.read_string(raw.decode())
    table = reference.sole_block().find_mmcif_category("_atom_site.")
    column = list(table.tags).index("_atom_site.pdbx_PDB_model_num")
    for atom in table:
        atom[column] = "1"
    reference = reduce(
        nibbler.cif.read(reference.as_string().encode()),
        model_id=1,
        stop_on_any_missing_hydrogen=True,
        required_hydrogen_target_keys=record["token_keys"].tobytes().decode(),
    )

    def without_model(document):
        return [
            {key: value for key, value in atom.items() if key != "pdbx_PDB_model_num"}
            for atom in rows(document)
        ]

    assert without_model(result.document) == without_model(reference.document)
    assert source.to_canonical() == before
