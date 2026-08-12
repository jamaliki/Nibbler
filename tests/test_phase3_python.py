"""Public Python coverage for Phase 3 schemas, validation, and writing."""

from __future__ import annotations

import gzip
import io
from pathlib import Path

import pytest

import nibbler
from nibbler import CifDocument, CifTable, SchemaError, ValidationError

FIXTURES = Path(__file__).parent / "fixtures"
CHEMISTRY = FIXTURES / "chemistry" / "ligand_ion_water.cif"


def test_schema_guided_projection_uses_dictionary_arrow_types() -> None:
    pyarrow = pytest.importorskip("pyarrow")
    atoms = nibbler.chomp(
        CHEMISTRY,
        category="atom_site",
        columns=[
            "pdbx_PDB_model_num",
            "label_comp_id",
            "Cartn_x",
            "Cartn_y",
            "Cartn_z",
        ],
        schema="pdbx",
    )
    assert isinstance(atoms, CifTable)

    table = pyarrow.RecordBatchReader.from_stream(atoms).read_all()
    assert table.schema.field("pdbx_PDB_model_num").type == pyarrow.int64()
    assert table.schema.field("label_comp_id").type == pyarrow.string()
    assert table.schema.field("Cartn_x").type == pyarrow.float64()
    assert table.column("Cartn_x").to_pylist()[:3] == [0.0, 3.8, 7.0]


def test_schema_guided_scan_retains_types_and_document_schema() -> None:
    pyarrow = pytest.importorskip("pyarrow")
    projected = nibbler.feast(
        [CHEMISTRY, CHEMISTRY],
        category="atom_site",
        columns=["Cartn_x"],
        schema="mmcif",
        workers=2,
    )
    assert [
        pyarrow.RecordBatchReader.from_stream(batch)
        .read_all()
        .schema.field("Cartn_x")
        .type
        for batch in projected
    ] == [pyarrow.float64(), pyarrow.float64()]

    documents = list(nibbler.feast([CHEMISTRY], schema="pdbx", workers=1))
    assert len(documents) == 1
    assert isinstance(documents[0], CifDocument)
    assert documents[0].schema == "pdbx"


def test_schema_projection_rejects_unknown_items_and_bad_numerics() -> None:
    with pytest.raises(SchemaError) as unknown:
        nibbler.chomp(
            CHEMISTRY,
            category="atom_site",
            columns=["not_a_dictionary_item"],
            schema="pdbx",
        )
    assert unknown.value.code == "CIF_SCHEMA_ITEM_UNKNOWN"

    invalid = b"data_bad\nloop_\n_atom_site.Cartn_x\nnot-a-number\n"
    with pytest.raises(SchemaError) as typed:
        nibbler.chomp(
            invalid,
            category="atom_site",
            columns=["Cartn_x"],
            schema="pdbx",
        )
    assert typed.value.code == "CIF_SCHEMA_TYPE"

    with pytest.raises(SchemaError) as scanned:
        nibbler.feast(
            [CHEMISTRY],
            category="atom_site",
            columns=["not_a_dictionary_item"],
            schema="pdbx",
        )
    assert scanned.value.code == "CIF_SCHEMA_ITEM_UNKNOWN"


def test_schema_numeric_missing_values_remain_distinct() -> None:
    pyarrow = pytest.importorskip("pyarrow")
    table = nibbler.chomp(
        b"data_missing\nloop_\n_atom_site.Cartn_x\n?\n.\n",
        category="atom_site",
        columns=["Cartn_x"],
        schema="pdbx",
    )
    assert isinstance(table, CifTable)
    arrow = pyarrow.RecordBatchReader.from_stream(
        table.with_missing("columns")
    ).read_all()
    assert arrow.column("Cartn_x").to_pylist() == [None, None]
    assert arrow.column("Cartn_x__missing_kind").to_pylist() == [1, 2]


def test_sniff_reports_pinned_schema_and_structured_diagnostics() -> None:
    document = nibbler.chomp(CHEMISTRY, schema="pdbx")
    assert isinstance(document, CifDocument)
    report = nibbler.sniff(document)
    assert report.is_valid
    assert report.schema == "pdbx"
    assert report.dictionary_version == "5.416"
    assert "parent-child-links" in report.coverage

    invalid = nibbler.chomp(
        b"data_bad\n_chem_comp.id BAD\n_chem_comp.formula_weight 0\n",
        schema="pdbx",
    )
    invalid_report = nibbler.sniff(invalid)
    assert not invalid_report.is_valid
    assert "CIF_SCHEMA_RANGE" in {
        diagnostic.code for diagnostic in invalid_report.diagnostics
    }
    assert invalid_report.diagnostics[0].context[0] == "block=bad"


def test_sniff_without_schema_is_a_syntax_report() -> None:
    document = nibbler.chomp(b"data_valid\n_entry.id valid\n")
    report = nibbler.sniff(document)
    assert report.is_valid
    assert report.schema is None
    assert report.coverage == ("cif-1.1-syntax",)

    with pytest.raises(TypeError, match="requires CifDocument"):
        nibbler.sniff(object())
    with pytest.raises(SchemaError) as raised:
        nibbler.sniff(document, schema="invented")
    assert raised.value.code == "CIF_SCHEMA_UNKNOWN"
    with pytest.raises(TypeError, match="schema must be a string"):
        nibbler.chomp(CHEMISTRY, schema=1)  # type: ignore[arg-type]


def test_modelcif_artifact_is_embedded_and_lock_verified() -> None:
    document = nibbler.chomp(CHEMISTRY, schema="modelcif")
    report = nibbler.sniff(document)
    assert report.schema == "modelcif"
    assert report.dictionary_version == "1.4.9"


def test_dump_writes_canonical_stream_and_transactional_gzip(tmp_path: Path) -> None:
    document = nibbler.chomp(CHEMISTRY, schema="pdbx")
    stream = io.BytesIO()
    nibbler.dump(document, stream, validate="dictionary")
    stream_document = nibbler.chomp(stream.getvalue())
    assert isinstance(stream_document, CifDocument)
    assert stream_document.to_canonical() == document.to_canonical()

    destination = tmp_path / "roundtrip.cif.gz"
    nibbler.dump(document, destination, validate="dictionary")
    payload = destination.read_bytes()
    assert payload.startswith(b"\x1f\x8b")
    assert gzip.decompress(payload) == document.to_canonical().encode("utf-8")
    reparsed = nibbler.chomp(destination, schema="pdbx")
    assert nibbler.sniff(reparsed).is_valid


def test_dictionary_dump_failure_leaves_existing_destination_unchanged(
    tmp_path: Path,
) -> None:
    invalid = nibbler.chomp(
        b"data_bad\n_chem_comp.id BAD\n_chem_comp.formula_weight 0\n",
        schema="pdbx",
    )
    destination = tmp_path / "existing.cif"
    destination.write_bytes(b"original bytes")

    with pytest.raises(ValidationError):
        nibbler.dump(invalid, destination, validate="dictionary")
    assert destination.read_bytes() == b"original bytes"
    assert list(tmp_path.iterdir()) == [destination]


def test_dump_rejects_unknown_modes() -> None:
    document = nibbler.chomp(b"data_valid\n_entry.id valid\n")
    with pytest.raises(ValueError, match="mode must be"):
        nibbler.dump(document, io.BytesIO(), mode="invented")  # type: ignore[arg-type]
