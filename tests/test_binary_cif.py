"""Public Python coverage for BinaryCIF and preserving output."""

from __future__ import annotations

import gzip
import io
from pathlib import Path

import pytest

import nibbler
from nibbler import CifDocument, CifTable

CHEMISTRY = Path(__file__).parent / "fixtures" / "chemistry" / "ligand_ion_water.cif"

SOURCE = b"""data_binary
_entry.id binary
loop_
_atom_site.label_comp_id
_atom_site.Cartn_x
ALA 1.25
? .
"""


def _project(source: object) -> CifTable:
    table = nibbler.chomp(
        source,
        category="atom_site",
        columns=["label_comp_id", "Cartn_x"],
    )
    assert isinstance(table, CifTable)
    return table


def test_binary_cif_content_detection_and_projection_equality() -> None:
    document = nibbler.chomp(SOURCE)
    assert isinstance(document, CifDocument)
    binary = document.to_binary()

    decoded = nibbler.chomp(binary)
    assert isinstance(decoded, CifDocument)
    assert decoded.block_count == document.block_count
    pyarrow = pytest.importorskip("pyarrow")
    binary_table = pyarrow.RecordBatchReader.from_stream(_project(binary)).read_all()
    source_table = pyarrow.RecordBatchReader.from_stream(_project(SOURCE)).read_all()
    assert binary_table.equals(source_table)


def test_dump_infers_binary_cif_and_gzip_from_destination(tmp_path: Path) -> None:
    document = nibbler.chomp(SOURCE)
    assert isinstance(document, CifDocument)

    destination = tmp_path / "roundtrip.bcif.gz"
    nibbler.dump(document, destination)
    payload = destination.read_bytes()
    assert payload.startswith(b"\x1f\x8b")
    assert gzip.decompress(payload).startswith(b"\x83")
    pyarrow = pytest.importorskip("pyarrow")
    destination_table = pyarrow.RecordBatchReader.from_stream(
        _project(destination)
    ).read_all()
    source_table = pyarrow.RecordBatchReader.from_stream(_project(SOURCE)).read_all()
    assert destination_table.equals(source_table)


def test_binary_stream_is_explicit_and_preserve_mode_retains_lexemes() -> None:
    document = nibbler.chomp(b"data_quotes\n_a.one 'alpha beta'\n_a.two \"gamma\"\n")
    assert isinstance(document, CifDocument)

    binary_stream = io.BytesIO()
    nibbler.dump(document, binary_stream, format="bcif")
    assert binary_stream.getvalue().startswith(b"\x83")

    preserve_stream = io.BytesIO()
    nibbler.dump(document, preserve_stream, mode="preserve")
    assert b"_a.one 'alpha beta'\n" in preserve_stream.getvalue()
    assert b'_a.two "gamma"\n' in preserve_stream.getvalue()

    with pytest.raises(ValueError, match="requires format='cif'"):
        nibbler.dump(document, io.BytesIO(), mode="preserve", format="bcif")


def test_malformed_binary_cif_reports_a_stable_error_code() -> None:
    with pytest.raises(nibbler.ParseError) as raised:
        nibbler.chomp(b"\x83")
    assert raised.value.code == "CIF_BINARY_CONTAINER"


def test_binary_cif_feeds_the_same_pdbx_semantic_model() -> None:
    document = nibbler.chomp(CHEMISTRY)
    assert isinstance(document, CifDocument)
    text_model = nibbler.mmcif.read(document)
    binary_model = nibbler.mmcif.read(document.to_binary())

    assert binary_model.entry_id == text_model.entry_id
    assert binary_model.entity_count == text_model.entity_count
    assert binary_model.component_ids == text_model.component_ids
    assert binary_model.atom_site_count == text_model.atom_site_count
