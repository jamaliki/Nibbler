"""Public Python coverage for Phase 2 reads and Arrow interchange."""

from __future__ import annotations

import gzip
import io
from pathlib import Path

import pytest

import nibbler
from nibbler import (
    BatchError,
    CifDocument,
    CifTable,
    MissingKind,
    ParseError,
    ProjectionError,
)

FIXTURES = Path(__file__).parent / "fixtures"
MISSING = FIXTURES / "syntax" / "missing_values.cif"
CHEMISTRY = FIXTURES / "chemistry" / "ligand_ion_water.cif"


@pytest.mark.parametrize(
    "source",
    [
        MISSING,
        MISSING.read_bytes(),
        bytearray(MISSING.read_bytes()),
        memoryview(MISSING.read_bytes()),
        io.BytesIO(MISSING.read_bytes()),
    ],
)
def test_chomp_accepts_every_frozen_source_kind(source: object) -> None:
    document = nibbler.chomp(source)  # type: ignore[arg-type]
    assert isinstance(document, CifDocument)
    assert document.block_count == 1
    reparsed = nibbler.chomp(document.to_canonical().encode())
    assert isinstance(reparsed, CifDocument)
    assert reparsed.block_count == 1


def test_chomp_detects_gzip_by_magic_not_filename(tmp_path: Path) -> None:
    compressed = gzip.compress(MISSING.read_bytes())
    no_suffix = tmp_path / "compressed.cif"
    no_suffix.write_bytes(compressed)
    assert isinstance(nibbler.chomp(no_suffix), CifDocument)

    misleading_suffix = tmp_path / "plain.cif.gz"
    misleading_suffix.write_bytes(MISSING.read_bytes())
    assert isinstance(nibbler.chomp(misleading_suffix), CifDocument)


def test_projection_filters_on_unreturned_columns() -> None:
    atoms = nibbler.chomp(
        CHEMISTRY,
        category="atom_site",
        columns=["label_comp_id", "Cartn_x"],
        where={"group_PDB": {"ATOM"}, "pdbx_PDB_model_num": 1},
    )
    assert isinstance(atoms, CifTable)
    assert atoms.category == "atom_site"
    assert atoms.columns == ("label_comp_id", "Cartn_x")
    assert len(atoms) == 2


def test_missing_predicates_and_structured_errors() -> None:
    unknown = nibbler.chomp(
        MISSING,
        category="nibbler_missing",
        columns=["id"],
        where={"unknown_value": MissingKind.UNKNOWN},
    )
    assert isinstance(unknown, CifTable)
    assert len(unknown) == 1

    with pytest.raises(ProjectionError) as projected:
        nibbler.chomp(
            MISSING,
            category="nibbler_missing",
            columns=["absent"],
        )
    assert projected.value.code == "CIF_PROJECTION_MISSING_COLUMN"

    malformed = FIXTURES / "syntax" / "malformed_loop.cif"
    with pytest.raises(ParseError) as parsed:
        nibbler.chomp(malformed)
    assert parsed.value.code == "CIF_LOOP_VALUE_COUNT"
    assert parsed.value.line == 6


def test_arrow_missing_policies() -> None:
    pyarrow = pytest.importorskip("pyarrow")
    polars = pytest.importorskip("polars")
    pandas = pytest.importorskip("pandas")
    table = nibbler.chomp(MISSING, category="nibbler_missing")
    assert isinstance(table, CifTable)

    collapsed = table.to_pyarrow(missing="collapse")
    assert isinstance(collapsed, pyarrow.Table)
    assert collapsed.column("unknown_value").null_count == 2
    metadata = collapsed.schema.field("unknown_value").metadata
    assert metadata[b"nibbler:cif_unknown_count"] == b"1"
    assert metadata[b"nibbler:cif_not_applicable_count"] == b"1"

    columns = table.to_pyarrow(missing="columns")
    assert columns.column("unknown_value__missing_kind").to_pylist() == [1, 2]

    extension = table.to_pyarrow(missing="extension")
    assert extension.column("unknown_value").to_pylist() == [
        {"value": None, "kind": 1},
        {"value": None, "kind": 2},
    ]

    polars_frame = table.to_polars(missing="columns")
    assert isinstance(polars_frame, polars.DataFrame)
    assert polars_frame["unknown_value__missing_kind"].to_list() == [1, 2]

    pandas_frame = table.to_pandas()
    assert isinstance(pandas_frame, pandas.DataFrame)
    assert pandas_frame["unknown_value"].isna().all()


def test_feast_keeps_only_one_worker_window_admitted() -> None:
    admitted = 0

    def sources() -> object:
        nonlocal admitted
        for _ in range(5):
            admitted += 1
            yield MISSING

    result = nibbler.feast(sources(), workers=2)
    assert admitted == 2
    assert isinstance(next(result), CifDocument)
    assert admitted == 3
    assert len(list(result)) == 4
    assert admitted == 5


def test_feast_collects_errors_without_silent_skips() -> None:
    malformed = FIXTURES / "syntax" / "malformed_loop.cif"
    result = nibbler.feast(
        [MISSING, malformed, MISSING],
        category="nibbler_missing",
        columns=["id"],
        workers=2,
        on_error="collect",
    )
    batches = list(result)
    assert len(batches) == 2
    assert all(isinstance(batch, CifTable) for batch in batches)
    assert len(result.errors) == 1
    assert result.errors[0].source_index == 1
    assert result.errors[0].code == "CIF_LOOP_VALUE_COUNT"
    assert result.errors.to_polars()["source_index"].to_list() == [1]

    pyarrow = pytest.importorskip("pyarrow")
    first = batches[0].to_pyarrow()
    assert isinstance(first, pyarrow.Table)
    assert "_nibbler_source" in first.column_names
    assert "_nibbler_block" in first.column_names


def test_feast_raises_the_first_error_in_input_order() -> None:
    malformed = FIXTURES / "syntax" / "malformed_loop.cif"
    result = nibbler.feast([MISSING, malformed, MISSING], workers=3)
    assert isinstance(next(result), CifDocument)
    with pytest.raises(BatchError) as raised:
        next(result)
    assert raised.value.source_index == 1
