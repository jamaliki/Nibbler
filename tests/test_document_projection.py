"""Project parsed and Reduce3 documents through the existing Arrow contract."""

from pathlib import Path

import pytest

import nibbler

pa = pytest.importorskip("pyarrow")


def arrow(table):
    return pa.RecordBatchReader.from_stream(table.with_missing("columns")).read_all()


def test_document_projection_matches_read_and_preserves_source():
    raw = (
        b"data_x\n_entry.id x\nloop_\n_atom_site.id\n_atom_site.Cartn_x\n"
        b"_atom_site.label_alt_id\n1 1.25 ?\n2 -2 .\n"
    )
    document = nibbler.cif.read(raw)
    before = document.to_preserving()
    columns = ["label_alt_id", "Cartn_x", "id"]
    actual = arrow(document.project("atom_site", columns=columns, schema="pdbx"))
    expected = arrow(
        nibbler.cif.read(raw, category="atom_site", columns=columns, schema="pdbx")
    )
    assert actual.equals(expected, check_metadata=True)
    assert actual["id"].type == pa.string()
    assert actual["Cartn_x"].type == pa.float64()
    assert actual["label_alt_id__missing_kind"].to_pylist() == [1, 2]
    assert document.to_preserving() == before


def test_scalar_categories_and_missing_literals():
    document = nibbler.cif.read(b"data_x\n_entry.id ?\n_entry.note '.'\n")
    result = arrow(nibbler.cif.project(document, "entry"))
    assert result["id__missing_kind"].to_pylist() == [1]
    assert result["note"].to_pylist() == ["."]
    assert result["note__missing_kind"].to_pylist() == [0]


def test_project_has_same_structured_column_errors():
    document = nibbler.cif.read(b"data_x\n_entry.id x\n")
    with pytest.raises(nibbler.ProjectionError) as error:
        nibbler.cif.project(document, "entry", columns=["missing"])
    assert error.value.code == "CIF_PROJECTION_MISSING_COLUMN"
    with pytest.raises(TypeError, match="requires CifDocument"):
        nibbler.cif.project(b"data_x\n", "entry")


def test_schema_requires_explicit_selection():
    document = nibbler.cif.read(
        b"data_x\n_atom_site.pdbx_PDB_model_num 1\n", schema="pdbx"
    )
    assert (
        arrow(document.project("atom_site"))["pdbx_PDB_model_num"].type == pa.string()
    )
    assert (
        arrow(nibbler.cif.project(document, "atom_site", schema="mmcif"))[
            "pdbx_PDB_model_num"
        ].type
        == pa.int64()
    )


def test_binary_documents_and_missing_categories():
    document = nibbler.cif.read(b"data_x\n_atom_site.id 1\n_atom_site.Cartn_x 1.25\n")
    decoded = nibbler.cif.read(document.to_binary())
    assert arrow(decoded.project("atom_site", schema="pdbx")).equals(
        arrow(document.project("atom_site", schema="pdbx")), check_metadata=True
    )
    empty = arrow(decoded.project("entity", columns=["id"]))
    assert empty.num_rows == 0
    assert empty.column_names == ["id", "id__missing_kind"]


def test_reduce3_result_projects_without_serializing_for_export():
    if not nibbler.reduce.available():
        pytest.skip("built without Reduce3")
    fixture = Path(__file__).parent / "reduce3" / "crambin_1-5.cif"
    try:
        result = nibbler.reduce.run(fixture, add_flip_movers=True)
    except nibbler.ChemistryError as error:
        if error.code == "REDUCE_MONOMER_LIBRARY_NOT_FOUND":
            pytest.skip("chem_data not found")
        raise
    columns = ["label_atom_id", "type_symbol", "Cartn_x", "Cartn_y", "Cartn_z"]
    actual = arrow(result.document.project("atom_site", columns=columns, schema="pdbx"))
    expected = arrow(
        nibbler.cif.read(
            result.document.to_canonical().encode(),
            category="atom_site",
            columns=columns,
            schema="pdbx",
        )
    )
    assert actual.equals(expected, check_metadata=True)
    assert actual.num_rows == 66
