"""Public Python coverage for running Reduce3 on documents in memory."""

from __future__ import annotations

import math
import os
from pathlib import Path

import pytest

import nibbler
from nibbler import ChemistryError, CifDocument
from nibbler.reduce import Reduction

CRAMBIN = Path(__file__).parent / "reduce3" / "crambin_1-5.cif"
CRAMBIN_ENTRY = Path(__file__).parent / "reduce3" / "1crn.cif"
MBO = Path(__file__).parent / "reduce3" / "mbo.cif"

pytestmark = pytest.mark.skipif(
    not nibbler.reduce.available(), reason="built without the reduce3 feature"
)


def _chem_data_or_skip() -> None:
    try:
        nibbler.reduce.run(CRAMBIN, approach="remove")
    except ChemistryError as error:
        if error.code == "REDUCE_MONOMER_LIBRARY_NOT_FOUND":
            pytest.skip("chem_data not found (set REDUCE3_CHEM_DATA)")
        raise


def _atom_count(document: CifDocument) -> int:
    table = nibbler.cif.read(document.to_canonical().encode(), category="atom_site")
    return len(table)


def test_hydrogens_are_added_to_a_parsed_document() -> None:
    _chem_data_or_skip()
    document = nibbler.chomp(CRAMBIN)

    result = nibbler.reduce.run(document, add_flip_movers=True)

    assert isinstance(result, Reduction)
    assert result.document.block_count == 1
    assert result.document.to_canonical().startswith("data_crn5\n")
    assert _atom_count(document) == 33
    assert _atom_count(result.document) == 66
    assert "MoverSingleHydrogenRotator" in result.report


def test_flips_are_considered_by_default() -> None:
    _chem_data_or_skip()
    default = nibbler.reduce.run(CRAMBIN_ENTRY)
    flips = nibbler.reduce.run(CRAMBIN_ENTRY, add_flip_movers=True)
    no_flips = nibbler.reduce.run(CRAMBIN_ENTRY, add_flip_movers=False)

    assert "MoverAmideFlip" in default.report
    assert "MoverAmideFlip" not in no_flips.report
    assert default.document.to_canonical() == flips.document.to_canonical()


def _acid_dihedral(document: CifDocument) -> float:
    """O=C-O-H of the mercuribenzoic acid fixture, in degrees from syn."""
    pyarrow = pytest.importorskip("pyarrow")
    table = nibbler.cif.read(document.to_canonical().encode(), category="atom_site")
    rows = pyarrow.RecordBatchReader.from_stream(table).read_all().to_pylist()
    site = {
        row["label_atom_id"]: tuple(
            float(row[k]) for k in ("Cartn_x", "Cartn_y", "Cartn_z")
        )
        for row in rows
    }

    def sub(a: tuple[float, ...], b: tuple[float, ...]) -> tuple[float, ...]:
        return tuple(x - y for x, y in zip(a, b))

    def cross(a: tuple[float, ...], b: tuple[float, ...]) -> tuple[float, ...]:
        return (
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        )

    def dot(a: tuple[float, ...], b: tuple[float, ...]) -> float:
        return sum(x * y for x, y in zip(a, b))

    p0, p1, p2, p3 = (site[n] for n in ("OZ1", "CZ", "OZ2", "HZ2"))
    b0, b1, b2 = sub(p1, p0), sub(p2, p1), sub(p3, p2)
    n0, n1 = cross(b0, b1), cross(b1, b2)
    m = cross(n0, b1)
    norm = math.sqrt(dot(b1, b1))
    return abs(math.degrees(math.atan2(dot(m, n1) / norm, dot(n0, n1))))


def test_acid_hydrogens_prefer_syn_unless_turned_off() -> None:
    _chem_data_or_skip()
    default = nibbler.reduce.run(MBO)
    off = nibbler.reduce.run(MBO, planar_hydroxyl_preference=0.0, acid_syn_preference=0.0)

    assert _acid_dihedral(default.document) < 10.0
    assert _acid_dihedral(off.document) > 90.0


def test_paths_documents_and_compat_mode_agree() -> None:
    _chem_data_or_skip()
    from_path = nibbler.reduce.run(CRAMBIN, compat=True)
    from_document = nibbler.reduce.run(nibbler.chomp(CRAMBIN), compat=True)

    assert from_path.document.to_canonical() == from_document.document.to_canonical()


def test_the_source_metadata_is_kept() -> None:
    _chem_data_or_skip()
    result = nibbler.reduce.run(CRAMBIN_ENTRY)
    output = result.document.to_canonical().encode()

    before = nibbler.cif.read(CRAMBIN_ENTRY, category="struct_conn")
    after = nibbler.cif.read(output, category="struct_conn")
    assert len(before) == len(after) == 3
    assert after.columns == before.columns
    model = nibbler.mmcif.read(result.document)
    assert model.entry_id == "1CRN"
    assert model.atom_site_count > 600

    # compat mode writes Reduce2's layout, which has no struct_conn
    reduce2_layout = nibbler.reduce.run(CRAMBIN_ENTRY, compat=True)
    canonical = reduce2_layout.document.to_canonical().encode()
    assert len(nibbler.cif.read(canonical, category="struct_conn")) == 0


def test_remove_strips_the_added_hydrogens() -> None:
    _chem_data_or_skip()
    added = nibbler.reduce.run(CRAMBIN)

    removed = nibbler.reduce.run(added.document, approach="remove")

    assert _atom_count(removed.document) == 33


def test_named_library_directory_must_exist(tmp_path: Path) -> None:
    with pytest.raises(ChemistryError) as raised:
        nibbler.reduce.run(CRAMBIN, chem_data=tmp_path)
    assert raised.value.code == "REDUCE_MONOMER_LIBRARY_NOT_FOUND"


def test_documents_without_a_model_are_rejected() -> None:
    _chem_data_or_skip()
    with pytest.raises(ChemistryError) as raised:
        nibbler.reduce.run(b"data_empty\n_entry.id empty\n")
    assert raised.value.code == "REDUCE_INVALID_MODEL"


@pytest.mark.parametrize(
    ("options", "error"),
    [
        ({"approach": "replace"}, ValueError),
        ({"n_terminal_charge": "always"}, ValueError),
        ({"model_id": 0}, ValueError),
        ({"bonded_neighbor_depth": -1}, ValueError),
        ({"add_flip_movers": 1}, TypeError),
        ({"preference_magnitude": "1"}, TypeError),
        ({"acid_syn_preference": "1"}, TypeError),
        ({"probe": {"probe_radius": True}}, TypeError),
        ({"probe": {"set_polar_hydrogen_radius": 1}}, TypeError),
        ({"probe": {"radius": 0.3}}, ValueError),
    ],
)
def test_options_are_validated(
    options: dict[str, object], error: type[Exception]
) -> None:
    with pytest.raises(error):
        nibbler.reduce.run(CRAMBIN, **options)  # type: ignore[arg-type]


def test_explicit_chem_data_directory(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    chem_data = os.environ.get("REDUCE3_CHEM_DATA")
    if chem_data is None:
        pytest.skip("REDUCE3_CHEM_DATA is not set")
    monkeypatch.delenv("REDUCE3_CHEM_DATA")

    result = nibbler.reduce.run(CRAMBIN, chem_data=Path(chem_data))

    assert _atom_count(result.document) == 66
