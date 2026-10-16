"""Check ligand selection and sparse sampling against independent area tests."""

import math

import numpy as np
import pytest

import nibbler


def document(rows, *, entities=True):
    header = "data_sasa\n"
    if entities:
        header += (
            "loop_\n_entity.id\n_entity.type\n"
            "1 polymer\n2 non-polymer\n3 water\n4 branched\n"
        )
    header += (
        "loop_\n"
        + "\n".join(
            "_atom_site." + column
            for column in (
                "type_symbol",
                "label_entity_id",
                "label_comp_id",
                "Cartn_x",
                "Cartn_y",
                "Cartn_z",
                "occupancy",
                "label_alt_id",
                "pdbx_PDB_model_num",
            )
        )
        + "\n"
    )
    return nibbler.cif.read((header + "\n".join(rows) + "\n").encode())


def points(count):
    i = np.arange(count)
    z = 1 - 2 * (i + 0.5) / count
    phi = math.pi * (3 - math.sqrt(5)) * i
    r = np.sqrt(1 - z * z)
    return np.stack([r * np.cos(phi), r * np.sin(phi), z], axis=1)


@pytest.mark.parametrize("count", [1, 31, 256, 1000])
def test_sparse_kernel_matches_brute_force(count):
    rng = np.random.default_rng(729)
    xyz = rng.uniform(-12, 12, (140, 3))
    radii = rng.uniform(1.2, 2.3, len(xyz))
    selected = [139, 4, 18, 0, 85]
    expected = []
    expanded = radii + 1.4
    for i in selected:
        surface = xyz[i] + expanded[i] * points(count)
        squared = np.sum((surface[:, None] - xyz[None, :]) ** 2, axis=-1)
        squared[:, i] = np.inf
        exposed = np.all(squared >= expanded[None, :] ** 2, axis=1)
        expected.append(4 * math.pi * expanded[i] ** 2 * exposed.mean())
    actual = nibbler.geometry.sasa(xyz, radii, selected, point_count=count)
    assert actual == pytest.approx(expected, abs=1e-10)


def test_matches_biotite_with_same_radii_and_sphere():
    import biotite.structure as struc

    rng = np.random.default_rng(21)
    array = struc.AtomArray(100)
    array.coord = rng.uniform(-10, 10, (100, 3))
    array.element[:] = "C"
    array.res_name[:] = "ALA"
    radii = np.full(100, 1.7, dtype=np.float32)
    selected = np.zeros(100, dtype=bool)
    selected[[2, 7, 10, 99]] = True
    reference = struc.sasa(
        array,
        atom_filter=selected,
        ignore_ions=False,
        vdw_radii=radii,
        point_number=256,
        point_distr=lambda n: points(n).astype(np.float32),
    )
    result = nibbler.geometry.sasa(array.coord, radii, np.flatnonzero(selected))
    assert result == pytest.approx(reference[selected], abs=2e-5)


@pytest.mark.parametrize("entry", ["1dpn", "1xan", "1cul", "1a3g"])
def test_real_ligands_match_biotite(entry, preparation_data):
    import biotite.structure as struc
    import pyarrow as pa
    from Bio.PDB.SASA import ATOMIC_RADII

    selected = nibbler.mmcif.select_conformer(
        nibbler.cif.read(preparation_data.entry(entry))
    ).document
    rows = (
        pa.RecordBatchReader.from_stream(selected.project("atom_site"))
        .read_all()
        .to_pylist()
    )
    kept = [
        i
        for i, row in enumerate(rows)
        if row["type_symbol"].upper() not in {"H", "D", "T"}
        and row["label_comp_id"].upper() not in {"HOH", "DOD", "WAT"}
        and float(row["occupancy"]) > 0
    ]
    array = struc.AtomArray(len(kept))
    array.coord = np.array(
        [
            [float(rows[i][key]) for key in ("Cartn_x", "Cartn_y", "Cartn_z")]
            for i in kept
        ],
        dtype=np.float32,
    )
    array.element = np.array([rows[i]["type_symbol"].upper() for i in kept])
    array.res_name[:] = "LIG"
    radii = np.array([ATOMIC_RADII[s] for s in array.element], dtype=np.float32)
    actual = nibbler.mmcif.ligand_sasa(selected)
    assert actual == nibbler.mmcif.ligand_sasa(nibbler.mmcif.read(selected))
    mask = np.array([i in actual for i in kept])
    reference = struc.sasa(
        array,
        atom_filter=mask,
        ignore_ions=False,
        vdw_radii=radii,
        point_number=256,
        point_distr=lambda n: points(n).astype(np.float32),
    )[mask]
    # Float32 reference coordinates can change one boundary-point decision.
    tolerance = 4 * math.pi * (radii[mask] + 1.4) ** 2 / 256 + 1e-4
    assert np.all(np.abs(np.array(list(actual.values())) - reference) <= tolerance)


def test_ligand_context_includes_protein_and_excludes_other_models_water_and_h():
    doc = document(
        [
            "C 2 LIG 0 0 0 1 . 1",
            "C 1 ALA 3.1 0 0 1 . 1",
            "O 3 HOH 0 0 0 1 . 1",
            "H 2 LIG 0 0 0 1 . 1",
            "D 2 LIG 0 0 0 1 . 1",
            "T 2 LIG 0 0 0 1 . 1",
            "C 2 LIG 0 0 0 0 . 1",
            "C 2 LIG 0 0 0 1 . 2",
            "O 4 NAG 100 0 0 1 . 1",
        ]
    )
    result = nibbler.mmcif.ligand_sasa(doc)
    assert set(result) == {0, 8}
    assert (
        result[0] == nibbler.geometry.sasa([[0, 0, 0], [3.1, 0, 0]], [1.7, 1.7], [0])[0]
    )
    assert result[8] == pytest.approx(4 * math.pi * 2.92**2)
    assert nibbler.mmcif.ligand_sasa(doc, model_number=2) == {
        7: pytest.approx(4 * math.pi * 3.1**2)
    }
    assert nibbler.mmcif.ligand_sasa(doc, atom_rows=[]) == {}
    with pytest.raises(ValueError, match="no atoms"):
        nibbler.mmcif.ligand_sasa(doc, model_number=3)


def test_explicit_rows_need_no_entity_metadata_and_keep_source_order():
    doc = document(["C 2 LIG 0 0 0 1 . 1", "O 2 LIG 100 0 0 1 . 1"], entities=False)
    assert list(nibbler.mmcif.ligand_sasa(doc, atom_rows=[1, 0])) == [0, 1]
    with pytest.raises(ValueError, match="entity"):
        nibbler.mmcif.ligand_sasa(doc)
    for rows in ([0, 0], [2]):
        with pytest.raises(ValueError, match="atom_rows"):
            nibbler.mmcif.ligand_sasa(doc, atom_rows=rows)


@pytest.mark.parametrize(
    "field,value,match",
    [
        (0, "X", "radius"),
        (3, "?", "coordinates"),
        (3, "NaN", "coordinates"),
        (6, "?", "occupancy"),
        (6, "-0.1", "occupancy"),
        (6, "1.1", "occupancy"),
        (7, "A", "conformer"),
        (8, "1.5", "model number"),
    ],
)
def test_bad_document_values_fail(field, value, match):
    row = "C 2 LIG 0 0 0 1 . 1".split()
    row[field] = value
    with pytest.raises(ValueError, match=match):
        nibbler.mmcif.ligand_sasa(document([" ".join(row)]))


def test_custom_radii_probe_and_resolution():
    doc = document(["Fe 2 HEM 0 0 0 1 . 1"])
    result = nibbler.mmcif.ligand_sasa(
        doc, radii={"fe": 2.0}, probe_radius=1.2, point_count=31
    )
    assert result == {0: pytest.approx(4 * math.pi * 3.2**2)}
    for radius in (0, -1, math.nan, math.inf):
        with pytest.raises(ValueError, match="radii"):
            nibbler.mmcif.ligand_sasa(doc, radii={"Fe": radius})


def test_competing_alternates_fail_and_selected_labels_are_allowed():
    original = nibbler.cif.read("tests/fixtures/chemistry/ligand_ion_water.cif")
    with pytest.raises(ValueError, match="conformer"):
        nibbler.mmcif.ligand_sasa(original)
    with pytest.raises(ValueError, match="conformer"):
        nibbler.mmcif.ligand_sasa(nibbler.mmcif.read(original))
    selected = nibbler.mmcif.select_conformer(original).document
    assert nibbler.mmcif.ligand_sasa(selected) == nibbler.mmcif.ligand_sasa(
        nibbler.mmcif.read(selected)
    )


@pytest.mark.parametrize(
    "options",
    [
        {"point_count": 0},
        {"point_count": 1_000_001},
        {"probe_radius": -1},
        {"probe_radius": math.nan},
        {"radii": {"X": 0}},
    ],
)
def test_document_and_model_reject_bad_options(options):
    selected = nibbler.mmcif.select_conformer(
        nibbler.cif.read("tests/fixtures/chemistry/ligand_ion_water.cif")
    ).document
    for source in (selected, nibbler.mmcif.read(selected)):
        with pytest.raises(ValueError):
            nibbler.mmcif.ligand_sasa(source, **options)


def test_no_ligands_and_document_model_file_and_binary_inputs(tmp_path):
    doc = document(["C 1 ALA 0 0 0 1 . 1"])
    assert nibbler.mmcif.ligand_sasa(doc) == {}
    original = nibbler.cif.read("tests/fixtures/chemistry/ligand_ion_water.cif")
    selected = nibbler.mmcif.select_conformer(original).document
    expected = nibbler.mmcif.ligand_sasa(selected)
    assert expected
    model = nibbler.mmcif.read(selected)
    assert nibbler.mmcif.ligand_sasa(model) == expected
    assert nibbler.mmcif.ligand_sasa(selected.to_binary()) == expected
    file = tmp_path / "selected.cif.gz"
    nibbler.cif.write(selected, file)
    assert nibbler.mmcif.ligand_sasa(file) == expected
