from pathlib import Path

import pytest

import nibbler

FIXTURE = Path("tests/fixtures/modelcif/prediction_with_qa.cif")


def test_modelcif_read_exposes_prediction_metadata() -> None:
    model = nibbler.mmcif.read(FIXTURE, profile="modelcif")

    assert model.profile == "modelcif"
    assert model.entry_id == "NIBBLER_LIGAND_ION_WATER"
    assert model.entity_count == 4
    assert model.prediction_model_count == 1
    assert model.target_entity_count == 1
    assert model.template_count == 1
    assert model.qa_metric_count == 3
    assert model.qa_value_count == 4
    assert model.software_names == ["Rosetta"]
    assert model.qa_metric_names == ["pLDDT", "pTM", "PAE"]
    assert model.qa_metric_modes == ["local", "global", "local-pairwise"]
    assert model.source_document().schema == "modelcif"


def test_modelcif_sniff_reports_pinned_profile_coverage() -> None:
    model = nibbler.mmcif.read(FIXTURE, profile=nibbler.Profile.MODELCIF)
    report = nibbler.sniff(model, profile="modelcif")

    assert report.is_valid
    assert report.schema == "modelcif"
    assert report.dictionary_version == "1.4.9"
    assert "confidence-b-factor-separation" in report.coverage


def test_model_profile_cannot_be_silently_changed() -> None:
    model = nibbler.mmcif.read(FIXTURE, profile="modelcif")

    with pytest.raises(ValueError, match="model profile is 'modelcif'"):
        nibbler.sniff(model, profile="pdbx")


def test_modelcif_dump_round_trips_and_defaults_to_unmirrored_b_values(
    tmp_path: Path,
) -> None:
    pyarrow = pytest.importorskip("pyarrow")
    model = nibbler.mmcif.read(FIXTURE, profile="modelcif")
    output = tmp_path / "prediction.cif"

    nibbler.dump(model, output, profile="modelcif")
    rebuilt = nibbler.mmcif.read(output, profile="modelcif")
    projected = nibbler.chomp(
        output,
        category="atom_site",
        columns=["id", "B_iso_or_equiv"],
        schema="modelcif",
    )
    atoms = pyarrow.RecordBatchReader.from_stream(projected).read_all()

    assert rebuilt.qa_metric_names == model.qa_metric_names
    assert atoms.column("B_iso_or_equiv").to_pylist()[:2] == [10.0, 12.0]


def test_local_qa_mirroring_is_explicit_and_preserves_qa(
    tmp_path: Path,
) -> None:
    pyarrow = pytest.importorskip("pyarrow")
    model = nibbler.mmcif.read(FIXTURE, profile="modelcif")
    output = tmp_path / "mirrored.cif"

    nibbler.dump(
        model,
        output,
        profile="modelcif",
        mirror_local_qa_metric=1,
    )
    rebuilt = nibbler.mmcif.read(output, profile="modelcif")
    projected = nibbler.chomp(
        output,
        category="atom_site",
        columns=["id", "B_iso_or_equiv"],
        schema="modelcif",
    )
    atoms = pyarrow.RecordBatchReader.from_stream(projected).read_all()

    assert atoms.column("B_iso_or_equiv").to_pylist()[:2] == [91.0, 73.0]
    assert rebuilt.qa_metric_names == ["pLDDT", "pTM", "PAE"]
    assert rebuilt.qa_value_count == 4


def test_mirroring_rejects_global_metric(tmp_path: Path) -> None:
    model = nibbler.mmcif.read(FIXTURE, profile="modelcif")

    with pytest.raises(nibbler.WriteError, match="MODELCIF_MIRROR_MODE"):
        nibbler.dump(
            model,
            tmp_path / "invalid.cif",
            profile="modelcif",
            mirror_local_qa_metric=2,
        )


def test_strict_dump_rejects_missing_prediction_provenance(
    tmp_path: Path,
) -> None:
    source = FIXTURE.read_text().replace(
        "1 Rosetta 'model building' 2.1.0",
        "1 Rosetta 'model building' .",
    )
    model = nibbler.mmcif.read(source.encode(), profile="modelcif")

    with pytest.raises(nibbler.ValidationError) as caught:
        nibbler.dump(model, tmp_path / "invalid.cif", profile="modelcif")

    assert any(
        diagnostic.code == "MODELCIF_SOFTWARE_VERSION"
        for diagnostic in caught.value.diagnostics
    )
