"""Public Python coverage for PDBx semantic models and strict writing."""

from __future__ import annotations

import io
from pathlib import Path

import pytest

import nibbler
from nibbler import (
    ChemistryError,
    MmcifModel,
    ValidationError,
)

FIXTURES = Path(__file__).parent / "fixtures"
CHEMISTRY = FIXTURES / "chemistry" / "ligand_ion_water.cif"


def test_mmcif_read_builds_source_backed_semantic_model() -> None:
    model = nibbler.mmcif.read(CHEMISTRY)

    assert isinstance(model, MmcifModel)
    assert isinstance(model, nibbler.mmcif.Model)
    assert model.entry_id == "NIBBLER_LIGAND_ION_WATER"
    assert model.entity_count == 4
    assert model.asym_unit_count == 4
    assert model.component_count == 5
    assert model.atom_site_count == 6
    assert model.connection_count == 1
    assert model.entity_kinds == ("polymer", "non-polymer", "non-polymer", "water")
    assert model.component_ids == ("ALA", "ATP", "HOH", "MSE", "ZN")


def test_model_can_be_built_from_document_and_profile_validated() -> None:
    document = nibbler.chomp(CHEMISTRY)
    model = nibbler.mmcif.read(document)

    document_report = nibbler.sniff(document, profile="pdbx")
    model_report = nibbler.sniff(model, profile="pdbx")
    assert document_report == model_report
    assert model_report.is_valid
    assert model_report.schema == "pdbx"
    assert model_report.dictionary_version == "5.416"
    assert "chemical-component-resolution" in model_report.coverage


def test_profile_dump_is_canonical_and_transactional(tmp_path: Path) -> None:
    model = nibbler.mmcif.read(CHEMISTRY)
    destination = tmp_path / "model.cif"
    nibbler.dump(model, destination, profile="pdbx")

    payload = destination.read_text()
    assert payload.index("_entity.id") < payload.index("_chem_comp.id")
    assert payload.index("_atom_site.group_PDB") < payload.index("_atom_site.Cartn_x")
    reparsed = nibbler.chomp(destination)
    assert nibbler.sniff(reparsed, profile="pdbx").is_valid

    stream = io.BytesIO()
    nibbler.mmcif.write(model, stream, profile="pdbx")
    assert stream.getvalue() == destination.read_bytes()


def test_unresolved_ligand_is_preserved_on_read_but_blocks_strict_output(
    tmp_path: Path,
) -> None:
    source = CHEMISTRY.read_text().replace(
        "ATP non-polymer . \"ADENOSINE-5'-TRIPHOSPHATE\" 'C10 H16 N5 O13 P3' 507.181\n",
        "",
    )
    model = nibbler.mmcif.read(source.encode())
    report = nibbler.sniff(model, profile="pdbx")
    assert "ATP" in model.component_ids
    assert not report.is_valid
    assert "PDBX_COMPONENT_UNRESOLVED" in {
        diagnostic.code for diagnostic in report.errors
    }

    destination = tmp_path / "existing.cif"
    destination.write_bytes(b"original")
    with pytest.raises(ValidationError):
        nibbler.dump(model, destination, profile="pdbx")
    assert destination.read_bytes() == b"original"
    assert list(tmp_path.iterdir()) == [destination]


def test_local_ccd_registry_resolves_unknown_ligand_explicitly(
    tmp_path: Path,
) -> None:
    source = CHEMISTRY.read_text().replace(
        "ATP non-polymer . \"ADENOSINE-5'-TRIPHOSPHATE\" 'C10 H16 N5 O13 P3' 507.181\n",
        "",
    )
    registry = nibbler.components.Registry.from_ccd_cache(
        b"data_ATP\n"
        b"_chem_comp.id ATP\n"
        b"_chem_comp.type non-polymer\n"
        b'_chem_comp.name "ADENOSINE-5\'-TRIPHOSPHATE"\n'
        b"_chem_comp.formula 'C10 H16 N5 O13 P3'\n"
        b"_chem_comp.formula_weight 507.181\n"
        b"loop_\n"
        b"_chem_comp_atom.comp_id\n"
        b"_chem_comp_atom.atom_id\n"
        b"_chem_comp_atom.type_symbol\n"
        b"ATP PG P\n"
    )
    model = nibbler.mmcif.read(source.encode(), registry=registry)

    assert len(registry) == 1
    assert nibbler.sniff(model, profile="pdbx").is_valid
    output = tmp_path / "local-ccd.cif"
    nibbler.dump(model, output, profile="pdbx")
    assert nibbler.sniff(nibbler.chomp(output), profile="pdbx").is_valid


def test_local_ccd_conflict_is_a_structured_error() -> None:
    registry = nibbler.components.Registry.from_ccd_cache(
        b"data_ATP\n_chem_comp.id ATP\n_chem_comp.name 'INCOMPATIBLE NAME'\n"
    )
    with pytest.raises(ChemistryError) as raised:
        nibbler.mmcif.read(CHEMISTRY, registry=registry)
    assert raised.value.code == "PDBX_COMPONENT_CONFLICT"


def test_model_build_failures_are_structured_chemistry_errors() -> None:
    document = nibbler.chomp(b"data_bad\n_atom_site.id 1\n")
    with pytest.raises(ChemistryError) as raised:
        nibbler.mmcif.read(document)
    assert raised.value.code == "PDBX_CATEGORY_REQUIRED"
    assert raised.value.context == ("category=entry",)


def test_semantic_api_rejects_ambiguous_inputs_and_preserves_source_document() -> None:
    document = nibbler.chomp(CHEMISTRY)
    with pytest.raises(TypeError, match="requires MmcifModel"):
        nibbler.dump(document, io.BytesIO(), profile="pdbx")
    stream = io.BytesIO()
    nibbler.dump(
        nibbler.mmcif.read(document),
        stream,
        profile="pdbx",
        mode="preserve",
    )
    assert stream.getvalue() == document.to_preserving().encode()
