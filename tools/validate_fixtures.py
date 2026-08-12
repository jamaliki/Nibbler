"""Validate profile fixtures with locked and independent implementations."""

from __future__ import annotations

import subprocess
import tempfile
from importlib import import_module
from pathlib import Path

from .fetch_schemas import fetch_all

ROOT = Path(__file__).resolve().parents[1]
PDBX_FIXTURES = (
    ROOT / "tests" / "fixtures" / "chemistry" / "ligand_ion_water.cif",
    ROOT / "tests" / "fixtures" / "chemistry" / "branched_glycan.cif",
)
MODELCIF_FIXTURES = (
    ROOT / "tests" / "fixtures" / "modelcif" / "prediction_with_qa.cif",
)
ATP_ROW = (
    "ATP non-polymer . \"ADENOSINE-5'-TRIPHOSPHATE\" 'C10 H16 N5 O13 P3' 507.181\n"
)
ATP_CCD = b"""data_ATP
_chem_comp.id ATP
_chem_comp.type non-polymer
_chem_comp.name "ADENOSINE-5'-TRIPHOSPHATE"
_chem_comp.formula 'C10 H16 N5 O13 P3'
_chem_comp.formula_weight 507.181
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
ATP PG P
"""


def validate_gemmi(dictionary_file: Path, fixtures: tuple[Path, ...]) -> None:
    """Raise if Gemmi rejects any fixture against the locked dictionary."""
    command = [
        "gemmi",
        "validate",
        "--quiet",
        "--ddl",
        str(dictionary_file),
        *(str(fixture) for fixture in fixtures),
    ]
    completed = subprocess.run(command, check=False, text=True, capture_output=True)
    if completed.returncode != 0:
        details = completed.stderr.strip() or completed.stdout.strip()
        raise RuntimeError(f"Gemmi dictionary validation failed:\n{details}")


def validate_pdbe(dictionary_file: Path, fixtures: tuple[Path, ...]) -> None:
    """Raise if the pinned PDBe validator reports any dictionary issue."""
    validator = import_module("validate_mmcif")
    for fixture in fixtures:
        errors: list[object] = list(validator.validate(dictionary_file, fixture))
        if not errors:
            continue
        details = "\n".join(
            f"{getattr(error, 'severity', 'error')}: "
            f"{getattr(error, 'item', '<unknown>')}: "
            f"{getattr(error, 'message', error)}"
            for error in errors
        )
        raise RuntimeError(
            f"PDBe dictionary validation failed for {fixture}:\n{details}"
        )


def validate_python_modelcif(fixtures: tuple[Path, ...]) -> None:
    """Raise if Python-ModelCIF cannot construct typed systems from a fixture."""
    reader = import_module("modelcif.reader")
    for fixture in fixtures:
        with fixture.open(encoding="utf-8") as handle:
            systems: list[object] = reader.read(handle, reject_old_file=True)
        if not systems:
            raise RuntimeError(f"Python-ModelCIF read no systems from {fixture}")


def main() -> int:
    """Validate source and Nibbler-written profile goldens externally."""
    import nibbler

    schemas = fetch_all()
    validate_gemmi(schemas["pdbx"], PDBX_FIXTURES)
    validate_pdbe(schemas["pdbx"], PDBX_FIXTURES)
    validate_gemmi(schemas["modelcif"], MODELCIF_FIXTURES)
    validate_pdbe(schemas["modelcif"], MODELCIF_FIXTURES)
    validate_python_modelcif(MODELCIF_FIXTURES)
    with tempfile.TemporaryDirectory(prefix="nibbler-profile-") as temporary_dir:
        generated_pdbx = tuple(
            Path(temporary_dir) / fixture.name for fixture in PDBX_FIXTURES
        )
        for source, destination in zip(PDBX_FIXTURES, generated_pdbx, strict=True):
            model = nibbler.mmcif.read(source)
            nibbler.dump(model, destination, profile="pdbx")
        validate_gemmi(schemas["pdbx"], generated_pdbx)
        validate_pdbe(schemas["pdbx"], generated_pdbx)
        local_ccd_output = Path(temporary_dir) / "local_ccd.cif"
        coordinate_source = PDBX_FIXTURES[0].read_text().replace(ATP_ROW, "")
        registry = nibbler.components.read(ATP_CCD)
        model = nibbler.mmcif.read(coordinate_source.encode(), registry=registry)
        nibbler.dump(model, local_ccd_output, profile="pdbx")
        validate_gemmi(schemas["pdbx"], (local_ccd_output,))
        validate_pdbe(schemas["pdbx"], (local_ccd_output,))
        generated_modelcif: list[Path] = []
        for source in MODELCIF_FIXTURES:
            model = nibbler.mmcif.read(source, profile="modelcif")
            nibbler.sniff(model, profile="modelcif").raise_for_errors()
            canonical_output = Path(temporary_dir) / source.name
            mirrored_output = Path(temporary_dir) / f"mirrored-{source.name}"
            nibbler.dump(model, canonical_output, profile="modelcif")
            nibbler.dump(
                model,
                mirrored_output,
                profile="modelcif",
                mirror_local_qa_metric=1,
            )
            generated_modelcif.extend((canonical_output, mirrored_output))
        generated_modelcif_tuple = tuple(generated_modelcif)
        validate_gemmi(schemas["modelcif"], generated_modelcif_tuple)
        validate_pdbe(schemas["modelcif"], generated_modelcif_tuple)
        validate_python_modelcif(generated_modelcif_tuple)
    print(
        f"validated {len(PDBX_FIXTURES)} PDBx source/writer pair(s) "
        f"and {len(MODELCIF_FIXTURES)} ModelCIF source/writer/mirror set(s); "
        "qualified with locked Nibbler, Gemmi, PDBe, Python-ModelCIF, and "
        "local-CCD materialization"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
