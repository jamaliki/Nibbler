"""Validate profile fixtures with locked dictionaries and Gemmi."""

from __future__ import annotations

import subprocess
from pathlib import Path

from .fetch_schemas import fetch_all

ROOT = Path(__file__).resolve().parents[1]
PDBX_FIXTURES = (ROOT / "tests" / "fixtures" / "chemistry" / "ligand_ion_water.cif",)


def validate_pdbx(dictionary_file: Path, fixtures: tuple[Path, ...]) -> None:
    """Raise if Gemmi rejects any PDBx fixture against the locked dictionary."""
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


def main() -> int:
    """Fetch locked schemas and validate all current profile fixtures."""
    schemas = fetch_all()
    validate_pdbx(schemas["pdbx"], PDBX_FIXTURES)
    print(f"validated {len(PDBX_FIXTURES)} PDBx fixture(s) with locked dictionary")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
