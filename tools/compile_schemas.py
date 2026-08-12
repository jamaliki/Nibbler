"""Compile locked DDL2 dictionaries into deterministic release artifacts."""

from __future__ import annotations

import os
import subprocess
import tempfile
from pathlib import Path

from .fetch_schemas import ROOT, SchemaLock, fetch_all, load_schema_locks

ARTIFACT_DIRECTORY = ROOT / "schemas" / "compiled"


def compile_schema(lock: SchemaLock, source: Path) -> Path:
    """Compile one verified dictionary and atomically replace its artifact."""
    ARTIFACT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    destination = ARTIFACT_DIRECTORY / lock["artifact"]
    descriptor, temporary_name = tempfile.mkstemp(
        dir=ARTIFACT_DIRECTORY, prefix=f".{lock['artifact']}.", suffix=".tmp"
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    command = [
        "cargo",
        "run",
        "--quiet",
        "--example",
        "compile_schema",
        "--",
        lock["name"],
        lock["dictionary_name"],
        lock["version"],
        lock["sha256"],
        str(source),
        str(temporary),
    ]
    try:
        subprocess.run(command, cwd=ROOT, check=True)
        temporary.chmod(0o644)
        temporary.replace(destination)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    return destination


def compile_all() -> dict[str, Path]:
    """Compile every locked schema from its verified cached source."""
    locks = load_schema_locks()
    if locks["artifact_format_version"] != 1:
        raise RuntimeError(
            f"unsupported artifact format version: {locks['artifact_format_version']}"
        )
    sources = fetch_all()
    return {
        lock["name"]: compile_schema(lock, sources[lock["name"]])
        for lock in locks["schema"]
    }


def main() -> int:
    """Compile all schemas and report their release locations."""
    for name, file in compile_all().items():
        print(f"{name}: {file}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
