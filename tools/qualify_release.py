"""Install one release archive into a fresh environment and smoke-test it."""

from __future__ import annotations

import argparse
import os
import subprocess
import tempfile
import venv
from collections.abc import Sequence
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SMOKE_TEST = ROOT / "tools/installed_smoke.py"
PYARROW = "pyarrow==25.0.0"


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse one wheel or source-distribution path."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    return parser.parse_args(argv)


def environment_python(environment: Path) -> Path:
    """Return the platform-specific Python executable in a virtual environment."""
    if os.name == "nt":
        return environment / "Scripts/python.exe"
    return environment / "bin/python"


def run(command: Sequence[str | Path], *, cwd: Path | None = None) -> None:
    """Run one qualification command with visible, reproducible arguments."""
    rendered = [str(argument) for argument in command]
    print("+", " ".join(rendered), flush=True)
    subprocess.run(rendered, cwd=cwd, check=True)


def main(argv: Sequence[str] | None = None) -> None:
    """Create an isolated environment, install the archive, and run the smoke test."""
    archive = parse_args(argv).archive.resolve()
    if archive.is_dir():
        archives = [*archive.glob("*.whl"), *archive.glob("*.tar.gz")]
        if not archives or len(archives) > 2:
            raise RuntimeError(
                f"expected one wheel and/or one source distribution in {archive}, "
                f"found {len(archives)} archive(s)"
            )
        for candidate in sorted(archives):
            main((str(candidate),))
        return
    if not archive.is_file() or not (
        archive.name.endswith(".whl") or archive.name.endswith(".tar.gz")
    ):
        raise RuntimeError(f"expected one wheel or source distribution, got {archive}")

    with tempfile.TemporaryDirectory(prefix="nibbler-release-") as temporary:
        temporary_root = Path(temporary)
        environment = temporary_root / "environment"
        venv.EnvBuilder(with_pip=True).create(environment)
        python = environment_python(environment)
        run(
            (
                python,
                "-m",
                "pip",
                "install",
                "--disable-pip-version-check",
                "--only-binary=:all:",
                PYARROW,
            )
        )
        run(
            (
                python,
                "-m",
                "pip",
                "install",
                "--disable-pip-version-check",
                "--no-deps",
                archive,
            )
        )
        run((python, "-I", SMOKE_TEST), cwd=temporary_root)


if __name__ == "__main__":
    main()
