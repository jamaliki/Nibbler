"""Build a native, corpus-trained PGO wheel without persistent Cargo flags."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
import zipfile
from collections.abc import Sequence
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CORPUS = ROOT / ".cache/pdb-stress"
DEFAULT_MANIFEST = ROOT / "benchmarks/pdb_corpus.toml"
DEFAULT_OUTPUT = ROOT / "dist-pgo"
DEFAULT_WORK_ROOT = ROOT / ".cache/pgo"
ENCODED_FLAG_SEPARATOR = "\x1f"
TRAINING_WORKER = ROOT / "tools/pgo_train.py"


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus", type=Path, default=DEFAULT_CORPUS)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--work-root", type=Path, default=DEFAULT_WORK_ROOT)
    parser.add_argument(
        "--target",
        help="native Rust target triple; defaults to the active rustc host",
    )
    return parser.parse_args(argv)


def run(
    command: Sequence[str | Path],
    *,
    env: dict[str, str] | None = None,
    capture_output: bool = False,
) -> subprocess.CompletedProcess[str]:
    rendered = [str(argument) for argument in command]
    print("+", " ".join(rendered), flush=True)
    return subprocess.run(
        rendered,
        cwd=ROOT,
        env=env,
        check=True,
        text=True,
        capture_output=capture_output,
    )


def command_output(command: Sequence[str | Path]) -> str:
    return run(command, capture_output=True).stdout.strip()


def profile_inputs_digest(manifest: Path) -> str:
    files = [
        ROOT / "Cargo.toml",
        ROOT / "Cargo.lock",
        ROOT / "pyproject.toml",
        manifest,
        TRAINING_WORKER,
    ]
    files.extend(sorted((ROOT / "src").rglob("*.rs")))
    optional_files = (
        ROOT / "build.rs",
        ROOT / "rust-toolchain",
        ROOT / "rust-toolchain.toml",
        ROOT / ".cargo/config",
        ROOT / ".cargo/config.toml",
    )
    files.extend(file for file in optional_files if file.exists())
    digest = hashlib.sha256()
    for file in files:
        label = (
            file.relative_to(ROOT).as_posix()
            if file.is_relative_to(ROOT)
            else file.name
        )
        digest.update(label.encode())
        digest.update(b"\0")
        digest.update(file.read_bytes())
        digest.update(b"\0")
    return digest.hexdigest()


def rustc_details(rustc: str) -> tuple[str, str, Path]:
    version = command_output((rustc, "-vV"))
    fields = dict(line.split(": ", 1) for line in version.splitlines() if ": " in line)
    host = fields.get("host")
    if host is None:
        raise RuntimeError("rustc -vV did not report a host target")
    sysroot = Path(command_output((rustc, "--print", "sysroot")))
    llvm_profdata = sysroot / "lib/rustlib" / host / "bin/llvm-profdata"
    if not llvm_profdata.is_file():
        raise RuntimeError(
            f"the active rustc has no matching llvm-profdata at {llvm_profdata}"
        )
    return version, host, llvm_profdata


def pristine_build_environment() -> dict[str, str]:
    conflicts = tuple(
        name
        for name in (
            "RUSTFLAGS",
            "CARGO_ENCODED_RUSTFLAGS",
            "CARGO_PROFILE_RELEASE_LTO",
            "CARGO_PROFILE_RELEASE_CODEGEN_UNITS",
        )
        if os.environ.get(name)
    )
    if conflicts:
        raise RuntimeError(
            "PGO builds require an uncontaminated release profile; unset "
            + ", ".join(conflicts)
        )
    environment = os.environ.copy()
    environment["CARGO_INCREMENTAL"] = "0"
    return environment


def build_wheel(
    *,
    maturin: str,
    target: str,
    target_dir: Path,
    wheel_dir: Path,
    rust_flags: Sequence[str],
    base_environment: dict[str, str],
) -> Path:
    wheel_dir.mkdir(parents=True, exist_ok=True)
    environment = base_environment | {
        "CARGO_TARGET_DIR": str(target_dir),
        "CARGO_ENCODED_RUSTFLAGS": ENCODED_FLAG_SEPARATOR.join(rust_flags),
    }
    command = (
        maturin,
        "build",
        "--release",
        "--target",
        target,
        "--out",
        wheel_dir,
    )
    run(command, env=environment)
    wheels = tuple(wheel_dir.glob("*.whl"))
    if len(wheels) != 1:
        raise RuntimeError(f"expected one wheel in {wheel_dir}, found {len(wheels)}")
    return wheels[0]


def extract_wheel(wheel: Path, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(wheel) as archive:
        archive.extractall(destination)


def main(argv: Sequence[str] | None = None) -> None:
    arguments = parse_args(argv)
    rustc = shutil.which("rustc")
    maturin = shutil.which("maturin")
    if rustc is None or maturin is None:
        raise RuntimeError("rustc and maturin must be available on PATH")
    rustc_version, host, llvm_profdata = rustc_details(rustc)
    target = host if arguments.target is None else arguments.target
    if target != host:
        raise RuntimeError(
            "PGO training must execute the target wheel: run this tool natively on "
            f"{target}, not the current host {host}"
        )
    corpus = arguments.corpus.resolve()
    manifest = arguments.manifest.resolve()
    if not corpus.is_dir() or not manifest.is_file():
        raise RuntimeError("the PDB corpus and its manifest must exist before training")
    environment = pristine_build_environment()
    inputs_digest = profile_inputs_digest(manifest)
    timestamp = time.strftime("%Y%m%d-%H%M%S", time.gmtime()) + f"-{os.getpid()}"
    work = arguments.work_root.resolve() / f"{target}-{timestamp}"
    raw_profiles = work / "profiles/raw"
    instrumented_wheels = work / "instrumented-wheel"
    instrumented_site = work / "instrumented-site"
    raw_profiles.mkdir(parents=True)

    instrumented_wheel = build_wheel(
        maturin=maturin,
        target=target,
        target_dir=work / "target-instrumented",
        wheel_dir=instrumented_wheels,
        rust_flags=(f"-Cprofile-generate={raw_profiles}",),
        base_environment=environment,
    )
    extract_wheel(instrumented_wheel, instrumented_site)
    training_environment = environment | {
        "LLVM_PROFILE_FILE": str(raw_profiles / "nibbler-%p-%m.profraw"),
        "PYTHONPATH": str(instrumented_site),
    }
    run(
        (
            sys.executable,
            TRAINING_WORKER,
            "--corpus",
            corpus,
            "--manifest",
            manifest,
        ),
        env=training_environment,
    )
    profiles = tuple(raw_profiles.glob("*.profraw"))
    if not profiles:
        raise RuntimeError("training produced no LLVM raw profiles")
    merged_profile = work / "profiles/merged.profdata"
    run((llvm_profdata, "merge", "-o", merged_profile, *profiles))

    if profile_inputs_digest(manifest) != inputs_digest:
        raise RuntimeError(
            "profile inputs changed during PGO training; refusing stale profiles"
        )
    current_rustc_version, current_host, current_profdata = rustc_details(rustc)
    if (current_rustc_version, current_host, current_profdata) != (
        rustc_version,
        host,
        llvm_profdata,
    ):
        raise RuntimeError("the Rust compiler changed during PGO training")

    optimized_wheel = build_wheel(
        maturin=maturin,
        target=target,
        target_dir=work / "target-optimized",
        wheel_dir=work / "optimized-wheel",
        rust_flags=(f"-Cprofile-use={merged_profile}",),
        base_environment=environment,
    )
    final_rustc_version, final_host, final_profdata = rustc_details(rustc)
    if profile_inputs_digest(manifest) != inputs_digest or (
        final_rustc_version,
        final_host,
        final_profdata,
    ) != (rustc_version, host, llvm_profdata):
        raise RuntimeError(
            "profile inputs or compiler changed during the optimized build; "
            "refusing publication"
        )
    output = arguments.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    final_wheel = output / optimized_wheel.name
    shutil.copy2(optimized_wheel, final_wheel)
    metadata = {
        "wheel": str(final_wheel),
        "target": target,
        "rustc": rustc_version,
        "llvm_profdata": str(llvm_profdata),
        "profile_inputs_sha256": inputs_digest,
        "corpus": str(corpus),
        "manifest": str(manifest),
        "raw_profile_count": len(profiles),
        "work_directory": str(work),
    }
    metadata_file = final_wheel.with_suffix(".pgo.json")
    metadata_file.write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
