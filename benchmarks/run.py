"""Run reproducible projected-CIF benchmarks in isolated processes."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.metadata
import json
import math
import os
import platform
import resource
import subprocess
import sys
import tempfile
import time
from collections.abc import Sequence
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from .baselines import ENGINES, ProjectedRows

DEFAULT_FIXTURE = Path("tests/fixtures/chemistry/ligand_ion_water.cif")
DEFAULT_COLUMNS = ("label_comp_id", "Cartn_x", "Cartn_y", "Cartn_z")
WORKLOADS = {
    "full-document": "Parse a complete logical CIF document without projection.",
    "projected-atom-site": (
        "Read selected atom_site columns without materializing other data."
    ),
    "chemistry-heavy-validation": (
        "Validate PDBx entity, component, scheme, and atom relationships."
    ),
    "canonical-round-trip": (
        "Parse, write canonically, and reparse without semantic loss."
    ),
}


@dataclass(frozen=True, slots=True)
class BenchmarkResult:
    """One isolated engine result with correctness and resource evidence."""

    engine: str
    version: str
    samples_ns: tuple[int, ...]
    median_ns: int
    p95_ns: int
    p99_ns: int
    decompressed_mb_per_second: float
    files_per_second: float
    peak_rss_bytes: int
    row_count: int
    semantic_digest: str
    correct: bool


@dataclass(frozen=True, slots=True)
class BenchmarkReport:
    """A complete benchmark run before JSON serialization."""

    workload: str
    command: str
    file: str
    input_sha256: str
    input_bytes: int
    decompressed_bytes: int
    columns: tuple[str, ...]
    warmups: int
    samples: int
    platform: str
    machine: str
    python: str
    cpu_count: int | None
    results: tuple[BenchmarkResult, ...]


def build_parser() -> argparse.ArgumentParser:
    """Create the benchmark command-line parser."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--list", action="store_true", help="list the frozen plan")
    parser.add_argument("--file", type=Path, default=DEFAULT_FIXTURE)
    parser.add_argument("--engine", choices=("all", *ENGINES), default="all")
    parser.add_argument("--warmups", type=int, default=2)
    parser.add_argument("--samples", type=int, default=10)
    parser.add_argument(
        "--synthetic-rows",
        type=int,
        default=0,
        help="benchmark a generated 24-column atom_site loop with this many rows",
    )
    parser.add_argument("--json", action="store_true", dest="as_json")
    parser.add_argument("--_worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--_expected-digest", default="", help=argparse.SUPPRESS)
    return parser


def render_plan() -> str:
    """Render the stable human-readable benchmark plan."""
    lines = ["Workloads:"]
    lines.extend(f"  {name}: {description}" for name, description in WORKLOADS.items())
    lines.append("Baselines:")
    lines.extend(
        f"  {name} ({engine.distribution}): {engine.purpose}"
        for name, engine in ENGINES.items()
        if name != "nibbler"
    )
    return "\n".join(lines)


def run_worker(
    engine: str,
    file: Path,
    warmups: int,
    samples: int,
    expected_digest: str,
) -> BenchmarkResult:
    """Run one engine in the current isolated process."""
    if warmups < 0 or samples < 1:
        raise ValueError("warmups must be non-negative and samples must be positive")
    if engine == "nibbler":
        require_release_nibbler()
    adapter = ENGINES[engine].project
    for _ in range(warmups):
        adapter(file, DEFAULT_COLUMNS)
    durations = []
    rows: ProjectedRows = ()
    for _ in range(samples):
        start_ns = time.perf_counter_ns()
        rows = adapter(file, DEFAULT_COLUMNS)
        durations.append(time.perf_counter_ns() - start_ns)
    digest = semantic_digest(rows)
    median_ns = percentile(durations, 0.50)
    decompressed_bytes = decompressed_size(file)
    seconds = median_ns / 1_000_000_000
    return BenchmarkResult(
        engine=engine,
        version=importlib.metadata.version(ENGINES[engine].distribution),
        samples_ns=tuple(durations),
        median_ns=median_ns,
        p95_ns=percentile(durations, 0.95),
        p99_ns=percentile(durations, 0.99),
        decompressed_mb_per_second=(decompressed_bytes / 1_000_000) / seconds,
        files_per_second=1.0 / seconds,
        peak_rss_bytes=peak_rss_bytes(),
        row_count=len(rows),
        semantic_digest=digest,
        correct=digest == expected_digest,
    )


def run_isolated(
    engine: str,
    file: Path,
    warmups: int,
    samples: int,
    expected_digest: str,
) -> BenchmarkResult:
    """Run one engine in a clean subprocess and decode its result."""
    command = [
        sys.executable,
        "-m",
        "benchmarks.run",
        "--_worker",
        "--engine",
        engine,
        "--file",
        str(file),
        "--warmups",
        str(warmups),
        "--samples",
        str(samples),
        "--_expected-digest",
        expected_digest,
    ]
    completed = subprocess.run(command, check=True, capture_output=True, text=True)
    values: dict[str, Any] = json.loads(completed.stdout)
    values["samples_ns"] = tuple(values["samples_ns"])
    return BenchmarkResult(**values)


def report(
    file: Path,
    warmups: int,
    samples: int,
    engines: Sequence[str],
) -> BenchmarkReport:
    """Build a machine-readable, correctness-qualified benchmark report."""
    require_release_nibbler()
    reference_rows = ENGINES["nibbler"].project(file, DEFAULT_COLUMNS)
    expected_digest = semantic_digest(reference_rows)
    results = [
        run_isolated(engine, file, warmups, samples, expected_digest)
        for engine in engines
    ]
    if not all(result.correct for result in results):
        mismatches = [result.engine for result in results if not result.correct]
        raise RuntimeError(f"semantic mismatch for engines: {', '.join(mismatches)}")
    return BenchmarkReport(
        workload="projected-atom-site",
        command=" ".join(sys.argv),
        file=str(file),
        input_sha256=hashlib.sha256(file.read_bytes()).hexdigest(),
        input_bytes=file.stat().st_size,
        decompressed_bytes=decompressed_size(file),
        columns=DEFAULT_COLUMNS,
        warmups=warmups,
        samples=samples,
        platform=platform.platform(),
        machine=platform.machine(),
        python=platform.python_version(),
        cpu_count=os.cpu_count(),
        results=tuple(results),
    )


def require_release_nibbler() -> None:
    """Refuse to compare an unoptimized Nibbler extension with release binaries."""
    from nibbler._core import build_profile

    profile = build_profile()
    if profile != "release":
        raise RuntimeError(
            "Nibbler benchmarks require a release extension, "
            f"but the installed extension is {profile!r}; "
            "run `make develop-release` before invoking benchmarks directly"
        )


def semantic_digest(rows: ProjectedRows) -> str:
    """Hash projected cells with unambiguous row and field separators."""
    digest = hashlib.sha256()
    for row in rows:
        for value in row:
            encoded = value.encode("utf-8")
            digest.update(len(encoded).to_bytes(8, "little"))
            digest.update(encoded)
        digest.update(b"\xff")
    return digest.hexdigest()


def percentile(samples: Sequence[int], fraction: float) -> int:
    """Return the nearest-rank percentile for positive integer samples."""
    ordered = sorted(samples)
    rank = max(1, math.ceil(fraction * len(ordered)))
    return ordered[rank - 1]


def decompressed_size(file: Path) -> int:
    """Return logical input bytes, detecting gzip by content."""
    content = file.read_bytes()
    if content.startswith(b"\x1f\x8b"):
        return len(gzip.decompress(content))
    return len(content)


def peak_rss_bytes() -> int:
    """Normalize `ru_maxrss` to bytes on macOS and Unix-like systems."""
    maximum = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(maximum if sys.platform == "darwin" else maximum * 1024)


def generate_synthetic(file: Path, row_count: int) -> None:
    """Write a deterministic wide atom_site projection workload."""
    if row_count < 1:
        raise ValueError("synthetic row count must be positive")
    tags = [
        "label_comp_id",
        "Cartn_x",
        "Cartn_y",
        "Cartn_z",
        *(f"extra_{index}" for index in range(1, 21)),
    ]
    with file.open("w", encoding="utf-8", newline="\n") as output:
        output.write("data_large\nloop_\n")
        output.writelines(f"_atom_site.{tag}\n" for tag in tags)
        for row_index in range(1, row_count + 1):
            output.write(
                f"ALA {row_index / 10:.3f} {row_index / 20:.3f} {row_index / 30:.3f}"
            )
            output.writelines(f" x{column_index}" for column_index in range(1, 21))
            output.write("\n")


def main() -> int:
    """List the plan, run one worker, or orchestrate all selected engines."""
    arguments = build_parser().parse_args()
    if arguments.list:
        print(render_plan())
        return 0
    engines = tuple(ENGINES) if arguments.engine == "all" else (arguments.engine,)
    if arguments._worker:
        if len(engines) != 1:
            raise ValueError("a benchmark worker requires exactly one engine")
        result = run_worker(
            engines[0],
            arguments.file,
            arguments.warmups,
            arguments.samples,
            arguments._expected_digest,
        )
        print(json.dumps(asdict(result)))
        return 0
    if arguments.synthetic_rows:
        with tempfile.TemporaryDirectory(prefix="nibbler-benchmark-") as directory:
            file = Path(directory) / "synthetic-wide.cif"
            generate_synthetic(file, arguments.synthetic_rows)
            output = report(file, arguments.warmups, arguments.samples, engines)
    else:
        output = report(arguments.file, arguments.warmups, arguments.samples, engines)
    print(
        json.dumps(asdict(output), indent=2)
        if arguments.as_json
        else render_results(output)
    )
    return 0


def render_results(report: BenchmarkReport) -> str:
    """Render concise comparable timings after semantic qualification."""
    lines = [
        f"Workload: {report.workload}",
        f"Input: {report.file} ({report.decompressed_bytes} bytes)",
        "engine       median ms    p95 ms    MB/s    peak RSS MiB    correct",
    ]
    lines.extend(
        (
            f"{result.engine:<12} "
            f"{result.median_ns / 1_000_000:>9.3f} "
            f"{result.p95_ns / 1_000_000:>9.3f} "
            f"{result.decompressed_mb_per_second:>7.1f} "
            f"{result.peak_rss_bytes / 1_048_576:>15.1f} "
            f"{result.correct!s:>10}"
        )
        for result in report.results
    )
    return "\n".join(lines)


if __name__ == "__main__":
    raise SystemExit(main())
