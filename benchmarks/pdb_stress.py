"""Run staged, correctness-qualified Nibbler benchmarks on real PDB structures."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import statistics
import subprocess
import sys
import time
from collections.abc import Callable, Sequence
from dataclasses import asdict, dataclass
from importlib import import_module
from pathlib import Path
from typing import Any, Protocol, TypeVar, cast

import nibbler
from nibbler import CifDocument, CifTable

from .corpus import DEFAULT_MANIFEST, Structure, load_structures, select_structures
from .run import peak_rss_bytes, require_release_nibbler

DEFAULT_CORPUS = Path(".cache/pdb-stress")
PROJECTED_COLUMNS = ("label_comp_id", "Cartn_x", "Cartn_y", "Cartn_z")
T = TypeVar("T")


@dataclass(frozen=True, slots=True)
class StageResult:
    median_ms: float
    minimum_ms: float
    input_mb_per_second: float | None


@dataclass(frozen=True, slots=True)
class FormatResult:
    format: str
    input_bytes: int
    logical_input_bytes: int
    sha256: str
    read: StageResult
    native_project: StageResult
    arrow_import: StageResult
    warm_cache_read_fraction: float
    full_document: FullDocumentResult | None


@dataclass(frozen=True, slots=True)
class FullDocumentResult:
    stage: StageResult
    block_count: int
    peak_rss_bytes: int
    canonical_bytes: int
    canonical_sha256: str


@dataclass(frozen=True, slots=True)
class StructureResult:
    pdb_id: str
    workload: str
    atom_rows: int
    formats: tuple[FormatResult, ...]
    projections_equal: bool


@dataclass(frozen=True, slots=True)
class StressReport:
    command: str
    platform: str
    machine: str
    python: str
    nibbler: str
    warmups: int
    samples: int
    columns: tuple[str, ...]
    structures: tuple[StructureResult, ...]


class ArrowTable(Protocol):
    def __len__(self) -> int: ...

    def equals(self, other: object) -> bool: ...


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--corpus", type=Path, default=DEFAULT_CORPUS)
    parser.add_argument("--ids", nargs="*", default=[])
    parser.add_argument(
        "--formats",
        choices=("all", "both", "cif", "bcif", "gzip"),
        default="both",
    )
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--samples", type=int, default=3)
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--full-document", action="store_true")
    parser.add_argument("--json", action="store_true", dest="as_json")
    parser.add_argument("--_document-worker", type=Path, help=argparse.SUPPRESS)
    return parser.parse_args(argv)


def timed(
    operation: Callable[[], T],
    *,
    warmups: int,
    samples: int,
    input_bytes: int | None = None,
) -> tuple[StageResult, T]:
    if warmups < 0 or samples < 1:
        raise ValueError("warmups must be non-negative and samples must be positive")
    for _ in range(warmups):
        operation()
    durations = []
    for _ in range(samples):
        start = time.perf_counter_ns()
        result = operation()
        durations.append(time.perf_counter_ns() - start)
    median_ns = statistics.median(durations)
    throughput = None
    if input_bytes is not None:
        throughput = input_bytes / 1_000_000 / (median_ns / 1_000_000_000)
    return (
        StageResult(
            median_ms=median_ns / 1_000_000,
            minimum_ms=min(durations) / 1_000_000,
            input_mb_per_second=throughput,
        ),
        result,
    )


def project(file: Path) -> CifTable:
    table = nibbler.chomp(
        file,
        category="atom_site",
        columns=PROJECTED_COLUMNS,
        schema="pdbx",
    )
    if not isinstance(table, CifTable):
        raise RuntimeError(f"{file}: projection returned a document")
    return table


def run_format(
    file: Path,
    file_format: str,
    logical_input_bytes: int,
    *,
    warmups: int,
    samples: int,
    full_document: bool,
) -> tuple[FormatResult, object]:
    pyarrow = import_module("pyarrow")

    input_bytes = file.stat().st_size
    content = file.read_bytes()
    digest = hashlib.sha256(content).hexdigest()
    read_result, _ = timed(
        file.read_bytes,
        warmups=warmups,
        samples=samples,
        input_bytes=input_bytes,
    )
    native_result, native_table = timed(
        lambda: project(file),
        warmups=warmups,
        samples=samples,
        input_bytes=logical_input_bytes,
    )

    def import_arrow() -> ArrowTable:
        return cast(
            ArrowTable,
            pyarrow.RecordBatchReader.from_stream(native_table).read_all(),
        )

    arrow_result, arrow_table = timed(import_arrow, warmups=warmups, samples=samples)
    read_fraction = read_result.median_ms / native_result.median_ms
    return (
        FormatResult(
            format=file_format,
            input_bytes=input_bytes,
            logical_input_bytes=logical_input_bytes,
            sha256=digest,
            read=read_result,
            native_project=native_result,
            arrow_import=arrow_result,
            warm_cache_read_fraction=read_fraction,
            full_document=(
                run_document_isolated(file, warmups=warmups, samples=samples)
                if full_document
                else None
            ),
        ),
        arrow_table,
    )


def selected_formats(value: str) -> tuple[str, ...]:
    if value == "all":
        return ("cif", "bcif", "cif.gz")
    if value == "both":
        return ("cif", "bcif")
    if value == "gzip":
        return ("cif.gz",)
    return (value,)


def run_structure(
    structure: Structure,
    corpus: Path,
    formats: tuple[str, ...],
    *,
    warmups: int,
    samples: int,
    full_document: bool,
) -> StructureResult:
    pdb_id = structure.pdb_id
    results = []
    arrows: dict[str, ArrowTable] = {}
    for file_format in formats:
        expected = structure.file(file_format)
        logical_input_bytes = (
            structure.cif.bytes if file_format == "cif.gz" else expected.bytes
        )
        result, arrow_table = run_format(
            corpus / f"{pdb_id}.{file_format}",
            file_format,
            logical_input_bytes,
            warmups=warmups,
            samples=samples,
            full_document=full_document,
        )
        if result.input_bytes != expected.bytes or result.sha256 != expected.sha256:
            raise RuntimeError(f"{pdb_id}.{file_format}: corpus manifest mismatch")
        arrow = cast(ArrowTable, arrow_table)
        if len(arrow) != structure.atom_rows:
            raise RuntimeError(f"{pdb_id}.{file_format}: atom-row count mismatch")
        results.append(result)
        arrows[file_format] = arrow
    equal = True
    reference_format = "cif" if "cif" in arrows else next(iter(arrows))
    reference = arrows[reference_format]
    for file_format, arrow_table in arrows.items():
        if file_format == reference_format:
            continue
        equal = equal and reference.equals(arrow_table)
    if not equal:
        raise RuntimeError(f"{pdb_id}: projected formats differ")
    return StructureResult(
        pdb_id=pdb_id,
        workload=structure.workload,
        atom_rows=structure.atom_rows,
        formats=tuple(results),
        projections_equal=equal,
    )


def render(results: Sequence[StructureResult]) -> str:
    lines = []
    for structure in results:
        lines.append(
            f"{structure.pdb_id} ({structure.workload}, {structure.atom_rows:,} atoms)"
        )
        for result in structure.formats:
            throughput = result.native_project.input_mb_per_second
            if throughput is None:
                raise RuntimeError("native projection throughput was not measured")
            lines.append(
                f"  {result.format:7} native "
                f"{result.native_project.median_ms:9.2f} ms  "
                f"{throughput:8.1f} MB/s  Arrow "
                f"{result.arrow_import.median_ms:7.3f} ms  "
                f"read/native {result.warm_cache_read_fraction:5.1%}"
            )
            if result.full_document is not None:
                full = result.full_document
                full_throughput = full.stage.input_mb_per_second
                if full_throughput is None:
                    raise RuntimeError("document throughput was not measured")
                lines.append(
                    f"       full {full.stage.median_ms:9.2f} ms  "
                    f"{full_throughput:8.1f} MB/s  "
                    f"peak RSS {full.peak_rss_bytes / 1_000_000:,.0f} MB"
                )
    return "\n".join(lines)


def report(arguments: argparse.Namespace) -> StressReport:
    require_release_nibbler()
    if arguments.download:
        command = [
            sys.executable,
            "-m",
            "tools.fetch_pdb_corpus",
            "--manifest",
            str(arguments.manifest),
            "--destination",
            str(arguments.corpus),
        ]
        if arguments.ids:
            command.extend(["--ids", *arguments.ids])
        subprocess.run(command, check=True)
    structures = select_structures(load_structures(arguments.manifest), arguments.ids)
    formats = selected_formats(arguments.formats)
    results = tuple(
        run_structure(
            structure,
            arguments.corpus,
            formats,
            warmups=arguments.warmups,
            samples=arguments.samples,
            full_document=arguments.full_document,
        )
        for structure in structures
    )
    return StressReport(
        command=" ".join(sys.argv),
        platform=platform.platform(),
        machine=platform.machine(),
        python=platform.python_version(),
        nibbler=nibbler.__version__,
        warmups=arguments.warmups,
        samples=arguments.samples,
        columns=PROJECTED_COLUMNS,
        structures=results,
    )


def run_document_worker(file: Path, warmups: int, samples: int) -> FullDocumentResult:
    input_bytes = file.stat().st_size
    logical_input_bytes = (
        file.with_suffix("").stat().st_size if file.suffix == ".gz" else input_bytes
    )
    stage, result = timed(
        lambda: nibbler.chomp(file),
        warmups=warmups,
        samples=samples,
        input_bytes=logical_input_bytes,
    )
    if not isinstance(result, CifDocument):
        raise RuntimeError(f"{file}: full parse returned a table")
    parser_peak_rss = peak_rss_bytes()
    canonical = result.to_canonical().encode()
    reparsed = nibbler.chomp(canonical)
    if not isinstance(reparsed, CifDocument):
        raise RuntimeError(f"{file}: canonical round trip returned a table")
    if reparsed.to_canonical().encode() != canonical:
        raise RuntimeError(f"{file}: canonical full-document round trip changed values")
    return FullDocumentResult(
        stage=stage,
        block_count=result.block_count,
        peak_rss_bytes=parser_peak_rss,
        canonical_bytes=len(canonical),
        canonical_sha256=hashlib.sha256(canonical).hexdigest(),
    )


def run_document_isolated(
    file: Path, *, warmups: int, samples: int
) -> FullDocumentResult:
    command = [
        sys.executable,
        "-m",
        "benchmarks.pdb_stress",
        "--_document-worker",
        str(file),
        "--warmups",
        str(warmups),
        "--samples",
        str(samples),
    ]
    completed = subprocess.run(command, check=True, capture_output=True, text=True)
    values: dict[str, Any] = json.loads(completed.stdout)
    return FullDocumentResult(
        stage=StageResult(**values["stage"]),
        block_count=int(values["block_count"]),
        peak_rss_bytes=int(values["peak_rss_bytes"]),
        canonical_bytes=int(values["canonical_bytes"]),
        canonical_sha256=str(values["canonical_sha256"]),
    )


def main(argv: Sequence[str] | None = None) -> int:
    arguments = parse_args(argv)
    if arguments._document_worker is not None:
        require_release_nibbler()
        worker_result = run_document_worker(
            arguments._document_worker, arguments.warmups, arguments.samples
        )
        print(json.dumps(asdict(worker_result)))
        return 0
    report_result = report(arguments)
    if arguments.as_json:
        print(json.dumps(asdict(report_result), indent=2, sort_keys=True))
    else:
        print(render(report_result.structures))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
