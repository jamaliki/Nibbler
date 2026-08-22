"""Check a PDB stress report against broad hosted-runner guardrails."""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

if sys.version_info >= (3, 11):
    import tomllib
else:  # pragma: no cover - Python 3.10
    import tomli as tomllib

DEFAULT_GUARDRAILS = Path("benchmarks/pdb_guardrails.toml")


@dataclass(frozen=True, slots=True)
class Guardrails:
    """Hardware-independent floors applied only to large logical inputs."""

    minimum_logical_input_bytes: int
    maximum_peak_rss_bytes: int
    peak_rss_fixed_allowance_bytes: int
    maximum_peak_rss_to_document_bytes: float
    minimum_projection: Mapping[str, float]
    minimum_full_document: Mapping[str, float]


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse the report, guardrail, and optional summary paths."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("report", type=Path)
    parser.add_argument("--guardrails", type=Path, default=DEFAULT_GUARDRAILS)
    parser.add_argument("--summary", type=Path)
    return parser.parse_args(argv)


def numeric_map(value: object, label: str) -> dict[str, float]:
    """Type-check one TOML table of format-specific throughput floors."""
    if not isinstance(value, dict):
        raise RuntimeError(f"{label} must be a table")
    result = {}
    for key, number in value.items():
        if (
            not isinstance(key, str)
            or isinstance(number, bool)
            or not isinstance(number, (int, float))
        ):
            raise RuntimeError(f"{label} contains a non-numeric entry")
        result[key] = float(number)
    return result


def load_guardrails(file: Path) -> Guardrails:
    """Load and type-check the versioned guardrail document."""
    with file.open("rb") as stream:
        raw: dict[str, Any] = tomllib.load(stream)
    if raw.get("version") != 2:
        raise RuntimeError(f"{file}: unsupported guardrail version")

    def integer(name: str) -> int:
        value = raw.get(name)
        if isinstance(value, bool) or not isinstance(value, int):
            raise RuntimeError(f"{file}: {name} must be an integer")
        return value

    def number(name: str) -> float:
        value = raw.get(name)
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise RuntimeError(f"{file}: {name} must be numeric")
        return float(value)

    return Guardrails(
        minimum_logical_input_bytes=integer("minimum_logical_input_bytes"),
        maximum_peak_rss_bytes=integer("maximum_peak_rss_bytes"),
        peak_rss_fixed_allowance_bytes=integer("peak_rss_fixed_allowance_bytes"),
        maximum_peak_rss_to_document_bytes=number("maximum_peak_rss_to_document_bytes"),
        minimum_projection=numeric_map(
            raw.get("minimum_projection_mb_per_second"),
            "minimum_projection_mb_per_second",
        ),
        minimum_full_document=numeric_map(
            raw.get("minimum_full_document_mb_per_second"),
            "minimum_full_document_mb_per_second",
        ),
    )


def load_report(file: Path) -> dict[str, Any]:
    """Load one benchmark JSON object."""
    value: Any = json.loads(file.read_text())
    if not isinstance(value, dict):
        raise RuntimeError(f"{file}: report must be an object")
    return value


def check(report: Mapping[str, Any], guardrails: Guardrails) -> tuple[str, ...]:
    """Return every correctness, throughput, or memory violation."""
    failures = []
    structures = report.get("structures")
    if not isinstance(structures, list) or not structures:
        return ("report contains no structures",)
    for structure in structures:
        if not isinstance(structure, dict):
            failures.append("report contains a malformed structure")
            continue
        pdb_id = str(structure.get("pdb_id"))
        if structure.get("projections_equal") is not True:
            failures.append(f"{pdb_id}: CIF/BinaryCIF/gzip projections differ")
        formats = structure.get("formats")
        if not isinstance(formats, list):
            failures.append(f"{pdb_id}: formats are missing")
            continue
        for result in formats:
            if not isinstance(result, dict):
                failures.append(f"{pdb_id}: malformed format result")
                continue
            file_format = str(result.get("format"))
            logical_bytes = result.get("logical_input_bytes")
            if not isinstance(logical_bytes, int):
                failures.append(f"{pdb_id}.{file_format}: logical size is missing")
                continue
            if logical_bytes < guardrails.minimum_logical_input_bytes:
                continue
            projection = result.get("native_project")
            full = result.get("full_document")
            if not isinstance(projection, dict) or not isinstance(full, dict):
                failures.append(f"{pdb_id}.{file_format}: measured stages are missing")
                continue
            projection_rate = projection.get("input_mb_per_second")
            full_stage = full.get("stage")
            peak_rss = full.get("peak_rss_bytes")
            canonical_bytes = full.get("canonical_bytes")
            if not isinstance(full_stage, dict):
                failures.append(f"{pdb_id}.{file_format}: full stage is missing")
                continue
            full_rate = full_stage.get("input_mb_per_second")
            expected_projection = guardrails.minimum_projection.get(file_format)
            expected_full = guardrails.minimum_full_document.get(file_format)
            if (
                not isinstance(projection_rate, (int, float))
                or expected_projection is None
            ):
                failures.append(f"{pdb_id}.{file_format}: projection rate is missing")
            elif projection_rate < expected_projection:
                failures.append(
                    f"{pdb_id}.{file_format}: projection {projection_rate:.1f} MB/s "
                    f"is below {expected_projection:.1f} MB/s"
                )
            if not isinstance(full_rate, (int, float)) or expected_full is None:
                failures.append(
                    f"{pdb_id}.{file_format}: full-document rate is missing"
                )
            elif full_rate < expected_full:
                failures.append(
                    f"{pdb_id}.{file_format}: full document {full_rate:.1f} MB/s "
                    f"is below {expected_full:.1f} MB/s"
                )
            if not isinstance(peak_rss, int):
                failures.append(f"{pdb_id}.{file_format}: peak RSS is missing")
            elif not isinstance(canonical_bytes, int):
                failures.append(
                    f"{pdb_id}.{file_format}: canonical document size is missing"
                )
            else:
                document_bytes = max(logical_bytes, canonical_bytes)
                relative_limit = (
                    guardrails.peak_rss_fixed_allowance_bytes
                    + guardrails.maximum_peak_rss_to_document_bytes * document_bytes
                )
                memory_limit = min(guardrails.maximum_peak_rss_bytes, relative_limit)
                if peak_rss <= memory_limit:
                    continue
                failures.append(
                    f"{pdb_id}.{file_format}: peak RSS {peak_rss / 1e9:.2f} GB "
                    f"exceeds the {memory_limit / 1e9:.2f} GB guardrail"
                )
    return tuple(failures)


def markdown(report: Mapping[str, Any], failures: Sequence[str]) -> str:
    """Render one compact GitHub summary from a validated report shape."""
    lines = [
        "## PDB regression",
        "",
        f"Platform: `{report.get('platform')}`; Python: `{report.get('python')}`; "
        f"Nibbler: `{report.get('nibbler')}`",
        "",
        "| PDB | Format | Projection MB/s | Full MB/s | Peak RSS MB |",
        "| --- | --- | ---: | ---: | ---: |",
    ]
    for structure in report.get("structures", []):
        for result in structure.get("formats", []):
            projection = result["native_project"]["input_mb_per_second"]
            full = result["full_document"]
            full_rate = full["stage"]["input_mb_per_second"]
            lines.append(
                f"| {structure['pdb_id']} | {result['format']} | {projection:.1f} | "
                f"{full_rate:.1f} | {full['peak_rss_bytes'] / 1e6:.0f} |"
            )
    lines.extend(["", "### Guardrails", ""])
    lines.extend(f"- {failure}" for failure in failures)
    if not failures:
        lines.append("- All correctness, throughput, and memory guardrails passed.")
    return "\n".join(lines) + "\n"


def main(argv: Sequence[str] | None = None) -> int:
    """Validate the report, write a summary, and fail on any violation."""
    arguments = parse_args(argv)
    report = load_report(arguments.report)
    failures = check(report, load_guardrails(arguments.guardrails))
    rendered = markdown(report, failures)
    print(rendered)
    if arguments.summary is not None:
        with arguments.summary.open("a") as stream:
            stream.write(rendered)
    if failures:
        raise RuntimeError("PDB regression guardrails failed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
