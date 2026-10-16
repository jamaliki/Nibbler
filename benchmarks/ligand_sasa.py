"""Measure ligand SASA and its added cost on document and semantic-model reads."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import statistics
import time
from collections.abc import Callable
from pathlib import Path
from typing import Any

import nibbler


def measure(
    operations: dict[str, Callable[[], Any]], warmups: int, samples: int
) -> dict[str, dict[str, float]]:
    """Alternate operation order and retain p50 and nearest-rank p90 times."""
    times: dict[str, list[float]] = {key: [] for key in operations}
    names = list(operations)
    for iteration in range(warmups + samples):
        for key in names[:: -1 if iteration % 2 else 1]:
            started = time.perf_counter_ns()
            result = operations[key]()
            elapsed = (time.perf_counter_ns() - started) / 1e6
            del result
            if iteration >= warmups:
                times[key].append(elapsed)
    return {
        key: {
            "p50_ms": statistics.median(values),
            "p90_ms": sorted(values)[math.ceil(len(values) * 0.9) - 1],
        }
        for key, values in times.items()
    }


def benchmark(
    file: Path, warmups: int, samples: int, point_count: int
) -> dict[str, Any]:
    """Select one conformer outside timing, then compare matched in-memory reads."""
    raw = file.read_bytes()
    document = nibbler.mmcif.select_conformer(nibbler.cif.read(raw)).document
    selected = document.to_canonical().encode()
    model = nibbler.mmcif.read(document)
    areas = nibbler.mmcif.ligand_sasa(model, point_count=point_count)
    assert nibbler.mmcif.ligand_sasa(document, point_count=point_count) == areas

    def area(source: nibbler.CifDocument | nibbler.MmcifModel) -> dict[int, float]:
        return nibbler.mmcif.ligand_sasa(source, point_count=point_count)

    def read_with_area(semantic: bool) -> tuple[object, dict[int, float]]:
        source = (
            nibbler.mmcif.read(selected) if semantic else nibbler.cif.read(selected)
        )
        return source, area(source)

    timings = measure(
        {
            "read_document": lambda: nibbler.cif.read(selected),
            "read_document_sasa": lambda: read_with_area(False),
            "sasa_document": lambda: area(document),
            "read_model": lambda: nibbler.mmcif.read(selected),
            "read_model_sasa": lambda: read_with_area(True),
            "sasa_model": lambda: area(model),
        },
        warmups,
        samples,
    )
    overhead = {
        kind: 100
        * (
            timings[f"read_{kind}_sasa"]["p50_ms"] / timings[f"read_{kind}"]["p50_ms"]
            - 1
        )
        for kind in ("document", "model")
    }
    return {
        "file": file.name,
        "source_sha256": hashlib.sha256(raw).hexdigest(),
        "timed_cif_sha256": hashlib.sha256(selected).hexdigest(),
        "timed_cif_bytes": len(selected),
        "ligand_atoms": len(areas),
        "total_sasa_a2": sum(areas.values()),
        "timings": timings,
        "added_read_percent": overhead,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--file", type=Path, action="append", required=True)
    parser.add_argument("--warmups", type=int, default=5)
    parser.add_argument("--samples", type=int, default=51)
    parser.add_argument("--point-count", type=int, default=256)
    args = parser.parse_args()
    if args.warmups < 0 or args.samples < 1:
        parser.error("warmups must be nonnegative and samples must be positive")
    rows = [
        benchmark(file, args.warmups, args.samples, args.point_count)
        for file in args.file
    ]
    print(
        json.dumps(
            {
                "schema": "nibbler.ligand-sasa-benchmark.v1",
                "platform": platform.platform(),
                "python": platform.python_version(),
                "point_count": args.point_count,
                "probe_radius_a": 1.4,
                "warmups": args.warmups,
                "samples": args.samples,
                "rows": rows,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
