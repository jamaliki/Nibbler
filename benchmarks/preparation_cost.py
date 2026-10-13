"""Measure explicit preparation stages on hash-pinned, memory-resident inputs."""

from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import math
import platform
import resource
import statistics
import sys
import time
from pathlib import Path
from typing import Any

import nibbler


def percentile(values: list[float], fraction: float) -> float:
    """Return the nearest-rank percentile."""
    return sorted(values)[max(0, math.ceil(len(values) * fraction) - 1)]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--file", type=Path, action="append", required=True)
    parser.add_argument(
        "--workload",
        choices=(
            "document",
            "model",
            "select",
            "pipeline-select",
            "assembly",
            "product",
        ),
        required=True,
    )
    parser.add_argument("--ccd", type=Path)
    parser.add_argument("--warmups", type=int, default=2)
    parser.add_argument("--samples", type=int, default=21)
    args = parser.parse_args()
    if args.warmups < 0 or args.samples < 1:
        parser.error("warmups must be nonnegative and samples must be positive")
    if args.workload == "product" and args.ccd is None:
        parser.error("product preparation requires --ccd")
    registry = nibbler.components.read(args.ccd) if args.ccd is not None else None
    results = []
    for file in args.file:
        raw = file.read_bytes()
        document = nibbler.cif.read(raw)
        model = (
            nibbler.mmcif.read(document)
            if args.workload in {"assembly", "product"}
            else None
        )

        def operate(
            raw: bytes = raw,
            document: nibbler.CifDocument = document,
            model: nibbler.MmcifModel | None = model,
        ) -> Any:
            if args.workload == "document":
                return nibbler.cif.read(raw)
            if args.workload == "model":
                return nibbler.mmcif.read(raw)
            if args.workload == "select":
                return nibbler.mmcif.select_conformer(document).document
            if args.workload == "pipeline-select":
                client = importlib.import_module("kaveh_data.preparation")
                return client.select_conformer(raw)[0]
            if args.workload == "assembly":
                assert model is not None
                return nibbler.mmcif.assembly(model)
            assert model is not None
            return nibbler.mmcif.prepare_chemistry(model, registry=registry)

        times = []
        output: Any = None
        error = None
        for sample in range(args.warmups + args.samples):
            output = None
            started = time.perf_counter_ns()
            try:
                output = operate()
            except nibbler.ChemistryError as rejected:
                error = rejected.code
            elapsed = (time.perf_counter_ns() - started) / 1e6
            if sample >= args.warmups:
                times.append(elapsed)
        row: dict[str, object] = {
            "file": str(file.resolve()),
            "source_sha256": hashlib.sha256(raw).hexdigest(),
            "source_bytes": len(raw),
            "p50_ms": statistics.median(times),
            "p90_ms": percentile(times, 0.9),
            "samples_ms": times,
            "status": "rejected" if error else "ok",
        }
        if error:
            row["code"] = error
        elif args.workload == "product":
            row.update(atoms=output.atom_count, bonds=output.bond_count)
            encoded = json.dumps(output.to_dict(), sort_keys=True).encode()
            row["output_sha256"] = hashlib.sha256(encoded).hexdigest()
        else:
            result_document = (
                output.source_document() if args.workload == "model" else output
            )
            encoded = result_document.to_canonical().encode()
            row["output_sha256"] = hashlib.sha256(encoded).hexdigest()
            if args.workload == "model":
                row.update(
                    atoms=output.atom_site_count, components=output.component_count
                )
        results.append(row)
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print(
        json.dumps(
            {
                "schema": "preparation-cost.v1",
                "workload": args.workload,
                "platform": platform.platform(),
                "python": platform.python_version(),
                "nibbler_module": nibbler.__file__,
                "warmups": args.warmups,
                "samples": args.samples,
                "peak_rss_bytes": rss if sys.platform == "darwin" else rss * 1024,
                "ccd_sha256": None
                if args.ccd is None
                else hashlib.sha256(args.ccd.read_bytes()).hexdigest(),
                "cases": results,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
