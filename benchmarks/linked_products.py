"""Audit linked-product cases against a caller-supplied PDB and CCD snapshot."""

from __future__ import annotations

import argparse
import hashlib
import json
import resource
import statistics
import sys
import time
from pathlib import Path
from typing import cast

import nibbler
from benchmarks.preparation_cost import percentile

CASES = (
    "1ivo",
    "4js1",
    "1qfe",
    "1a3g",
    "1j8z",
    "1a1e",
    "1d9d",
    "1rcq",
    "4ndz",
    "6h9v",
    "1twr",
    "6q9t",
    "1dpn",
    "3n95",
    "4aah",
    "6dp5",
    "2odp",
    "1n4e",
    "8qia",
    "6n0a",
    "6w13",
    "1vq7",
    "1cul",
    "1ppv",
    "1u9x",
    "2zz6",
    "8r2z",
    "5lcl",
    "1q29",
    "8q2p",
    "7p1d",
)


def read_components(raw: bytes) -> set[str]:
    return set(nibbler.mmcif.read(raw).component_ids)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--ccd-root", type=Path, required=True)
    parser.add_argument("--repeats", type=int, default=7)
    parser.add_argument("--case", action="append")
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("--repeats must be positive")
    cases = args.case or CASES
    sources = {}
    for case in cases:
        nested = args.archive / case[1:3] / f"{case}.cif.gz"
        sources[case] = (
            nested if nested.is_file() else args.archive / f"{case}.cif.gz"
        ).read_bytes()
    components = set().union(*(read_components(raw) for raw in sources.values()))
    files = []
    for component in sorted(components):
        nested = args.ccd_root / component[0].lower() / f"data_{component}.cif"
        files.append(nested if nested.is_file() else args.ccd_root / f"{component}.cif")
    available = [file for file in files if file.is_file()]
    ccd = b"\n".join(file.read_bytes() for file in available)
    start = time.perf_counter()
    registry = nibbler.components.read(ccd)
    registry_seconds = time.perf_counter() - start
    results = []
    for case, raw in sources.items():
        parse_times = []
        preparation_times = []
        row: dict[str, object] = {
            "pdb_id": case,
            "source_sha256": hashlib.sha256(raw).hexdigest(),
        }
        for _ in range(args.repeats):
            start = time.perf_counter()
            model = nibbler.mmcif.read(raw)
            parse_times.append(time.perf_counter() - start)
            start = time.perf_counter()
            try:
                graph = nibbler.mmcif.prepare_chemistry(model, registry=registry)
            except nibbler.ChemistryError as error:
                row.update(status="rejected", code=error.code, error=str(error))
                break
            preparation_times.append(time.perf_counter() - start)
        if preparation_times:
            data = graph.to_dict()
            row.update(
                status="ok",
                atoms=graph.atom_count,
                bonds=graph.bond_count,
                edits=len(cast(list[object], data["edits"])),
                diagnostics=data["diagnostics"],
                prepare_p50_ms=1000 * statistics.median(preparation_times),
                prepare_p90_ms=1000 * percentile(preparation_times, 0.9),
            )
        row["parse_p50_ms"] = 1000 * statistics.median(parse_times)
        results.append(row)
    print(
        json.dumps(
            {
                "schema": "linked-product-audit.v1",
                "ccd_sha256": hashlib.sha256(ccd).hexdigest(),
                "registry_seconds": registry_seconds,
                "components": len(available),
                "missing_components": [
                    file.stem for file in files if not file.is_file()
                ],
                "peak_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
                * (1 if sys.platform == "darwin" else 1024),
                "cases": results,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
