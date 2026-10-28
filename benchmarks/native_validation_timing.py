"""Time complete native and reference validation in persistent workers."""

import argparse
import hashlib
import importlib.util
import json
import os
import platform
import statistics
import time
from pathlib import Path

import psutil

from nibbler.geometry import validate


def load_reference(script):
    spec = importlib.util.spec_from_file_location("geometry_reference", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.validate


def checked_result(result):
    required = ("bonds", "angles", "chiralities", "planarities", "hydrogens", "clashes")
    if any(
        result.get("checks", {}).get(name, {}).get("status") != "checked"
        for name in required
    ):
        raise ValueError(f"Cannot time incomplete validation: {result.get('checks')}")
    if result.get("complete") is False or result.get("near_cutoff"):
        raise ValueError("Cannot time unsupported coverage or a cutoff ambiguity")
    return {name: result["checks"][name].get("count") for name in required}


def timed(function):
    process = psutil.Process()
    load = os.getloadavg()
    cpu_start = time.process_time()
    started = time.perf_counter()
    result = function()
    elapsed = time.perf_counter() - started
    cpu = time.process_time() - cpu_start
    counts = checked_result(result)
    return dict(
        wall_seconds=elapsed,
        cpu_seconds=cpu,
        load_before=load,
        load_after=os.getloadavg(),
        rss_bytes=process.memory_info().rss,
        counts=counts,
        failure_count=len(result.get("failures", [])),
    )


def summaries(records):
    result = {}
    for name in ("reference", "native"):
        rows = [row for row in records if row["backend"] == name]
        result[name] = {
            metric: dict(
                median=statistics.median(row[metric] for row in rows),
                minimum=min(row[metric] for row in rows),
                maximum=max(row[metric] for row in rows),
            )
            for metric in ("wall_seconds", "cpu_seconds", "rss_bytes")
        }
    return result


def benchmark(args):
    if args.output.exists():
        raise FileExistsError(args.output)
    os.environ["KAVEH_CHEM_DATA"] = str(args.chem_data.resolve())
    reference = load_reference(args.reference_script)
    functions = {
        "reference": lambda: reference(args.input, nuclear=args.nuclear, sigma=5.0),
        "native": lambda: validate(
            args.input, chem_data=args.chem_data, nuclear=args.nuclear
        ),
    }
    for function in functions.values():
        checked_result(function())
    records = []
    for repeat in range(args.repeats):
        order = ("reference", "native") if repeat % 2 == 0 else ("native", "reference")
        for name in order:
            row = timed(functions[name])
            records.append(dict(repeat=repeat, backend=name, **row))
    medians = {
        name: statistics.median(
            row["wall_seconds"] for row in records if row["backend"] == name
        )
        for name in functions
    }
    summary = summaries(records)
    report = dict(
        schema="nibbler.native-validation.timing.v1",
        records=records,
        median_seconds=medians,
        reference_over_native=medians["reference"] / medians["native"],
        summary=summary,
        cpu_reference_over_native=(
            summary["reference"]["cpu_seconds"]["median"]
            / summary["native"]["cpu_seconds"]["median"]
        ),
        input=dict(
            file=str(args.input.resolve()),
            sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),
            nuclear=args.nuclear,
            use_cdl=True,
            sigma_cutoff=5.0,
        ),
        reference_script_sha256=hashlib.sha256(
            args.reference_script.read_bytes()
        ).hexdigest(),
        machine=platform.machine(),
        cpu_count=os.cpu_count(),
        worker_threads={
            name: os.environ.get(name)
            for name in (
                "RAYON_NUM_THREADS",
                "OMP_NUM_THREADS",
                "OPENBLAS_NUM_THREADS",
                "MKL_NUM_THREADS",
            )
        },
        scope="Warm workers include input reads, parsing, interpretation, geometry, "
        "contacts and normal reporting. Reference checks retain the scan's "
        "unchanged-atom guard. Detailed graph capture is outside timed regions. "
        "Atom-level parity must pass separately; timing alone does not prove it.",
    )
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--chem-data", type=Path, required=True)
    parser.add_argument("--reference-script", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--nuclear", action="store_true")
    parser.add_argument("--repeats", type=int, default=10)
    options = parser.parse_args()
    if options.repeats < 4:
        parser.error("Use at least four paired repeats")
    print(json.dumps(benchmark(options), indent=2))
