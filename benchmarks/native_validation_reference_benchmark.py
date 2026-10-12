"""Measure fresh-process native and CCTBX capture in alternating order."""

import argparse
import json
import statistics
import subprocess
import time
from pathlib import Path


def run(command, output, timeout):
    started = time.perf_counter()
    result = subprocess.run(
        list(map(str, command)), capture_output=True, text=True, timeout=timeout
    )
    elapsed = time.perf_counter() - started
    report = json.loads(output.read_text()) if output.exists() else None
    if result.returncode or report is None or not report.get("complete"):
        raise ValueError(
            f"Capture failed: {result.stderr[-4096:]} {report and report.get('error')}"
        )
    return dict(
        seconds=elapsed,
        output_bytes=output.stat().st_size,
        capture_seconds=report.get("elapsed_seconds"),
        command=list(map(str, command)),
    )


def benchmark(args):
    args.output.mkdir(parents=True, exist_ok=False)
    reference = Path(__file__).with_name("native_validation_reference_capture.py")
    records = []
    for repeat in range(args.repeats):
        for kind in (
            ("reference", "native") if repeat % 2 == 0 else ("native", "reference")
        ):
            output = args.output / f"{repeat}-{kind}.json"
            interpreter = (
                args.reference_python if kind == "reference" else args.native_python
            )
            script = reference if kind == "reference" else args.native_capture
            command = [
                interpreter,
                script,
                "--input",
                args.input,
                "--chem-data",
                args.chem_data,
                "--output",
                output,
            ]
            if args.nuclear:
                command.append("--nuclear")
            row = run(command, output, args.timeout)
            row.update(repeat=repeat, kind=kind)
            records.append(row)
    medians = {
        kind: statistics.median(
            row["seconds"] for row in records if row["kind"] == kind
        )
        for kind in ("reference", "native")
    }
    report = dict(
        records=records,
        median_seconds=medians,
        reference_over_native_ratio=medians["reference"] / medians["native"],
        scope="Fresh processes include imports, setup, validation and JSON output. "
        "AB/BA order limits order effects. Reference capture adds instrumentation. "
        "Validate parity first and report output sizes. This measures startup too; "
        "it does not measure persistent-worker speed.",
    )
    (args.output / "timing.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--chem-data", type=Path, required=True)
    parser.add_argument("--reference-python", type=Path, required=True)
    parser.add_argument("--native-python", type=Path, required=True)
    parser.add_argument("--native-capture", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--nuclear", action="store_true")
    parser.add_argument("--repeats", type=int, default=6)
    parser.add_argument("--timeout", type=int, default=120)
    args = parser.parse_args()
    if args.repeats < 2:
        parser.error("At least two repeats are required for AB/BA order")
    print(json.dumps(benchmark(args), indent=2))
