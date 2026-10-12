"""Capture the independently built native graph for the reference harness."""

import argparse
import hashlib
import json
import time
from pathlib import Path

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--chem-data", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--nuclear", action="store_true")
    parser.add_argument("--no-cdl", action="store_true")
    parser.add_argument("--dictionary", type=Path, action="append", default=[])
    args = parser.parse_args()
    started = time.perf_counter()
    import nibbler._core
    import nibbler.geometry

    report = nibbler.geometry.validate(
        args.input.read_bytes(),
        chem_data=args.chem_data,
        nuclear=args.nuclear,
        use_cdl=not args.no_cdl,
        include_graph=True,
        restraint_dictionaries=args.dictionary,
    )
    report.update(
        elapsed_seconds=time.perf_counter() - started,
        input=dict(
            file=str(args.input.resolve()),
            sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),
        ),
        native_binary=str(nibbler._core.__file__),
        native_binary_sha256=hashlib.sha256(
            Path(nibbler._core.__file__).read_bytes()
        ).hexdigest(),
    )
    with args.output.open("x") as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
    print(
        json.dumps(
            dict(complete=report["complete"], elapsed_seconds=report["elapsed_seconds"])
        )
    )
