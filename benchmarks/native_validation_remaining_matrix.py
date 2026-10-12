"""Compare the remaining chemistry paths with unmodified CCTBX and Probe."""

import argparse
import hashlib
import json
from pathlib import Path

from native_validation_reference_capture import capture
from native_validation_reference_compare import compare

from nibbler.geometry import validate


def run(args):
    args.output.mkdir(parents=True, exist_ok=False)
    cases = [
        (f"zinc-cys{count}", args.fixtures / f"zinc-cys{count}.cif", [])
        for count in (1, 4)
    ]
    cases.append(
        (
            "modified-rna",
            args.fixtures / "modified-rna.cif",
            [args.fixtures / "modified-rna-dictionary.cif"],
        )
    )
    cases.append(("symmetry-disulfide", args.fixtures / "symmetry-disulfide.cif", []))
    cases.extend((file.stem, file, []) for file in args.extra)
    results = []
    for name, source, dictionaries in cases:
        for nuclear in (False, True):
            label = f"{name}-{'nuclear' if nuclear else 'electron'}"
            reference = capture(
                source, args.chem_data, nuclear=nuclear, dictionaries=dictionaries
            )
            native = validate(
                source,
                chem_data=args.chem_data,
                nuclear=nuclear,
                include_graph=True,
                restraint_dictionaries=dictionaries,
            )
            native["input"] = dict(
                sha256=hashlib.sha256(source.read_bytes()).hexdigest()
            )
            native["restraint_dictionaries"] = [
                dict(sha256=hashlib.sha256(file.read_bytes()).hexdigest())
                for file in dictionaries
            ]
            comparison = compare(reference, native)
            for kind, report in (
                ("reference", reference),
                ("native", native),
                ("comparison", comparison),
            ):
                with (args.output / f"{label}-{kind}.json").open("x") as stream:
                    json.dump(report, stream, indent=2, allow_nan=False)
            row = dict(
                name=label,
                passed=comparison.get("all_required_parity", False),
                reference_complete=reference["complete"],
                native_complete=native["complete"],
                error=reference.get("error"),
            )
            print(json.dumps(row), flush=True)
            results.append(row)
    with (args.output / "summary.json").open("x") as stream:
        json.dump(results, stream, indent=2)
    return all(row["passed"] for row in results)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixtures", type=Path, required=True)
    parser.add_argument("--chem-data", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--extra", type=Path, action="append", default=[])
    raise SystemExit(0 if run(parser.parse_args()) else 1)
