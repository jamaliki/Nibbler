"""Compare declared chemistry extensions with independent CCTBX and Probe."""

import argparse
import hashlib
import json
from pathlib import Path

from native_validation_reference_capture import capture
from native_validation_reference_compare import compare

from nibbler.geometry import validate


def matrix(args):
    args.output.mkdir(parents=True, exist_ok=False)
    fixtures = args.fixtures
    cases = [
        (f"nucleotide-{base}", fixtures / f"nucleotide-{base}.cif", [], [])
        for base in ("A", "C", "G", "U", "DA", "DC", "DG", "DT", "A-pair")
    ]
    cases.extend(
        [
            ("swapped-val", fixtures / "crambin-swapped-val.cif", [], []),
            (
                "custom-ZZZ",
                fixtures / "custom-ZZZ-model.cif",
                [fixtures / "custom-ZZZ.cif"],
                [],
            ),
            (
                "enol-isolated",
                fixtures / "enol-isolated.cif",
                [fixtures / "enol-ALA.cif"],
                [],
            ),
            (
                "enol-linked",
                fixtures / "enol-linked.cif",
                [fixtures / "enol-THR.cif"],
                [],
            ),
            (
                "enol-override",
                fixtures / "enol-isolated.cif",
                [fixtures / "enol-target-override.cif"],
                [],
            ),
        ]
    )
    for name, folder in (
        ("cross-link", args.link),
        ("modified-link", args.modified_link),
    ):
        if folder:
            cases.append(
                (
                    name,
                    folder / "model.cif",
                    [folder / "dictionary.cif"],
                    [("XLINK", 0, 3)],
                )
            )
    records = []
    for name, source, dictionaries, links in cases:
        for nuclear in (False, True):
            label = f"{name}-{'nuclear' if nuclear else 'electron'}"
            reference = capture(
                source,
                args.chem_data,
                nuclear=nuclear,
                dictionaries=dictionaries,
                links=links,
            )
            native = validate(
                source,
                chem_data=args.chem_data,
                nuclear=nuclear,
                include_graph=True,
                restraint_dictionaries=dictionaries,
                restraint_links=links,
            )
            native["input"] = dict(
                file=str(source.resolve()),
                sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
            )
            native["restraint_dictionaries"] = [
                dict(
                    file=str(value.resolve()),
                    sha256=hashlib.sha256(value.read_bytes()).hexdigest(),
                )
                for value in dictionaries
            ]
            comparison = compare(reference, native)
            for kind, report in (
                ("reference", reference),
                ("native", native),
                ("comparison", comparison),
            ):
                with (args.output / f"{label}-{kind}.json").open("x") as stream:
                    json.dump(report, stream, indent=2, allow_nan=False)
            passed = comparison.get("all_required_parity", False)
            record = dict(
                name=label,
                passed=passed,
                reference_complete=reference["complete"],
                native_complete=native["complete"],
                error=reference.get("error"),
            )
            records.append(record)
            print(json.dumps(record), flush=True)
    with (args.output / "summary.json").open("x") as stream:
        json.dump(
            dict(
                cases=records,
                passed=sum(row["passed"] for row in records),
                total=len(records),
            ),
            stream,
            indent=2,
        )
    return all(row["passed"] for row in records)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixtures", type=Path, required=True)
    parser.add_argument("--chem-data", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--link", type=Path)
    parser.add_argument("--modified-link", type=Path)
    raise SystemExit(0 if matrix(parser.parse_args()) else 1)
