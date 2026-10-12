"""Pack complete extension captures beside their unchanged CIF fixtures."""

import argparse
import gzip
import hashlib
import json
from pathlib import Path

from native_validation_reference_pack import pack


def freeze(args):
    records = []
    for file in sorted(args.matrix.glob("*-reference.json")):
        name = file.name.removesuffix("-reference.json")
        family = name.removesuffix("-electron").removesuffix("-nuclear")
        reference = json.loads(file.read_text())
        comparison = json.loads(file.with_name(f"{name}-comparison.json").read_text())
        if not reference["complete"] or not comparison.get("all_required_parity"):
            raise ValueError(f"Cannot freeze incomplete comparison {name}")
        source = args.output / (
            f"{family}.cif"
            if family in ("cross-link", "modified-link")
            else Path(reference["input"]["file"]).name
        )
        if (
            hashlib.sha256(source.read_bytes()).hexdigest()
            != reference["input"]["sha256"]
        ):
            raise ValueError(f"The copied source changed for {name}")
        dictionaries = []
        for row in reference.get("restraint_dictionaries", []):
            dictionary = (
                f"{family}-dictionary.cif"
                if family in ("cross-link", "modified-link")
                else Path(row["file"]).name
            )
            if (
                hashlib.sha256((args.output / dictionary).read_bytes()).hexdigest()
                != row["sha256"]
            ):
                raise ValueError(f"The copied dictionary changed for {name}")
            dictionaries.append(dict(file=dictionary, sha256=row["sha256"]))
        target = args.output / f"{name}-reference.json.gz"
        with gzip.open(target, "xt") as stream:
            json.dump(pack(reference), stream, separators=(",", ":"), allow_nan=False)
        records.append(
            dict(
                name=name,
                input=source.name,
                input_sha256=reference["input"]["sha256"],
                reference=target.name,
                reference_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
                dictionaries=dictionaries,
                links=reference.get("requested_links", []),
                nuclear=reference["nuclear"],
                use_cdl=reference["use_cdl"],
            )
        )
    with (args.output / "manifest.json").open("x") as stream:
        json.dump(
            dict(
                cases=records, reference="CCTBX 2026.9 and Probe2", count=len(records)
            ),
            stream,
            indent=2,
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    freeze(parser.parse_args())
