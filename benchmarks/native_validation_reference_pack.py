"""Replace repeated CCTBX atom metadata with immutable atom-table indices."""

import argparse
import hashlib
import json
from pathlib import Path

from native_validation_reference_atoms import atom_key


def pack(report):
    result = json.loads(json.dumps(report))
    indices = {atom_key(atom): atom["index"] for atom in report.get("atoms", [])}
    if len(indices) != len(report.get("atoms", [])):
        raise ValueError("The processed reference atom IDs are not unique")
    for category in result.get("restraint_results", {}).values():
        for row in category["rows"]:
            atoms = row.pop("atoms_info", [])
            row["atom_indices"] = []
            for atom in atoms:
                renamed = dict(
                    atom,
                    model=atom.get("model_id", ""),
                    chain=atom.get("chain_id", ""),
                    component=atom.get("resname", ""),
                )
                if atom.get("symop") not in (None, ""):
                    raise ValueError("Packing needs explicit symmetry-site metadata")
                key = atom_key(renamed)
                if key not in indices:
                    raise ValueError(
                        f"An evaluated atom lacks a processed atom-table match: {key}"
                    )
                row["atom_indices"].append(indices[key])
    result["restraint_result_atom_encoding"] = (
        "atom_indices reference the processed atoms table; "
    )
    result["restraint_result_atom_encoding"] += (
        "all CCTBX scalar and result fields remain unchanged."
    )
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = pack(json.loads(args.input.read_text()))
    report["full_capture_sha256"] = hashlib.sha256(args.input.read_bytes()).hexdigest()
    with args.output.open("x") as stream:
        json.dump(report, stream, separators=(",", ":"), allow_nan=False)
    print(
        json.dumps(
            dict(
                bytes_before=args.input.stat().st_size,
                bytes_after=args.output.stat().st_size,
            )
        )
    )
