"""Retain a direct CCTBX plane control for its signed-wrapper discrepancy."""

import argparse
import hashlib
import json
import math
from pathlib import Path

from cctbx import geometry_restraints
from native_validation_reference_atoms import atom_key
from scitbx.array_family import flex


def build(source, output):
    reference = json.loads(source.read_text())
    lookup = {atom_key(atom): atom["index"] for atom in reference["atoms"]}
    wrappers = {}
    for row in reference["restraint_results"]["planarities"]["rows"]:
        atoms = [
            dict(
                atom,
                model=atom.get("model_id", ""),
                chain=atom.get("chain_id", ""),
                component=atom.get("resname", ""),
            )
            for atom in row["atoms_info"]
        ]
        indices = tuple(sorted(lookup[atom_key(atom)] for atom in atoms))
        wrappers[indices] = row
    sites = flex.vec3_double([atom["xyz"] for atom in reference["atoms"]])
    choices = []
    for proxy in reference["graph"]["planarities"]:
        indices, weights = proxy["i_seqs"], proxy["weights"]
        wrapper = wrappers.get(tuple(sorted(indices)))
        if wrapper is None:
            continue
        restraint = geometry_restraints.planarity(
            sites_cart=sites,
            proxy=geometry_restraints.planarity_proxy(
                flex.size_t(indices), flex.double(weights)
            ),
        )
        deltas = list(restraint.deltas())
        score = max(abs(d) * math.sqrt(w) for d, w in zip(deltas, weights, strict=True))
        if score - wrapper["score"] > 1:
            choices.append(
                dict(
                    sites=[list(sites[index]) for index in indices],
                    weights=weights,
                    score=score,
                    signed_wrapper_score=wrapper["score"],
                    deltas=deltas,
                    residual=restraint.residual(),
                    rms_deltas=restraint.rms_deltas(),
                    atoms=[reference["atoms"][index] for index in indices],
                )
            )
    if not choices:
        raise ValueError("The capture has no signed-plane wrapper discrepancy")
    control = max(choices, key=lambda row: row["score"] - row["signed_wrapper_score"])
    control.update(
        capture_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        input_sha256=reference["input"]["sha256"],
        primitive="cctbx.geometry_restraints.planarity",
        scope="Absolute standardized deviation; signed wrapper kept separately",
    )
    with output.open("x") as stream:
        json.dump(control, stream, indent=2, allow_nan=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    build(args.reference, args.output)
