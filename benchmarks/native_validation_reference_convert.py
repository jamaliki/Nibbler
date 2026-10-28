"""Convert reference PDB inputs to mmCIF without chemical preparation."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any


def convert(source: Path, output: Path) -> dict[str, Any]:
    import iotbx.pdb

    hierarchy = iotbx.pdb.input(file_name=str(source)).construct_hierarchy()

    def identity(atom: Any) -> tuple[str, ...]:
        labels = atom.fetch_labels()
        return tuple(
            str(value).strip()
            for value in (
                labels.model_id,
                labels.chain_id,
                labels.resseq,
                labels.icode,
                labels.resname,
                labels.altloc,
                atom.name,
                atom.element,
            )
        )

    before = [(identity(atom), tuple(atom.xyz)) for atom in hierarchy.atoms()]
    text = hierarchy.as_mmcif_string()
    after = iotbx.pdb.input(
        source_info="reference conversion", lines=text
    ).construct_hierarchy()
    observed = [(identity(atom), tuple(atom.xyz)) for atom in after.atoms()]
    if before != observed:
        raise ValueError(
            "PDB-to-CIF conversion changed chemical atom IDs or coordinates"
        )
    with output.open("x") as stream:
        stream.write(text)
    report = dict(
        source=str(source),
        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        output=str(output),
        output_sha256=hashlib.sha256(output.read_bytes()).hexdigest(),
        atoms=len(before),
        coordinates_and_ids_unchanged=True,
        identity_scope="Model/chain/residue/component/altloc/name/element, stripping "
        "fixed-width label padding; coordinates must match exactly.",
        scope="Serialization only; no added H, restraint processing or repair.",
    )
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(convert(args.input, args.output), indent=2))
