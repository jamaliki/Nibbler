"""Create explicit CB-HB1 bond cutoff controls without coordinate rounding."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any


def build(source: Path, folder: Path) -> list[dict[str, Any]]:
    import iotbx.pdb
    from scitbx.array_family import flex

    folder.mkdir(exist_ok=False)
    entries = []
    for offset in (-1e-7, 0.0, 1e-7):
        hierarchy = iotbx.pdb.input(file_name=str(source)).construct_hierarchy()
        atoms = {atom.name.strip(): atom for atom in hierarchy.atoms()}
        carbon, hydrogen = atoms["CB"], atoms["HB1"]
        direction = [h - c for h, c in zip(hydrogen.xyz, carbon.xyz, strict=True)]
        norm = math.sqrt(sum(value * value for value in direction))
        distance = 0.97 + 5 * 0.02 + offset
        hydrogen.xyz = tuple(
            c + distance * value / norm
            for c, value in zip(carbon.xyz, direction, strict=True)
        )
        block = hierarchy.as_cif_block()
        for i, axis in enumerate(("x", "y", "z")):
            block[f"_atom_site.Cartn_{axis}"] = flex.std_string(
                [format(atom.xyz[i], ".17g") for atom in hierarchy.atoms()]
            )
        file = folder / f"cb-hb1-{offset:+.0e}.cif"
        with file.open("x") as stream:
            stream.write("data_cutoff\n" + str(block))
        entries.append(
            dict(
                file=str(file),
                target=0.97,
                sigma=0.02,
                requested_distance=distance,
                offset=offset,
            )
        )
    return entries


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(build(args.input, args.output), indent=2))
