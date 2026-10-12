"""Build complete FES and F3S controls with declared Cys partner sites."""

import argparse
from pathlib import Path

import numpy as np
from native_validation_extension_fixtures import dictionary_sites
from native_validation_metal_fixtures import place, write
from native_validation_remaining_fixtures import zinc


def cluster(source, chemistry, name):
    root = dictionary_sites(chemistry / f"mon_lib/f/{name}.cif", name)
    model = root.only_model()
    model.only_chain().id = "Z"
    group = model.only_chain().only_residue_group().only_atom_group()
    center = np.mean([atom.xyz for atom in group.atoms()], axis=0)
    metals = [atom for atom in group.atoms() if atom.element.strip().upper() == "FE"]
    donors = zinc(source, 4).only_model().chains()[:4]
    for index, metal in enumerate(metals):
        chain = donors[index].detached_copy()
        donor = chain.only_residue_group().only_atom_group()
        direction = np.array(metal.xyz) - center
        direction /= np.linalg.norm(direction)
        place(donor, "SG", "CB", direction, 2.305 if name == "FES" else 2.318)
        for atom in chain.atoms():
            atom.xyz = tuple(np.array(atom.xyz) + np.array(metal.xyz))
        model.append_chain(chain)
    return root


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protein", type=Path, required=True)
    parser.add_argument("--chem-data", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    for name in ("FES", "F3S"):
        write(
            cluster(args.protein, args.chem_data, name),
            args.output,
            name,
            "MCL cluster site naming and weighted external bonds/angles",
        )
