"""Build independent MCL partner-selection controls from declared atom sites."""

import argparse
import hashlib
import json
from pathlib import Path

import iotbx.pdb
import numpy as np
from cctbx import crystal
from native_validation_extension_fixtures import dictionary_sites
from native_validation_remaining_fixtures import rotation, zinc


def place(group, anchor, neighbor, direction, distance):
    origin = np.array(group.get_atom(anchor).xyz)
    axis = np.array(group.get_atom(neighbor).xyz) - origin
    matrix = rotation(axis, direction)
    for atom in group.atoms():
        atom.xyz = tuple(
            round(float(value), 3)
            for value in matrix @ (np.array(atom.xyz) - origin)
            + np.array(direction) * distance
        )


def mixed_zinc(source, chemistry, cysteines, histidines, nitrogens):
    root = zinc(source, cysteines + histidines)
    model = root.only_model()
    template = dictionary_sites(chemistry / "mon_lib/h/HIS.cif", "HIS")
    histidine = (
        template.only_model().only_chain().only_residue_group().only_atom_group()
    )
    for index in range(cysteines, cysteines + histidines):
        residue = model.chains()[index].only_residue_group()
        original = residue.only_atom_group()
        direction = np.array(original.get_atom("SG").xyz)
        direction /= np.linalg.norm(direction)
        group = histidine.detached_copy()
        nitrogen = nitrogens[(index - cysteines) % len(nitrogens)]
        remove = "HD1" if nitrogen == "ND1" else "HE2"
        group.remove_atom(group.get_atom(remove))
        place(group, nitrogen, "CG" if nitrogen == "ND1" else "CD2", direction, 2.04)
        residue.remove_atom_group(original)
        residue.append_atom_group(group)
    return root


def write(root, output, name, scope, symmetry=None):
    for index, atom in enumerate(root.atoms(), 1):
        atom.serial = f"{index:5d}"
        atom.xyz = tuple(round(value, 3) for value in atom.xyz)
    file = output / f"{name}.cif"
    with file.open("x") as stream:
        stream.write(root.as_mmcif_string(crystal_symmetry=symmetry))
    return dict(
        file=file.name,
        sha256=hashlib.sha256(file.read_bytes()).hexdigest(),
        scope=scope,
    )


def build(args):
    args.output.mkdir(parents=True, exist_ok=False)
    records = []
    for cysteines, histidines, nitrogens in (
        (3, 1, ["ND1"]),
        (2, 2, ["ND1"]),
        (1, 3, ["NE2"]),
        (2, 2, ["ND1", "NE2"]),
        (0, 4, ["NE2"]),
        (0, 1, ["ND1"]),
    ):
        name = f"zinc-cys{cysteines}-his{histidines}-{'-'.join(nitrogens)}"
        records.append(
            write(
                mixed_zinc(
                    args.protein, args.chem_data, cysteines, histidines, nitrogens
                ),
                args.output,
                name,
                "MCL Cys/His counts and distinct nitrogen angle keys",
            )
        )
    for distance in (2.999, 3.001):
        root = zinc(args.protein, 1)
        group = root.only_model().chains()[0].only_residue_group().only_atom_group()
        sg = np.array(group.get_atom("SG").xyz)
        shift = sg / np.linalg.norm(sg) * distance - sg
        for atom in group.atoms():
            atom.xyz = tuple(
                round(float(value), 6) for value in np.array(atom.xyz) + shift
            )
        # The non-coordinated control must retain its thiol hydrogen.
        if distance > 3:
            import iotbx.pdb.hierarchy as hierarchy

            atom = hierarchy.atom()
            atom.name, atom.element, atom.occ, atom.b = " HG ", "H", 1, 10
            sg = np.array(group.get_atom("SG").xyz)
            cb = np.array(group.get_atom("CB").xyz)
            atom.xyz = tuple(sg + (sg - cb) / np.linalg.norm(sg - cb) * 1.2)
            group.append_atom(atom)
        records.append(
            write(
                root,
                args.output,
                f"zinc-cutoff-{distance}",
                "Zn neighbour distance on each side of the 3 A cutoff",
            )
        )
    root = mixed_zinc(args.protein, args.chem_data, 0, 1, ["ND1"])
    group = root.only_model().chains()[0].only_residue_group().only_atom_group()
    # Both nitrogens share an atom group. The reference selects only the first
    # nonbonded donor in its distance-minus-radius ordering.
    records.append(
        write(
            root,
            args.output,
            "zinc-competing-histidine",
            "Only one selected nitrogen per atom group",
        )
    )
    root = zinc(args.protein, 1)
    residue = root.only_model().chains()[0].only_residue_group()
    original = residue.only_atom_group()
    original.altloc = "A"
    duplicate = original.detached_copy()
    duplicate.altloc = "B"
    for atom in original.atoms():
        atom.occ = 0.5
    for atom in duplicate.atoms():
        atom.occ = 0.5
        atom.xyz = tuple(np.array(atom.xyz) + np.array([0.1, 0, 0]))
    residue.append_atom_group(duplicate)
    records.append(
        write(
            root,
            args.output,
            "zinc-altloc",
            "Separate atom-group donors and conformer exclusions",
        )
    )
    for name, cell in (
        ("zinc-symmetry", (20, 40, 40, 90, 90, 90)),
        ("zinc-oblique-symmetry", (20, 40, 40, 80, 75, 65)),
    ):
        root = zinc(args.protein, 1)
        symmetry = crystal.symmetry(unit_cell=cell, space_group_symbol="P 1")
        shift = np.array(symmetry.unit_cell().orthogonalize((1, 0, 0)))
        for atom in root.only_model().chains()[0].atoms():
            atom.xyz = tuple(np.array(atom.xyz) + shift)
        records.append(
            write(
                root,
                args.output,
                name,
                "Zn uses the nearest translated partner; source sites stay distant",
                symmetry,
            )
        )
    if args.iron_sulfur:
        root = iotbx.pdb.input(file_name=str(args.iron_sulfur)).construct_hierarchy()
        symmetry = crystal.symmetry(
            unit_cell=(40, 40, 40, 90, 90, 90), space_group_symbol="P 1"
        )
        for atom in root.only_model().chains()[1].atoms():
            atom.xyz = tuple(np.array(atom.xyz) + np.array([40, 0, 0]))
        import iotbx.pdb.hierarchy as hierarchy

        for group in root.atom_groups():
            if group.resname != "CYS":
                continue
            atom = hierarchy.atom()
            atom.name, atom.element, atom.occ, atom.b = " HG ", "H", 1, 10
            sg = np.array(group.get_atom("SG").xyz)
            cb = np.array(group.get_atom("CB").xyz)
            atom.xyz = tuple(sg + (sg - cb) / np.linalg.norm(sg - cb) * 1.2)
            group.append_atom(atom)
        records.append(
            write(
                root,
                args.output,
                "cluster-symmetry-reset",
                "Any cluster symmetry candidate clears the reference MCL cluster list",
                symmetry,
            )
        )
    with (args.output / "manifest.json").open("x") as stream:
        json.dump(records, stream, indent=2)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protein", type=Path, required=True)
    parser.add_argument("--chem-data", type=Path, required=True)
    parser.add_argument("--iron-sulfur", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    build(parser.parse_args())
