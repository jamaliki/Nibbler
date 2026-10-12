"""Build metal and modified-nucleotide controls from declared source sites."""

import argparse
import hashlib
import json
from pathlib import Path

import iotbx.cif
import iotbx.pdb
import iotbx.pdb.hierarchy as hierarchy
import numpy as np
from cctbx import crystal
from native_validation_extension_fixtures import dictionary_sites, write
from scitbx.array_family import flex


def rotation(first, second):
    first, second = first / np.linalg.norm(first), second / np.linalg.norm(second)
    cross = np.cross(first, second)
    skew = np.array(
        [[0, -cross[2], cross[1]], [cross[2], 0, -cross[0]], [-cross[1], cross[0], 0]]
    )
    cosine = first @ second
    if cosine <= -0.999999:
        raise ValueError("The fixture rotation needs a non-antiparallel axis")
    return np.eye(3) + skew + skew @ skew / (1 + cosine)


def zinc(source, count):
    original = iotbx.pdb.input(file_name=str(source)).construct_hierarchy()
    template = next(group for group in original.atom_groups() if group.resname == "CYS")
    root, model = hierarchy.root(), hierarchy.model(id="")
    root.append_model(model)
    directions = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]) / np.sqrt(
        3
    )
    sg = np.array(template.get_atom("SG").xyz)
    cb = np.array(template.get_atom("CB").xyz)
    for index, direction in enumerate(directions[:count]):
        chain = hierarchy.chain(id=chr(65 + index))
        residue = hierarchy.residue_group(resseq="1", icode=" ")
        group = template.detached_copy()
        group.altloc = ""
        matrix = rotation(cb - sg, direction)
        for atom in group.atoms():
            atom.xyz = tuple(
                round(float(value), 3)
                for value in matrix @ (np.array(atom.xyz) - sg) + direction * 2.33
            )
        hg = group.get_atom("HG")
        if hg is not None:
            group.remove_atom(hg)
        residue.append_atom_group(group)
        chain.append_residue_group(residue)
        model.append_chain(chain)
    chain = hierarchy.chain(id="Z")
    residue = hierarchy.residue_group(resseq="1", icode=" ")
    group = hierarchy.atom_group(resname="ZN")
    atom = hierarchy.atom()
    atom.name, atom.element, atom.xyz = " ZN ", "ZN", (0, 0, 0)
    atom.occ, atom.b = 1, 10
    group.append_atom(atom)
    residue.append_atom_group(group)
    chain.append_residue_group(residue)
    model.append_chain(chain)
    for index, atom in enumerate(root.atoms(), 1):
        atom.serial = f"{index:5d}"
    return root


def renamed_nucleotide(source, dictionary, output):
    root = iotbx.pdb.input(file_name=str(source)).construct_hierarchy()
    for group in root.atom_groups():
        group.resname = "XRU"
    model = iotbx.cif.reader(file_path=str(dictionary)).model()
    renamed = iotbx.cif.model.cif()
    for name, block in model.items():
        for tag, value in list(block.items()):
            if isinstance(value, str):
                block[tag] = "XRU" if value == "U" else value
            else:
                block[tag] = flex.std_string(
                    ["XRU" if cell == "U" else cell for cell in value]
                )
        renamed[name.replace("comp_U", "comp_XRU")] = block
    with (output / "modified-rna-dictionary.cif").open("x") as stream:
        stream.write(str(renamed))
    return root


def symmetry_disulfide(source):
    root = zinc(source, 1)
    model = root.only_model()
    model.remove_chain(model.chains()[-1])
    first = model.only_chain()
    second = first.detached_copy()
    second.id = "B"
    first_sg = np.array(first.only_residue_group().only_atom_group().get_atom("SG").xyz)
    for chain, destination in (
        (first, np.array([1.0, 10, 10])),
        (second, np.array([18.969, 10, 10])),
    ):
        delta = destination - first_sg
        for atom in chain.atoms():
            atom.xyz = tuple(
                round(float(value), 3) for value in np.array(atom.xyz) + delta
            )
    model.append_chain(second)
    return root


def iron_sulfur(source):
    original = iotbx.pdb.input(file_name=str(source)).construct_hierarchy()
    clusters = [group for group in original.atom_groups() if group.resname == "SF4"]
    if len(clusters) != 1:
        raise ValueError("The control needs exactly one SF4 cluster")
    metals = [
        np.array(atom.xyz)
        for atom in clusters[0].atoms()
        if atom.element.strip().upper() == "FE"
    ]
    partners = [
        group
        for group in original.atom_groups()
        if group.resname == "CYS"
        and any(
            np.linalg.norm(np.array(group.get_atom("SG").xyz) - metal) < 3.5
            for metal in metals
        )
    ]
    root, model = hierarchy.root(), hierarchy.model(id="")
    root.append_model(model)
    for index, group in enumerate([clusters[0], *partners]):
        chain = hierarchy.chain(id=chr(65 + index))
        residue = hierarchy.residue_group(resseq="1", icode=" ")
        residue.append_atom_group(group.detached_copy())
        chain.append_residue_group(residue)
        model.append_chain(chain)
    return root


def glycosidic(chemistry):
    root = dictionary_sites(chemistry / "mon_lib/n/NAG.cif", "NAG")
    chain = root.only_model().only_chain()
    first = chain.only_residue_group().only_atom_group()
    second_residue = chain.only_residue_group().detached_copy()
    second_residue.resseq = "2"
    second = second_residue.only_atom_group()
    oxygen = np.array(first.get_atom("O4").xyz)
    carbon = np.array(first.get_atom("C4").xyz)
    c1, o1 = np.array(second.get_atom("C1").xyz), np.array(second.get_atom("O1").xyz)
    direction = (oxygen - carbon) / np.linalg.norm(oxygen - carbon)
    matrix = rotation(o1 - c1, -direction)
    for atom in second.atoms():
        atom.xyz = tuple(
            round(float(value), 3)
            for value in matrix @ (np.array(atom.xyz) - c1) + oxygen + direction * 1.439
        )
    for name in ("O1", "HO1"):
        second.remove_atom(second.get_atom(name))
    first.remove_atom(first.get_atom("HO4"))
    chain.append_residue_group(second_residue)
    return root


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protein", type=Path, required=True)
    parser.add_argument("--rna", type=Path, required=True)
    parser.add_argument("--rna-dictionary", type=Path, required=True)
    parser.add_argument("--iron-sulfur", type=Path)
    parser.add_argument("--chem-data", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    records = [
        write(
            args.output / f"zinc-cys{count}.cif",
            zinc(args.protein, count),
            {"scope": "Metal Coordination Library targets and omitted thiol H"},
        )
        for count in (1, 4)
    ]
    records.append(
        write(
            args.output / "modified-rna.cif",
            renamed_nucleotide(args.rna, args.rna_dictionary, args.output),
            {
                "scope": "A supplied RNA component uses the generic pucker mod"
            },
        )
    )
    symmetry = symmetry_disulfide(args.protein)
    file = args.output / "symmetry-disulfide.cif"
    symmetry_data = crystal.symmetry(
        unit_cell=(20, 40, 40, 90, 90, 90), space_group_symbol="P 1"
    )
    with file.open("x") as stream:
        stream.write(symmetry.as_mmcif_string(crystal_symmetry=symmetry_data))
    records.append(
        dict(
            file=file.name,
            sha256=hashlib.sha256(file.read_bytes()).hexdigest(),
            scope="Disulfide across a cell translation; source sites stay distant",
        )
    )
    if args.iron_sulfur:
        records.append(
            write(
                args.output / "iron-sulfur.cif",
                iron_sulfur(args.iron_sulfur),
                {
                    "scope": "Original SF4 and its four Cys ligands; no heme"
                },
            )
        )
    if args.chem_data:
        records.append(
            write(
                args.output / "glycosidic.cif",
                glycosidic(args.chem_data),
                {
                    "scope": "NAG C1-O4 link with the named leaving groups absent"
                },
            )
        )
    with (args.output / "manifest.json").open("x") as stream:
        json.dump(records, stream, indent=2)
