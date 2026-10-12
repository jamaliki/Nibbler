"""Build declared dictionary-site and label-swap fixtures for validation."""

import argparse
import hashlib
import json
from pathlib import Path

import iotbx.cif
import iotbx.pdb
import iotbx.pdb.hierarchy as hierarchy
from scitbx.array_family import flex
from scitbx.matrix import col


def dictionary_sites(file, component, allowed=None):
    model = iotbx.cif.reader(file_path=str(file)).model()
    block = model[f"comp_{component}"]
    root, chain = hierarchy.root(), hierarchy.chain(id="A")
    entry = hierarchy.model(id="")
    root.append_model(entry)
    entry.append_chain(chain)
    residue = hierarchy.residue_group(resseq="1", icode=" ")
    group = hierarchy.atom_group(altloc="", resname=component)
    chain.append_residue_group(residue)
    residue.append_atom_group(group)
    columns = [
        block[f"_chem_comp_atom.{name}"]
        for name in ("atom_id", "type_symbol", "x", "y", "z")
    ]
    for index, (name, element, x, y, z) in enumerate(zip(*columns, strict=True), 1):
        if component == "DT":
            name = {"C7": "C5M", "H71": "H5M1", "H72": "H5M2", "H73": "H5M3"}.get(
                name, name
            )
        if allowed is not None and name not in allowed:
            continue
        atom = hierarchy.atom()
        atom.serial = f"{index:5d}"
        atom.name, atom.element = f"{name:>4}", f"{element:>2}"
        atom.xyz, atom.occ, atom.b = (float(x), float(y), float(z)), 1.0, 10.0
        group.append_atom(atom)
    return root


def write(output, root, provenance):
    if output.exists():
        raise FileExistsError(output)
    output.write_text(root.as_mmcif_string())
    return dict(
        file=output.name,
        sha256=hashlib.sha256(output.read_bytes()).hexdigest(),
        atoms=root.atoms_size(),
        **provenance,
    )


def custom_component(args):
    file = args.chem_data / "geostd/a/data_ALA.cif"
    original = iotbx.cif.reader(file_path=str(file)).model()
    modified = iotbx.cif.model.cif()
    for name, block in original.items():
        for tag, values in list(block.items()):
            if isinstance(values, str):
                block[tag] = "ZZZ" if values == "ALA" else values
            else:
                block[tag] = flex.std_string(
                    ["ZZZ" if value == "ALA" else value for value in values]
                )
        modified[name.replace("comp_ALA", "comp_ZZZ")] = block
    dictionary = args.output / "custom-ZZZ.cif"
    if dictionary.exists():
        raise FileExistsError(dictionary)
    dictionary.write_text(str(modified))
    root = iotbx.pdb.input(file_name=str(args.alanine)).construct_hierarchy()
    for group in root.atom_groups():
        group.resname = "ZZZ"
    return write(
        args.output / "custom-ZZZ-model.cif",
        root,
        dict(
            dictionary_sha256=hashlib.sha256(dictionary.read_bytes()).hexdigest(),
            source_sha256=hashlib.sha256(args.alanine.read_bytes()).hexdigest(),
            scope="Rename the validated alanine control and its full dictionary",
        ),
    )


def enol_fixture(source, output, linked=False):
    root = iotbx.pdb.input(file_name=str(source)).construct_hierarchy()
    groups = list(root.atom_groups())
    group = groups[0]
    if linked:
        group.resname = "ENT"
    oxygen, carbon = group.get_atom("O"), group.get_atom("C")
    atom = hierarchy.atom()
    atom.serial, atom.name, atom.element = f"{root.atoms_size() + 1:5d}", " HNO", " H"
    atom.xyz = tuple(
        round(value, 3)
        for value in col(oxygen.xyz)
        + (col(oxygen.xyz) - col(carbon.xyz)).normalize() * 0.97
    )
    atom.occ, atom.b = 1.0, 10.0
    group.append_atom(atom)
    if linked:
        hydrogen = groups[1].get_atom("H")
        if hydrogen is None:
            raise ValueError("The linked enol fixture needs a removable N-H")
        groups[1].remove_atom(hydrogen)
    return write(
        output,
        root,
        dict(
            source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
            scope="Add HNO; remove the following N-H only for a linked enol",
        ),
    )


def enol_dictionary(root, component, output):
    model = iotbx.cif.reader(
        file_path=str(root / component[0].lower() / f"data_{component}.cif")
    ).model()
    block = next(
        block for block in model.values() if "_chem_comp_atom.atom_id" in block
    )

    def append(category, values):
        loop = block.get_loop(f"_{category}")
        loop.add_row([values.get(tag.split(".")[1], ".") for tag in loop.keys()])

    append(
        "chem_comp_atom",
        dict(comp_id=component, atom_id="HNO", type_symbol="H", type_energy="HOH1"),
    )
    append(
        "chem_comp_bond",
        dict(
            comp_id=component,
            atom_id_1="HNO",
            atom_id_2="O",
            type="single",
            value_dist=".97",
            value_dist_neutron=".97",
            value_dist_esd=".02",
        ),
    )
    for index, (a, b) in enumerate(
        zip(
            block["_chem_comp_bond.atom_id_1"],
            block["_chem_comp_bond.atom_id_2"],
            strict=True,
        )
    ):
        if {a, b} == {"C", "O"}:
            for tag, value in (
                ("value_dist", "1.43"),
                ("value_dist_neutron", "1.43"),
                ("value_dist_esd", ".02"),
            ):
                block[f"_chem_comp_bond.{tag}"][index] = value
    append(
        "chem_comp_angle",
        dict(
            comp_id=component,
            atom_id_1="HNO",
            atom_id_2="O",
            atom_id_3="C",
            value_angle="120",
            value_angle_esd="3",
        ),
    )
    if component == "THR":
        renamed = iotbx.cif.model.cif()
        for name, block in model.items():
            for tag, values in list(block.items()):
                if isinstance(values, str):
                    block[tag] = "ENT" if values == "THR" else values
                else:
                    block[tag] = flex.std_string(
                        ["ENT" if value == "THR" else value for value in values]
                    )
            renamed[name.replace("THR", "ENT")] = block
        model = renamed
    with output.open("x") as stream:
        stream.write(str(model))


def nucleotide_pair(source, output):
    root = iotbx.pdb.input(file_name=str(source)).construct_hierarchy()
    chain = root.only_model().only_chain()
    first = chain.only_residue_group().only_atom_group()
    c1, outbound, oxygen = [
        col(first.get_atom(name).xyz) for name in ("C1'", "N9", "O3'")
    ]
    axis = (outbound - c1).normalize()
    radial = oxygen - c1 - axis * (oxygen - c1).dot(axis)
    destination = oxygen - radial.normalize() * 1.6
    residue = chain.only_residue_group().detached_copy()
    residue.resseq = "2"
    shift = destination - col(residue.only_atom_group().get_atom("P").xyz)
    for atom in residue.atoms():
        atom.xyz = tuple(round(value, 3) for value in col(atom.xyz) + shift)
    chain.append_residue_group(residue)
    for index, atom in enumerate(root.atoms(), 1):
        atom.serial = f"{index:5d}"
    return write(
        output,
        root,
        dict(
            source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
            scope="The next phosphate tests the RNA pair's 2-prime pucker selection",
        ),
    )


def conflicting_enol_dictionary(source, output):
    model = iotbx.cif.reader(file_path=str(source)).model()
    block = model["comp_ALA"]
    for row, (a, b) in enumerate(
        zip(
            block["_chem_comp_bond.atom_id_1"],
            block["_chem_comp_bond.atom_id_2"],
            strict=True,
        )
    ):
        if {a, b} in ({"C", "O"}, {"HNO", "O"}):
            for column in ("value_dist", "value_dist_neutron"):
                block[f"_chem_comp_bond.{column}"][row] = "1.2"
    for row, atom in enumerate(block["_chem_comp_angle.atom_id_1"]):
        if atom == "HNO":
            block["_chem_comp_angle.value_angle"][row] = "90"
    with output.open("x") as stream:
        stream.write(str(model))


def fixtures(args):
    args.output.mkdir(parents=True, exist_ok=True)
    records = []
    if args.alanine:
        records.append(custom_component(args))
        records.append(enol_fixture(args.alanine, args.output / "enol-isolated.cif"))
        records.append(
            enol_fixture(args.crambin, args.output / "enol-linked.cif", linked=True)
        )
        enol_dictionary(args.chem_data / "geostd", "ALA", args.output / "enol-ALA.cif")
        enol_dictionary(args.chem_data / "geostd", "THR", args.output / "enol-THR.cif")
        conflicting_enol_dictionary(
            args.output / "enol-ALA.cif", args.output / "enol-target-override.cif"
        )
    for component in ("A", "C", "G", "U", "DA", "DC", "DG", "DT"):
        file = args.chem_data / "mon_lib" / component[0].lower() / f"{component}.cif"
        if not file.exists():
            continue
        lookup = {"DA": "AD", "DC": "CD", "DG": "GD", "DT": "TD"}.get(
            component, component
        )
        effective = args.chem_data / "geostd" / lookup[0].lower() / f"data_{lookup}.cif"
        allowed = None
        if effective.exists():
            dictionary = iotbx.cif.reader(file_path=str(effective)).model()
            definitions = [
                block
                for block in dictionary.values()
                if "_chem_comp_atom.atom_id" in block
            ]
            if len(definitions) != 1:
                raise ValueError("Expected one effective component dictionary")
            allowed = set(definitions[0]["_chem_comp_atom.atom_id"])
        records.append(
            write(
                args.output / f"nucleotide-{component}.cif",
                dictionary_sites(file, component, allowed),
                dict(
                    dictionary=str(file),
                    dictionary_sha256=hashlib.sha256(file.read_bytes()).hexdigest(),
                    effective_dictionary_sha256=(
                        hashlib.sha256(effective.read_bytes()).hexdigest()
                        if effective.exists()
                        else None
                    ),
                    scope=(
                        "Declared monomer sites restricted to the effective "
                        "dictionary's atom set; no hydrogen preparation"
                    ),
                ),
            )
        )
    root = iotbx.pdb.input(file_name=str(args.crambin)).construct_hierarchy()
    records.append(
        nucleotide_pair(
            args.output / "nucleotide-A.cif", args.output / "nucleotide-A-pair.cif"
        )
    )
    changes = []
    for group in root.atom_groups():
        if group.resname != "VAL":
            continue
        for left, right in (
            ("CG1", "CG2"),
            ("HG11", "HG21"),
            ("HG12", "HG22"),
            ("HG13", "HG23"),
        ):
            a, b = group.get_atom(left), group.get_atom(right)
            if a is None or b is None:
                raise ValueError("VAL fixture lacks a swap pair")
            a.xyz, b.xyz = b.xyz, a.xyz
            a.b, b.b = 12.0, 34.0
            changes.append([a.id_str(), b.id_str()])
        break
    records.append(
        write(
            args.output / "crambin-swapped-val.cif",
            root,
            dict(
                source_sha256=hashlib.sha256(args.crambin.read_bytes()).hexdigest(),
                deliberate_swaps=changes,
            ),
        )
    )
    manifest = args.output / "manifest.json"
    if manifest.exists():
        raise FileExistsError(manifest)
    manifest.write_text(json.dumps(records, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chem-data", type=Path, required=True)
    parser.add_argument("--crambin", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--alanine", type=Path)
    fixtures(parser.parse_args())
