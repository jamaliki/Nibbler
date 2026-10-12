"""Build a declared cross-component chirality control for both validators."""

import argparse
from pathlib import Path

import iotbx.cif.model
import iotbx.pdb.hierarchy as hierarchy


def table(block, category, columns, rows):
    block.add_loop(
        iotbx.cif.model.loop(
            header=[f"_{category}.{column}" for column in columns],
            data=[value for row in rows for value in row],
        )
    )


def build(output, modification=False):
    output.mkdir(parents=True, exist_ok=False)
    root = hierarchy.root()
    model = hierarchy.model(id="")
    root.append_model(model)
    definitions = iotbx.cif.model.cif()
    for component, chain_id, atoms, bonds, angles in [
        (
            "XAA",
            "A",
            [
                ("C1", "C", (0, 0, 0)),
                ("C2", "C", (1.5, 0, 0)),
                ("C3", "C", (0, 1.5, 0)),
            ],
            [
                ["XAA", "C1", "C2", "single", "1.5", ".02"],
                ["XAA", "C1", "C3", "single", "1.5", ".02"],
            ],
            [["XAA", "C2", "C1", "C3", "90", "3"]],
        ),
        (
            "XBB",
            "B",
            [("C4", "C", (0, 0, 1.5)), ("H4", "H", (0, 0, 2.47))],
            [["XBB", "C4", "H4", "single", ".97", ".02"]],
            [],
        ),
    ]:
        chain = hierarchy.chain(id=chain_id)
        residue = hierarchy.residue_group(resseq="1", icode=" ")
        group = hierarchy.atom_group(resname=component, altloc="")
        model.append_chain(chain)
        chain.append_residue_group(residue)
        residue.append_atom_group(group)
        for name, element, xyz in atoms:
            atom = hierarchy.atom()
            atom.serial = f"{root.atoms_size() + 1:5d}"
            atom.name, atom.element = f"{name:>4}", f"{element:>2}"
            atom.xyz, atom.occ, atom.b, atom.hetero = xyz, 1.0, 10.0, True
            group.append_atom(atom)
        block = iotbx.cif.model.block()
        table(block, "chem_comp", ["id", "group"], [[component, "non-polymer"]])
        table(
            block,
            "chem_comp_atom",
            ["comp_id", "atom_id", "type_symbol", "type_energy"],
            [
                [component, name, element, "HCH1" if element == "H" else "C"]
                for name, element, _ in atoms
            ],
        )
        table(
            block,
            "chem_comp_bond",
            [
                "comp_id",
                "atom_id_1",
                "atom_id_2",
                "type",
                "value_dist",
                "value_dist_esd",
            ],
            bonds,
        )
        if angles:
            table(
                block,
                "chem_comp_angle",
                [
                    "comp_id",
                    "atom_id_1",
                    "atom_id_2",
                    "atom_id_3",
                    "value_angle",
                    "value_angle_esd",
                ],
                angles,
            )
        definitions[f"comp_{component}"] = block
    block = iotbx.cif.model.block()
    table(
        block,
        "chem_link",
        [
            "id",
            "comp_id_1",
            "mod_id_1",
            "group_comp_1",
            "comp_id_2",
            "mod_id_2",
            "group_comp_2",
            "name",
        ],
        [
            [
                "XLINK",
                "XAA",
                ".",
                ".",
                "XBB",
                "XBB-H" if modification else ".",
                ".",
                "XLINK",
            ]
        ],
    )
    definitions["link_list"] = block
    block = iotbx.cif.model.block()
    table(
        block,
        "chem_link_bond",
        [
            "link_id",
            "atom_1_comp_id",
            "atom_id_1",
            "atom_2_comp_id",
            "atom_id_2",
            "value_dist",
            "value_dist_esd",
        ],
        [["XLINK", "1", "C1", "2", "C4", "1.5", ".02"]],
    )
    table(
        block,
        "chem_link_angle",
        [
            "link_id",
            "atom_1_comp_id",
            "atom_id_1",
            "atom_2_comp_id",
            "atom_id_2",
            "atom_3_comp_id",
            "atom_id_3",
            "value_angle",
            "value_angle_esd",
        ],
        [
            ["XLINK", "1", name, "1", "C1", "2", "C4", "90", "3"]
            for name in ("C2", "C3")
        ],
    )
    table(
        block,
        "chem_link_chir",
        [
            "link_id",
            "atom_centre_comp_id",
            "atom_id_centre",
            "atom_1_comp_id",
            "atom_id_1",
            "atom_2_comp_id",
            "atom_id_2",
            "atom_3_comp_id",
            "atom_id_3",
            "volume_sign",
        ],
        [["XLINK", "1", "C1", "1", "C2", "1", "C3", "2", "C4", "positive"]],
    )
    definitions["link_XLINK"] = block
    if modification:
        block = iotbx.cif.model.block()
        table(
            block,
            "chem_mod_atom",
            ["mod_id", "function", "atom_id", "new_type_energy"],
            [["XBB-H", "change", "H4", "HCH3"]],
        )
        table(
            block,
            "chem_mod_bond",
            [
                "mod_id",
                "function",
                "atom_id_1",
                "atom_id_2",
                "new_value_dist",
                "new_value_dist_esd",
            ],
            [["XBB-H", "change", "C4", "H4", ".95", ".01"]],
        )
        definitions["mod_XBB-H"] = block
    with (output / "dictionary.cif").open("x") as stream:
        stream.write(str(definitions))
    structure = root.as_cif_block()
    table(
        structure,
        "struct_conn",
        [
            "id",
            "conn_type_id",
            "ptnr1_auth_asym_id",
            "ptnr1_auth_seq_id",
            "ptnr1_label_comp_id",
            "ptnr1_label_atom_id",
            "ptnr2_auth_asym_id",
            "ptnr2_auth_seq_id",
            "ptnr2_label_comp_id",
            "ptnr2_label_atom_id",
            "ptnr1_symmetry",
            "ptnr2_symmetry",
        ],
        [
            [
                "XLINK",
                "covale",
                "A",
                "1",
                "XAA",
                "C1",
                "B",
                "1",
                "XBB",
                "C4",
                "1_555",
                "1_555",
            ]
        ],
    )
    cif = iotbx.cif.model.cif()
    cif["model"] = structure
    with (output / "model.cif").open("x") as stream:
        stream.write(str(cif))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--modification", action="store_true")
    args = parser.parse_args()
    build(args.output, args.modification)
