"""Pack ligand tokens from the table the caller passes, as the element list did before it."""

import json
import os
from pathlib import Path

import numpy as np
import pytest

import nibbler
from nibbler.observations import ObservationPacker
from token_tables import frozen, packed, token_table

ROOT = Path(__file__).parent
RECORD = json.loads(frozen("token_table_parity.json.gz"))
EXTENDED = "SI PD TI SB TE RE CR GE SN TL RH".split()
CHEM_CONFIGS = {"chem_data", "molecular"}
LIGAND = ["C", *EXTENDED, "ZR"]


def packer(table):
    result = ObservationPacker(
        frozen("atom32_schema.json.gz"), table, RECORD["excluded"]
    )
    result.set_molecular_schema(frozen("molecular_schema.json.gz"))
    return result


@pytest.mark.parametrize("case", sorted(RECORD["cases"]))
def test_the_66_token_table_reproduces_the_element_list_packer(case):
    # Recorded with Nibbler 9b8041e from the 43 elements and hydrogen; see record_token_parity.py.
    # Only molecular crop_positions digests were corrected for the no-bridge packing bug.
    chem_data = os.environ.get("REDUCE3_CHEM_DATA")
    expected = RECORD["cases"][case]
    actual = packed(
        packer(token_table(*RECORD["ligand_elements"])),
        (ROOT / case).read_bytes(),
        chem_data,
    )
    if chem_data is None:
        expected = {
            name: value for name, value in expected.items() if name not in CHEM_CONFIGS
        }
    assert actual.keys() == expected.keys()
    for name in expected:
        assert actual[name] == expected[name], name


def test_the_parity_record_packs_ligands_metals_and_the_molecular_path():
    outcomes = [value for case in RECORD["cases"].values() for value in case.values()]
    assert sum("error" not in value for value in outcomes) >= 70
    packed_cases = {
        case
        for case, value in RECORD["cases"].items()
        if "error" not in value["molecular"]
    }
    assert {"fixtures/chemistry/branched_glycan.cif", "reduce3/mbo.cif"} <= packed_cases


def document():
    rows = "".join(
        f"HETATM {element}1 {element} G LIG 1 . 101 . {3 * index} 0 0 1\n"
        for index, element in enumerate(LIGAND)
    )
    return nibbler.cif.read(
        (
            "data_ligand\nloop_\n_atom_site.group_PDB\n_atom_site.label_atom_id\n"
            "_atom_site.type_symbol\n_atom_site.label_asym_id\n_atom_site.label_comp_id\n"
            "_atom_site.label_entity_id\n_atom_site.label_seq_id\n_atom_site.auth_seq_id\n"
            "_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n"
            "_atom_site.pdbx_PDB_model_num\n" + rows
        ).encode()
    )


def aatype(table):
    buffers, keys = packer(table).pack(document(), validate=True)
    dtype, shape, buffer = buffers["aatype"]
    return np.frombuffer(buffer, dtype=dtype).reshape(shape).tolist(), [
        key[4] for key in keys
    ]


def test_the_77_token_table_gives_the_eleven_elements_66_to_76():
    elements = RECORD["ligand_elements"]
    extended = token_table(*elements)
    extended["elements"] |= {
        element: 66 + index for index, element in enumerate(EXTENDED)
    }
    tokens, names = aatype(extended)
    assert names == [f"{element}1" for element in LIGAND]
    assert tokens == [21, *range(66, 77), 64]
    assert aatype(token_table(*elements))[0] == [21, *[64] * len(EXTENDED), 64]


@pytest.mark.parametrize(
    "table,message",
    [
        (
            dict(elements={"C": 21, "O": 21}, other=64, hydrogen=65),
            "given to both C and O",
        ),
        (
            dict(elements={"C": 64}, other=64, hydrogen=65),
            "given to both C and other elements",
        ),
        (
            dict(elements={"C": 65}, other=64, hydrogen=65),
            "given to both C and hydrogen",
        ),
        (
            dict(elements={"C": 20}, other=64, hydrogen=65),
            "not above the residue tokens",
        ),
        (
            dict(elements={"C": 21}, other=0, hydrogen=65),
            "not above the residue tokens",
        ),
        (dict(elements={"H": 66}, other=64, hydrogen=65), "Hydrogen takes"),
        (dict(elements={"Si": 66}, other=64, hydrogen=65), "not an upper-case symbol"),
        (dict(elements={"C": 21}, other=64), "no hydrogen token"),
        (dict(elements={"C": 21}, hydrogen=65), "no other token"),
        (dict(other=64, hydrogen=65), "no elements"),
        (
            dict(elements={"C": 21}, other=64, hydrogen=65, unknown=66),
            "Unknown ligand token table field",
        ),
    ],
)
def test_tables_that_break_the_token_contract_are_refused(table, message):
    with pytest.raises(ValueError, match=message):
        packer(table)


def test_an_element_list_is_not_a_table():
    with pytest.raises(TypeError):
        packer([*RECORD["ligand_elements"], "H"])
