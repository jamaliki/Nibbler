"""Annotate a known symmetry bond without changing its source atom sites."""

import argparse
from pathlib import Path

import iotbx.cif


def build(source, output):
    document = iotbx.cif.reader(file_path=str(source)).model()
    block = next(block for block in document.values() if "_atom_site.id" in block)
    columns = {
        key: list(block[f"_atom_site.{key}"])
        for key in ("label_atom_id", "label_comp_id", "label_asym_id", "auth_seq_id")
    }
    indices = [columns["label_atom_id"].index(name) for name in ("SG", "ZN")]
    loop = iotbx.cif.model.loop()
    for field, value in (("id", "known_image"), ("conn_type_id", "metalc")):
        loop.add_column(f"_struct_conn.{field}", [value])
    for side, index in enumerate(indices, 1):
        for field in columns:
            loop.add_column(f"_struct_conn.ptnr{side}_{field}", [columns[field][index]])
        loop.add_column(
            f"_struct_conn.ptnr{side}_symmetry", ["1_555" if side == 1 else "1_655"]
        )
    block.add_loop(loop)
    with output.open("x") as stream:
        stream.write(str(document))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    build(args.input, args.output)
