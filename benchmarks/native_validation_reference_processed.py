"""Freeze CCTBX-processed sites while retaining all original CIF categories."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any


def prepare(source: Path, reference: dict[str, Any], output: Path) -> dict[str, Any]:
    import iotbx.cif
    from scitbx.array_family import flex

    source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
    if reference["input"]["sha256"] != source_hash:
        raise ValueError("The source differs from the captured CCTBX input")
    document = iotbx.cif.reader(file_path=str(source)).model()
    block = next(block for block in document.values() if "_atom_site.Cartn_x" in block)

    def normalize(value: object) -> str:
        return "" if str(value).strip() in ("", ".", "?") else str(value).strip()

    keys = (
        "model",
        "chain",
        "resseq",
        "icode",
        "component",
        "altloc",
        "name",
        "element",
    )
    lookup = {
        tuple(normalize(atom[field]) for field in keys): atom
        for atom in reference["atoms"]
    }
    if len(lookup) != len(reference["atoms"]):
        raise ValueError("Processed atom identities are not unique")
    fields = (
        "pdbx_PDB_model_num",
        "auth_asym_id",
        "auth_seq_id",
        "pdbx_PDB_ins_code",
        "label_comp_id",
        "label_alt_id",
        "label_atom_id",
        "type_symbol",
    )
    columns = [list(block["_atom_site." + field]) for field in fields]
    if set(columns[0]) == {"1"}:
        columns[0] = [""] * len(columns[0])
    rows, used, changed = [], set(), []
    for index, values in enumerate(zip(*columns, strict=True)):
        key = tuple(normalize(value) for value in values)
        if key not in lookup or key in used:
            raise ValueError(
                f"Original CIF row has no unique processed atom match: {key}"
            )
        used.add(key)
        xyz = lookup[key]["xyz"]
        before = [
            float(block["_atom_site.Cartn_" + axis][index]) for axis in ("x", "y", "z")
        ]
        if before != xyz:
            changed.append(dict(atom=key, before=before, after=xyz))
        rows.append(xyz)
    if used != lookup.keys():
        raise ValueError(
            "CCTBX added or removed atoms; this site-only fixture cannot encode it"
        )
    for i, axis in enumerate(("x", "y", "z")):
        block["_atom_site.Cartn_" + axis] = flex.std_string(
            [format(row[i], ".17g") for row in rows]
        )
    with output.open("x") as stream:
        stream.write(str(document))
    return dict(
        source_sha256=source_hash,
        output_sha256=hashlib.sha256(output.read_bytes()).hexdigest(),
        atoms=len(rows),
        site_differences=changed,
        all_source_categories_retained=True,
        scope="The fixture uses recorded CCTBX sites to compare evaluators. "
        "It does not repair native coordinates or change the source structure. "
        "The original input and comparison remain available.",
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    args = parser.parse_args()
    result = prepare(args.input, json.loads(args.reference.read_text()), args.output)
    with args.manifest.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
