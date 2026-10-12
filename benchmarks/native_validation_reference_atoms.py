"""Normalize fixed-width atom labels for identity comparisons."""


def atom_key(atom):
    return tuple(
        str(atom.get(field, "")).strip()
        for field in (
            "model",
            "chain",
            "resseq",
            "icode",
            "component",
            "altloc",
            "name",
            "element",
        )
    )
