"""Normalize fixed-width atom labels for identity comparisons."""

from __future__ import annotations

from typing import Any


def atom_key(atom: dict[str, Any]) -> tuple[str, ...]:
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
