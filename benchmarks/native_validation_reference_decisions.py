"""Compare whole-entry restraint decisions and reported scalar values."""

from __future__ import annotations

import math
from collections import Counter
from collections.abc import Iterable
from typing import Any

from native_validation_reference_atoms import atom_key


def failure_key(
    check: str, atoms: Iterable[dict[str, Any]]
) -> tuple[str, tuple[tuple[str, ...], ...]]:
    keys = tuple(atom_key(atom) for atom in atoms)
    if check in ("bonds", "angles"):
        keys = min(keys, tuple(reversed(keys)))
    if check == "planarities":
        keys = tuple(sorted(keys))
    return check, keys


def expected_rows(
    reference: dict[str, Any], cutoff: float
) -> dict[tuple[str, tuple[tuple[str, ...], ...]], dict[str, Any]]:
    flags = {
        tuple(proxy["i_seqs"]): proxy.get("pseudochiral")
        for proxy in reference["graph"].get("chiralities", [])
    }
    atom_indices = {atom_key(atom): atom["index"] for atom in reference["atoms"]}
    rows = {}
    for check, category in reference.get("restraint_results", {}).items():
        for row in category["rows"]:
            if "atom_indices" in row:
                atoms = [reference["atoms"][index] for index in row["atom_indices"]]
            else:
                atoms = [
                    dict(
                        atom,
                        model=atom.get("model_id", ""),
                        chain=atom.get("chain_id", ""),
                        component=atom.get("resname", ""),
                    )
                    for atom in row["atoms_info"]
                ]
            if check == "chiralities":
                indices = tuple(atom_indices[atom_key(atom)] for atom in atoms)
                if flags.get(indices):
                    continue
            if row["score"] > cutoff:
                key = failure_key(check, atoms)
                if key in rows:
                    raise ValueError("Reference contains duplicate restraint findings")
                rows[key] = row
    return rows


def compare_decisions(
    reference: dict[str, Any], native: dict[str, Any], cutoff: float = 5.0
) -> dict[str, Any]:
    if "restraint_results" not in reference or "failures" not in native:
        return dict(decision_scalar_parity=False, unavailable=True)
    old = expected_rows(reference, cutoff)
    new = {
        failure_key(
            row["check"], [native["atoms"][index] for index in row["indices"]]
        ): row
        for row in native["failures"]
        if row["check"] in reference["restraint_results"]
    }
    native_count = sum(
        row["check"] in reference["restraint_results"] for row in native["failures"]
    )
    if len(new) != native_count:
        raise ValueError("Native report contains duplicate restraint findings")
    maximum_errors: Counter[str]
    mismatches, maximum_errors = [], Counter()
    for key in old.keys() & new.keys():
        fields = {
            "score": "score",
            "residual": "residual",
            "model": "model",
            "delta": "delta",
            "target": "target",
            "sigma": "sigma",
        }
        if key[0] == "planarities":
            fields = {
                "score": "score",
                "residual": "residual",
                "rms_deltas": "model",
                "delta_max": "delta",
            }
        for before, after in fields.items():
            if before not in old[key]:
                continue
            if after not in new[key]:
                mismatches.append(
                    dict(atoms=key, field=before, reason="Native scalar is missing")
                )
                continue
            a, b = old[key][before], new[key][after]
            maximum_errors[before] = max(maximum_errors[before], abs(a - b))
            if not math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-10):
                mismatches.append(dict(atoms=key, field=before, reference=a, native=b))
    missing, extra = list(old.keys() - new.keys()), list(new.keys() - old.keys())
    return dict(
        cutoff=cutoff,
        reference_outliers=len(old),
        native_outliers=len(new),
        reference_counts=dict(Counter(key[0] for key in old)),
        native_counts=dict(Counter(key[0] for key in new)),
        missing=missing,
        extra=extra,
        scalar_mismatches=mismatches,
        maximum_absolute_errors=dict(maximum_errors),
        decision_scalar_parity=not missing and not extra and not mismatches,
        scope="Strict score>5 decisions on reported required restraint categories; "
        "pseudochiral centers excluded. Plane score is maximum standardized "
        "atom deviation, not RMS. Comparison uses the captured processed sites.",
    )
