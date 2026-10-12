"""Compare independently built native targets with processed CCTBX targets."""

import argparse
import gzip
import json
import math
from pathlib import Path

from native_validation_reference_atoms import atom_key


def load_report(file):
    opener = gzip.open if file.suffix == ".gz" else Path.open
    with opener(file, "rt") as stream:
        return json.load(stream)


def graph_rows(report, reference):
    atoms = [atom_key(row) for row in report["atoms"]]
    graph, result = report["graph"], {}
    for category in ("bonds", "angles", "chiralities", "planes"):
        rows = []
        for proxy in graph.get(
            "planarities" if reference and category == "planes" else category, []
        ):
            if reference:
                indices = proxy["i_seqs"]
                values = (
                    {
                        "bonds": lambda proxy=proxy: [
                            proxy["distance_ideal"],
                            proxy["weight"],
                            proxy.get("slack", 0),
                            proxy.get("origin_id", 0) == 0,
                        ],
                        "angles": lambda proxy=proxy: [
                            proxy["angle_ideal"],
                            proxy["weight"],
                        ],
                        "chiralities": lambda proxy=proxy: [
                            proxy["volume_ideal"],
                            proxy["weight"],
                            proxy["both_signs"],
                        ],
                        "planes": lambda proxy=proxy: list(proxy["weights"]),
                    }[category]
                )()
            else:
                arity = {"bonds": 2, "angles": 3, "chiralities": 4}.get(category)
                indices = proxy[0] if category == "planes" else proxy[:arity]
                values = proxy[1] if category == "planes" else proxy[arity:]
                if category == "chiralities":
                    values = values[:3]
            keys = tuple(atoms[index] for index in indices)
            if category in ("bonds", "angles"):
                keys = min(keys, tuple(reversed(keys)))
            if category == "planes":
                ordered = sorted(zip(keys, values, strict=True))
                keys, values = tuple(x[0] for x in ordered), [x[1] for x in ordered]
            rows.append((keys, list(values)))
        result[category] = rows
    return result


def compare_graph(reference, native):
    result = dict(
        reference_complete=reference.get("complete"),
        native_complete=native.get("complete"),
    )
    if not reference.get("complete") or "graph" not in native:
        return dict(
            result, comparable=False, errors=native.get("errors"), graph_parity=False
        )
    ref_atoms = {atom_key(row): row for row in reference["atoms"]}
    new_atoms = {atom_key(row): row for row in native["atoms"]}
    missing_atoms, extra_atoms = (
        sorted(ref_atoms.keys() - new_atoms.keys()),
        sorted(new_atoms.keys() - ref_atoms.keys()),
    )
    max_coordinate_error = max(
        (
            max(
                abs(x - y)
                for x, y in zip(
                    ref_atoms[key]["xyz"], new_atoms[key]["xyz"], strict=True
                )
            )
            for key in ref_atoms.keys() & new_atoms.keys()
        ),
        default=0,
    )
    result.update(
        comparable=True,
        missing_atoms=missing_atoms,
        extra_atoms=extra_atoms,
        max_coordinate_error=max_coordinate_error,
        categories={},
    )
    ref, new = graph_rows(reference, True), graph_rows(native, False)
    for category in ref:
        old, candidate = dict(ref[category]), dict(new[category])
        mismatch = []
        for key in old.keys() & candidate.keys():
            a, b = old[key], candidate[key]
            if len(a) != len(b) or any(
                (
                    x != y
                    if isinstance(x, bool) or isinstance(y, bool)
                    else not math.isclose(x, y, rel_tol=1e-9, abs_tol=1e-10)
                )
                for x, y in zip(a, b, strict=False)
            ):
                mismatch.append(dict(atoms=key, reference=a, native=b))
        result["categories"][category] = dict(
            reference_count=len(ref[category]),
            native_count=len(new[category]),
            duplicate_reference_keys=len(ref[category]) - len(old),
            duplicate_native_keys=len(new[category]) - len(candidate),
            missing=list(old.keys() - candidate.keys()),
            extra=list(candidate.keys() - old.keys()),
            target_weight_mismatches=mismatch,
        )
    result["unsupported_asu_proxy_count"] = len(reference["graph"].get("bonds_asu", []))
    result["graph_parity"] = (
        not missing_atoms
        and not extra_atoms
        and max_coordinate_error <= 1e-10
        and not result["unsupported_asu_proxy_count"]
        and all(
            not any(
                row[field]
                for field in (
                    "missing",
                    "extra",
                    "target_weight_mismatches",
                    "duplicate_reference_keys",
                    "duplicate_native_keys",
                )
            )
            for row in result["categories"].values()
        )
    )
    result["graph_target_weight_parity"] = result["graph_parity"]
    old_flags = {
        tuple(atom_key(reference["atoms"][i]) for i in proxy["i_seqs"]): proxy.get(
            "pseudochiral"
        )
        for proxy in reference["graph"].get("chiralities", [])
    }
    new_flags = {
        tuple(atom_key(native["atoms"][i]) for i in proxy[:4]): proxy[7]
        for proxy in native["graph"].get("chiralities", [])
    }
    result["missing_reference_pseudochiral_flags"] = sum(
        value is None for value in old_flags.values()
    )
    result["pseudochiral_flag_mismatches"] = [
        dict(atoms=key, reference=value, native=new_flags.get(key))
        for key, value in old_flags.items()
        if value is not None and value != new_flags.get(key)
    ]
    result["graph_parity"] &= (
        not result["missing_reference_pseudochiral_flags"]
        and not result["pseudochiral_flag_mismatches"]
    )
    return result


def compare_typing(reference, native):
    result = {}
    initial = reference.get("probe2", {}).get("initial_atom_info_calls", [])
    contact = native.get("contacts") or {}
    if initial and "initial_atom_info" in contact:
        old_info = {
            atom_key(initial[0]["atoms"][row["index"]]): row
            for row in initial[0]["extra_atom_info"]
        }
        new_info = {
            atom_key(native["atoms"][row["index"]]): row
            for row in contact["initial_atom_info"]
        }
        changed = []
        for key in old_info.keys() & new_info.keys():
            differences = {}
            for field in (
                "vdw_radius",
                "is_acceptor",
                "is_donor",
                "is_dummy_hydrogen",
                "is_ion",
                "charge",
                "alt",
            ):
                old_value, new_value = old_info[key][field], new_info[key][field]
                if field == "alt":
                    old_value = ord(old_value.strip()[:1]) if old_value.strip() else 0
                if old_value != new_value:
                    differences[field] = [old_value, new_value]
            if differences:
                changed.append(dict(atom=key, differences=differences))
        result["initial_probe_typing"] = dict(
            missing=list(old_info.keys() - new_info.keys()),
            extra=list(new_info.keys() - old_info.keys()),
            mismatches=changed,
        )
        result["initial_probe_typing_parity"] = not any(
            result["initial_probe_typing"].values()
        )
    else:
        result["initial_probe_typing_parity"] = False
        result["initial_probe_typing"] = dict(unavailable=True)
    return result


def compare_contacts(reference, native):
    result, probe, contact = {}, reference.get("probe2"), native.get("contacts") or {}
    if not probe:
        return dict(clash_pair_parity=False, clash_pairs=dict(unavailable=True))
    old_pairs = {}
    for row in probe["pairs"]:
        keys = tuple(
            tuple(
                str(row.get(prefix + field, "") or "").strip()
                for field in (
                    "model_id",
                    "chain_id",
                    "resseq",
                    "icode",
                    "resname",
                    "altloc",
                    "name",
                )
            )
            for prefix in ("", "target_")
        )
        old_pairs[tuple(sorted(keys))] = round(row["overlap"], 3)
    new_pairs = {}
    for row in contact.get("clashes", []):
        keys = tuple(atom_key(native["atoms"][index])[:-1] for index in row["indices"])
        new_pairs[tuple(sorted(keys))] = round(row["overlap"], 3)
    result["clash_pairs"] = dict(
        reference_count=len(old_pairs),
        native_count=len(new_pairs),
        missing=list(old_pairs.keys() - new_pairs.keys()),
        extra=list(new_pairs.keys() - old_pairs.keys()),
        overlap_mismatches=[
            dict(atoms=key, reference=old_pairs[key], native=new_pairs[key])
            for key in old_pairs.keys() & new_pairs.keys()
            if old_pairs[key] != new_pairs[key]
        ],
    )
    result["clashscore"] = dict(
        reference=probe["clashscore"], native=contact.get("clashscore")
    )
    result["clash_pair_parity"] = (
        not any(
            result["clash_pairs"][field]
            for field in ("missing", "extra", "overlap_mismatches")
        )
        and contact.get("clashscore") == probe["clashscore"]
    )
    return result


def compare(reference, native):
    result = compare_graph(reference, native)
    if not result["comparable"]:
        return result
    result.update(
        compare_typing(reference, native), **compare_contacts(reference, native)
    )
    from native_validation_reference_decisions import compare_decisions

    result["restraint_decisions"] = compare_decisions(reference, native)
    result["native_coverage"] = native.get("coverage", [])
    result["native_checks"] = native.get("checks", {})
    source = reference.get("input", {}).get("sha256")
    candidate_source = native.get("input", {}).get("sha256")
    result["source_bytes_match"] = (
        source == candidate_source if source and candidate_source else None
    )
    result["mode_match"] = all(
        reference.get(name) == native.get(name) for name in ("nuclear", "use_cdl")
    )
    result["dictionary_sources_match"] = (
        [row["sha256"] for row in reference.get("restraint_dictionaries", [])]
        == [row["sha256"] for row in native.get("restraint_dictionaries", [])]
    )
    result["requested_links_match"] = (
        [list(row) for row in reference.get("requested_links", [])]
        == [list(row) for row in native.get("requested_links", [])]
    )
    result["all_required_parity"] = bool(
        reference.get("complete")
        and native.get("complete")
        and result["graph_parity"]
        and result["initial_probe_typing_parity"]
        and result["clash_pair_parity"]
        and result["source_bytes_match"]
        and result["mode_match"]
        and result["dictionary_sources_match"]
        and result["requested_links_match"]
        and result["restraint_decisions"]["decision_scalar_parity"]
    )
    result["scope"] = (
        "Identity-mapped graph targets and weights; chirality order remains strict. "
    )
    result["scope"] += (
        "Initial Probe typing and rounded clash pairs use separate identity maps. "
    )
    result["scope"] += (
        "No untested contact-dot or coverage claim follows from these comparisons."
    )
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = compare(load_report(args.reference), load_report(args.native))
    with args.output.open("x") as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
    print(
        json.dumps(
            dict(
                graph_parity=report["graph_parity"],
                native_complete=report["native_complete"],
            )
        )
    )
