"""Tests for the versioned fixture corpus."""

from __future__ import annotations

from benchmarks.workloads import workload_names
from tools.fetch_schemas import load_schema_locks
from tools.verify_corpus import load_manifest, verify_corpus


def test_fixture_manifest_is_complete_and_immutable() -> None:
    assert verify_corpus() == ()


def test_corpus_prioritizes_chemistry_risks() -> None:
    entries = load_manifest()["fixture"]
    features = {feature for entry in entries for feature in entry["features"]}
    assert {"ligand", "ion", "water", "modified-residue", "glycan"} <= features


def test_benchmark_workload_names_are_stable() -> None:
    assert workload_names() == (
        "full-document",
        "projected-atom-site",
        "chemistry-heavy-validation",
        "canonical-round-trip",
    )


def test_schema_versions_and_digests_are_locked() -> None:
    locks = load_schema_locks()
    versions = {entry["name"]: entry["version"] for entry in locks["schema"]}
    assert versions == {"pdbx": "5.416", "modelcif": "1.4.9"}
    assert all(len(entry["sha256"]) == 64 for entry in locks["schema"])
