"""Frozen user-visible workloads used to evaluate parser performance."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class Workload:
    """A benchmark operation and its required corpus features."""

    name: str
    operation: str
    required_features: tuple[str, ...]
    description: str


WORKLOADS = (
    Workload(
        name="full-document",
        operation="read",
        required_features=("scalar", "loop"),
        description="Parse a complete logical CIF document without projection.",
    ),
    Workload(
        name="projected-atom-site",
        operation="read",
        required_features=("polymer", "ligand", "ion", "water"),
        description="Read selected atom_site columns without materializing other data.",
    ),
    Workload(
        name="chemistry-heavy-validation",
        operation="validate",
        required_features=("modified-residue", "ligand", "ion", "water"),
        description="Validate PDBx entity, component, scheme, and atom relationships.",
    ),
    Workload(
        name="canonical-round-trip",
        operation="round-trip",
        required_features=("unknown", "not-applicable", "present-empty"),
        description="Parse, write canonically, and reparse without semantic loss.",
    ),
)


def workload_names() -> tuple[str, ...]:
    """Return workload names in their stable reporting order."""
    return tuple(workload.name for workload in WORKLOADS)
