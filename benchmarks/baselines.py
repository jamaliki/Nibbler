"""Reference engines used for correctness-first benchmark comparisons."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class Baseline:
    """A reference engine and the distribution that supplies it."""

    name: str
    distribution: str
    purpose: str


BASELINES = (
    Baseline("gemmi", "gemmi", "fast general CIF parsing"),
    Baseline("biotite", "biotite", "columnar structural-biology workflows"),
    Baseline("biopython", "biopython", "Bio.PDB compatibility workloads"),
)
