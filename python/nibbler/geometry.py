"""Evaluate an explicit restraint graph without repeating chemical interpretation.

The caller must supply complete, chemically assigned restraints with atom IDs
and flags in coordinate order. This API does not infer chemical topology,
expand crystallographic symmetry, or check Probe clashes. Unsupported topology
must use the reference validator; an empty category means no supplied proxies,
not proof that the structure needs none.

Bond rows contain ``(i, j, target, weight, slack, report)``. Angle rows contain
``(i, j, k, target_degrees, weight)``. Chirality rows contain
``(centre, i, j, k, target_volume, weight, both_signs, pseudochiral)``. Plane rows
contain ``(indices, weights)``. Weights equal inverse variance. All source
proxies count toward coverage, including non-reporting bonds and pseudochiral
centres. Hydrogen parent checks use every supplied bond.

Evaluation preserves the strict sigma cutoff. ``near_cutoff`` marks a score
within ``1e-8 * max(1, sigma_cutoff)`` on either side; callers that need the
reference decision must run the reference validator for those graphs.
"""

from ._core import RestraintGraph

__all__ = ["RestraintGraph"]
