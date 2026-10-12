"""Native model validation and reusable explicit restraint evaluation.

For ``RestraintGraph``, supply chemically assigned restraints with atom IDs
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

__all__ = ["RestraintGraph", "validate"]


def validate(
    source, *, chem_data, nuclear=False, use_cdl=True, sigma_cutoff=5.0,
    include_graph=False,
):
    """Check existing coordinates with native restraints and Probe contacts.

    ``complete`` requires full chemical coverage. Unsupported checks and scores
    near a cutoff need reference validation. This API neither adds hydrogens nor
    decides whether training should exclude an entry. ``chem_data`` must remain
    unchanged while the process runs because native dictionaries are cached.

    ``include_graph`` returns atom identities and the independently built graph
    for reference comparisons. The native work releases the GIL.
    """
    import json
    from os import fspath

    from . import _core, cif
    from ._core import CifDocument, MmcifModel

    if not hasattr(_core, "validate_geometry_document"):
        raise RuntimeError("this Nibbler build lacks native model validation")
    for name, flag in (("nuclear", nuclear), ("use_cdl", use_cdl),
                       ("include_graph", include_graph)):
        if not isinstance(flag, bool):
            raise TypeError(f"{name} must be a bool")
    if isinstance(sigma_cutoff, bool) or not isinstance(sigma_cutoff, (int, float)):
        raise TypeError("sigma_cutoff must be a number")
    if isinstance(source, MmcifModel):
        document = source.source_document()
    elif isinstance(source, CifDocument):
        document = source
    else:
        document = cif.read(source)
    return json.loads(_core.validate_geometry_document(
        document, fspath(chem_data), nuclear, use_cdl, sigma_cutoff, include_graph,
    ))
