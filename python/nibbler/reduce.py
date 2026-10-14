"""Hydrogen addition and optimization with Reduce3, on documents in memory."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
import json
from os import PathLike, fspath
from typing import Literal

from . import _core, cif
from ._core import CifDocument, MmcifModel
from ._input import Source
from ._native import raise_chemistry_error

Approach = Literal["add", "remove", "optimize"]
NTerminalCharge = Literal["residue_one", "first_in_chain", "no_charge"]

_APPROACHES = frozenset({"add", "remove", "optimize"})
_N_TERMINAL_CHARGES = frozenset({"residue_one", "first_in_chain", "no_charge"})
_PROBE_NUMBERS = frozenset(
    {
        "probe_radius",
        "density",
        "worse_clash_cutoff",
        "clash_cutoff",
        "contact_cutoff",
        "uncharged_hydrogen_cutoff",
        "charged_hydrogen_cutoff",
        "bump_weight",
        "hydrogen_bond_weight",
        "gap_weight",
    }
)
_PROBE_FLAGS = frozenset(
    {
        "allow_weak_hydrogen_bonds",
        "ignore_ion_interactions",
        "set_polar_hydrogen_radius",
    }
)


@dataclass(frozen=True, slots=True)
class Reduction:
    """The output model of one Reduce3 run and its report."""

    document: CifDocument
    report: str
    target_uncertainties: tuple[dict, ...] = ()


def _require_number(name: str, value: object) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{name} must be a number")
    return float(value)


def available() -> bool:
    """Return whether this Nibbler build includes Reduce3."""
    return hasattr(_core, "reduce_document")


def hydrogen_parents(
    document: CifDocument, *, chem_data: str | PathLike[str] | None = None,
    model_id: int = 1,
) -> dict[str, dict[str, str]]:
    """Return hydrogen parents from the dictionaries Reduce3 selects."""
    if not isinstance(document, CifDocument):
        raise TypeError("hydrogen_parents requires a parsed CIF document")
    if isinstance(model_id, bool) or not isinstance(model_id, int) or model_id < 1:
        raise ValueError("model_id must be a one-based model number")
    if not hasattr(_core, "reduce_hydrogen_parents"):
        raise RuntimeError("This Nibbler build has no Reduce3 hydrogen-parent helper")
    return _core.reduce_hydrogen_parents(document,
        chem_data=None if chem_data is None else fspath(chem_data), model=model_id)


def required_parent_sources(
    document: CifDocument, expected_keys_json: str, *, model_id: int = 1,
    chem_data: str | PathLike[str] | None = None,
) -> list[int]:
    """Map checked target keys to original heavy-atom source rows."""
    if not isinstance(document, CifDocument):
        raise TypeError("required_parent_sources requires a parsed CIF document")
    if not isinstance(expected_keys_json, str):
        raise TypeError("expected_keys_json must be a JSON string")
    if isinstance(model_id, bool) or not isinstance(model_id, int) or model_id < 1:
        raise ValueError("model_id must be a one-based model number")
    if not hasattr(_core, "reduce_required_parent_sources"):
        raise RuntimeError("This Nibbler build has no required-parent source helper")
    return _core.reduce_required_parent_sources(document, expected_keys_json, model=model_id,
        chem_data=None if chem_data is None else fspath(chem_data))


def run(
    source: Source | CifDocument | MmcifModel,
    *,
    approach: Approach = "add",
    add_flip_movers: bool = True,
    compat: bool = False,
    chem_data: str | PathLike[str] | None = None,
    n_terminal_charge: NTerminalCharge = "residue_one",
    keep_existing_h: bool = False,
    exclude_water: bool = True,
    use_neutron_distances: bool = False,
    preference_magnitude: float = 1.0,
    non_flip_preference: float = 0.5,
    planar_hydroxyl_preference: float = 1.0,
    acid_syn_preference: float = 1.0,
    skip_bond_fix_up: bool = False,
    set_flip_states: str | None = None,
    model_id: int | None = None,
    alt_id: str | None = None,
    bonded_neighbor_depth: int = 4,
    stop_on_any_missing_hydrogen: bool = False,
    ignore_missing_restraints: bool = False,
    allow_missing_restraints_for: Sequence[str] = (),
    required_hydrogen_target_keys: str | None = None,
    allow_undefined_ligand_targets: bool = False,
    allow_unresolved_bond_targets: bool = False,
    verbosity: int = 2,
    probe: Mapping[str, float | bool] | None = None,
) -> Reduction:
    """Add hydrogens and optimize movable groups with Reduce3, in memory.

    Reduce3 reads the first data block with an ``_atom_site`` loop straight from the
    parsed document and builds the result as a new :class:`CifDocument`; no text is
    written or parsed again. The options are Reduce2's parameters with the same
    defaults, except that ``add_flip_movers`` is True (Reduce2: False), so Asn, Gln
    and His flips are considered unless it is set to False. ``compat=True`` reproduces
    Reduce2 exactly; pass ``add_flip_movers=False`` too to match a default Reduce2
    run. Without ``compat``, a hydroxyl hydrogen on a planar atom (a phenol, an enol,
    a carboxylic acid) prefers that atom's plane by ``planar_hydroxyl_preference``,
    and an acid's hydrogen prefers syn to its carbonyl oxygen by
    ``acid_syn_preference`` (Probe score units; 0 turns either off). ``chem_data`` is
    the cctbx monomer library directory; by default it is found through
    ``REDUCE3_CHEM_DATA``, ``CHEM_DATA``, or the active conda environment. The run
    releases the GIL.

    ``allow_missing_restraints_for`` keeps named components as contact context
    when they lack a dictionary. It adds no hydrogens to those components and
    records them in the report. Known dictionaries and strict checks for all
    other components and unparameterized hydrogens stay unchanged. This option
    requires ``compat=False``.

    ``required_hydrogen_target_keys`` contains the retained target-key JSON.
    Native code resolves its heavy parent rows in the selected model. Strict
    fixed-mode addition may remove and report an unplaced H only when its
    unique parent lies outside those targets. Valid context hydrogens and all
    heavy atoms remain unchanged. This option requires explicit ``model_id``.

    ``allow_undefined_ligand_targets`` retains incomplete ligand fragments and
    reports only undefined H targets. It also retains validated bridging
    hydrides as observed atoms or reports their missing coordinates.
    ``allow_unresolved_bond_targets`` reports target uncertainty for deposited
    bonds that the placement graph cannot support. Both options require strict
    fixed-mode addition and retained keys. ``target_uncertainties`` names the
    affected original heavy atoms; pass it to ObservationPacker when packing
    training targets. Ordinary protein H placement remains strict.
    """
    if not available():
        raise RuntimeError(
            "this Nibbler build does not include Reduce3; rebuild it with the "
            "'reduce3' Cargo feature"
        )
    if isinstance(source, MmcifModel):
        document = source.source_document()
    elif isinstance(source, CifDocument):
        document = source
    else:
        document = cif.read(source)
    if approach not in _APPROACHES:
        raise ValueError("approach must be 'add', 'remove', or 'optimize'")
    if n_terminal_charge not in _N_TERMINAL_CHARGES:
        raise ValueError(
            "n_terminal_charge must be 'residue_one', 'first_in_chain', or 'no_charge'"
        )
    if (isinstance(allow_missing_restraints_for, (str, bytes))
            or not isinstance(allow_missing_restraints_for, Sequence)
            or any(not isinstance(name, str) for name in allow_missing_restraints_for)):
        raise TypeError("allow_missing_restraints_for must be a sequence of component identifiers")
    if required_hydrogen_target_keys is not None and not isinstance(required_hydrogen_target_keys, str):
        raise TypeError("required_hydrogen_target_keys must be a JSON string")
    flags = {
        "add_flip_movers": add_flip_movers,
        "compat": compat,
        "keep_existing_h": keep_existing_h,
        "exclude_water": exclude_water,
        "use_neutron_distances": use_neutron_distances,
        "skip_bond_fix_up": skip_bond_fix_up,
        "stop_on_any_missing_hydrogen": stop_on_any_missing_hydrogen,
        "ignore_missing_restraints": ignore_missing_restraints,
        "allow_undefined_ligand_targets": allow_undefined_ligand_targets,
        "allow_unresolved_bond_targets": allow_unresolved_bond_targets,
    }
    for name, flag in flags.items():
        if not isinstance(flag, bool):
            raise TypeError(f"{name} must be a bool")
    if model_id is not None and (
        isinstance(model_id, bool) or not isinstance(model_id, int) or model_id < 1
    ):
        raise ValueError("model_id must be a one-based model number")
    if (
        isinstance(bonded_neighbor_depth, bool)
        or not isinstance(bonded_neighbor_depth, int)
        or bonded_neighbor_depth < 0
    ):
        raise ValueError("bonded_neighbor_depth must be a non-negative integer")
    if isinstance(verbosity, bool) or not isinstance(verbosity, int):
        raise TypeError("verbosity must be an integer")
    for name, number in (
        ("preference_magnitude", preference_magnitude),
        ("non_flip_preference", non_flip_preference),
        ("planar_hydroxyl_preference", planar_hydroxyl_preference),
        ("acid_syn_preference", acid_syn_preference),
    ):
        _require_number(name, number)
    probe_numbers: dict[str, float] = {}
    probe_flags: dict[str, bool] = {}
    for name, value in (probe or {}).items():
        if name in _PROBE_FLAGS:
            if not isinstance(value, bool):
                raise TypeError(f"probe[{name!r}] must be a bool")
            probe_flags[name] = value
        elif name in _PROBE_NUMBERS:
            probe_numbers[name] = _require_number(f"probe[{name!r}]", value)
        else:
            raise ValueError(f"unknown probe parameter {name!r}")
    options = {
        **flags,
        "chem_data": None if chem_data is None else fspath(chem_data),
        "approach": approach,
        "n_terminal_charge": n_terminal_charge,
        "preference_magnitude": float(preference_magnitude),
        "non_flip_preference": float(non_flip_preference),
        "planar_hydroxyl_preference": float(planar_hydroxyl_preference),
        "acid_syn_preference": float(acid_syn_preference),
        "set_flip_states": set_flip_states,
        "model_id": model_id,
        "alt_id": alt_id,
        "bonded_neighbor_depth": bonded_neighbor_depth,
        "verbosity": verbosity,
        "probe": probe_numbers,
        "probe_flags": probe_flags,
        "allow_missing_restraints_for": list(allow_missing_restraints_for),
        "required_hydrogen_target_keys": required_hydrogen_target_keys,
    }
    try:
        output, report, statuses = _core.reduce_document(document, options)
    except ValueError as error:
        raise_chemistry_error(error)
    return Reduction(document=output, report=report, target_uncertainties=tuple(json.loads(statuses)))
