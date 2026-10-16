"""Conventional PDBx/mmCIF and ModelCIF semantic API."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Literal

from . import _core, cif
from ._core import CifDocument, MmcifModel
from ._core import ConformerSelection as ConformerSelection
from ._core import PreparedChemistry as PreparedChemistry
from ._input import Source
from ._native import (
    raise_chemistry_error,
    raise_schema_error,
    raise_write_error,
    validation_report,
)
from .cif import Destination
from .components import Registry
from .contracts import Profile, ValidationReport


def ligand_sasa(
    source: Source | CifDocument | MmcifModel,
    *,
    atom_rows: Sequence[int] | None = None,
    model_number: int = 1,
    probe_radius: float = 1.4,
    point_count: int = 256,
    radii: Mapping[str, float] | None = None,
) -> dict[int, float]:
    """Return ligand SASA in A^2 keyed by zero-based source atom_site row.

    Only ligand atoms receive sphere samples. All occupied non-water heavy atoms
    in the chosen model can occlude them. Automatic selection includes non-polymer
    and branched entities, including ions; explicit atom_rows overrides selection.
    Select alternate locations with select_conformer before calling this function.
    Assemblies and symmetry contacts count only when present as explicit atoms.

    Uses elemental Bondi radii and a Fibonacci Shrake-Rupley sphere. radii overrides
    element symbols; unknown elements raise ValueError. H, D, T and water do not
    count. This heavy-atom convention does not infer missing H or use ProtOr radii.
    Existing documents need no second parse, CCD lookup, or model construction.
    The native calculation releases the GIL. Sum the values for total ligand SASA.
    """
    if isinstance(model_number, bool) or not isinstance(model_number, int):
        raise TypeError("model_number must be an integer")
    if isinstance(point_count, bool) or not isinstance(point_count, int):
        raise TypeError("point_count must be an integer")
    overrides = (
        None if radii is None else {key.upper(): value for key, value in radii.items()}
    )
    if isinstance(source, MmcifModel):
        rows, areas = _core.ligand_sasa_model(
            source, atom_rows, model_number, probe_radius, point_count, overrides
        )
    else:
        document = source if isinstance(source, CifDocument) else cif.read(source)
        rows, areas = _core.ligand_sasa_document(
            document, atom_rows, model_number, probe_radius, point_count, overrides
        )
    return dict(zip(rows, areas, strict=True))


def read(
    source: Source | CifDocument,
    *,
    profile: Profile | str = Profile.PDBX,
    registry: Registry | None = None,
) -> MmcifModel:
    """Build an immutable PDBx or ModelCIF semantic model."""
    selected = Profile(profile)
    document = (
        source
        if isinstance(source, CifDocument)
        else cif.read(source, schema=selected.value)
    )
    try:
        return _core.build_mmcif_model(document, selected.value, registry)
    except ValueError as error:
        raise_chemistry_error(error)


def prepare_chemistry(
    source: Source | CifDocument | MmcifModel,
    *,
    registry: Registry | None = None,
    model_number: int = 1,
    alternate_id: str | None = None,
    chem_data: str | None = None,
    use_observed_hydrogens: bool = True,
) -> PreparedChemistry:
    """Resolve authored covalent products without changing the source model.

    This separate operation returns an immutable native graph. ``to_dict()``
    exports its atoms, bond orders, hydrogen counts and repair records on demand.
    Missing template sites have no coordinates; inactive leaving atoms retain
    their source row. Metal coordination does not count toward covalent valence.
    ``model_number`` is the deposited number. Without ``alternate_id``, selection
    uses the highest summed occupancy per residue; selection is not clash-aware.
    A supplied alternate takes precedence where present. Unknown linked chemistry
    raises :class:`nibbler.ChemistryError`. This does not assign a pH or place atoms.
    ``chem_data`` supplies cached CCD templates and cannot accompany ``registry``.
    Set ``use_observed_hydrogens=False`` to derive expected H counts from the
    templates and product rules, independent of H atoms already in the source.
    """
    if isinstance(model_number, bool) or not isinstance(model_number, int):
        raise TypeError("model_number must be an integer")
    if alternate_id is not None and not isinstance(alternate_id, str):
        raise TypeError("alternate_id must be a string")
    model = source if isinstance(source, MmcifModel) else read(source)
    try:
        return _core.prepare_mmcif_chemistry(
            model, registry, model_number, alternate_id, chem_data, use_observed_hydrogens
        )
    except ValueError as error:
        raise_chemistry_error(error)


def select_conformer(
    source: Source | CifDocument | MmcifModel,
    *,
    model_number: int | None = None,
    seed: int = 0,
) -> ConformerSelection:
    """Select clash-free residue alternatives without changing chemical identity.

    Selection prefers higher occupancy and uses the seed to order ties. Nearby
    residues on different chains prefer matching letters on equal occupancy.
    The native search uses sparse spatial contacts and a 100,000-choice limit;
    an unresolved conflict raises ``ChemistryError`` rather than returning clashes.
    ``atom_rows`` indexes the source table. ``document`` creates the selected view
    only when requested. This never substitutes a modified residue's parent.
    """
    if model_number is not None and (
        isinstance(model_number, bool) or not isinstance(model_number, int)
    ):
        raise TypeError("model_number must be an integer")
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise TypeError("seed must be a nonnegative integer")
    try:
        if isinstance(source, MmcifModel):
            return _core.select_mmcif_conformer(source, model_number, seed)
        document = source if isinstance(source, CifDocument) else cif.read(source)
        return _core.select_document_conformer(document, model_number, seed)
    except ValueError as error:
        raise_chemistry_error(error)


def select_atom_conformers(
    coordinates: list[list[float]],
    occupancy: list[float],
    residue_ids: list[int],
    chain_ids: list[str],
    sequence_ids: list[int | None],
    atom_names: list[str],
    alternate_ids: list[str | None],
    *,
    bonds: list[list[int]],
    clash_mask: list[bool],
    seed: int = 0,
) -> ConformerSelection:
    """Select already parsed atom arrays without a second CIF read.

    Residue IDs name occurrences; bond indices refer to input rows. ``clash_mask``
    must exclude H, D, and T from contact tests; their occupancies still count.
    Both selection APIs apply the same rules. The array result has no source document.
    """
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise TypeError("seed must be a nonnegative integer")
    try:
        return _core.select_conformer_arrays(
            coordinates,
            occupancy,
            residue_ids,
            chain_ids,
            sequence_ids,
            atom_names,
            alternate_ids,
            bonds,
            clash_mask,
            seed,
        )
    except ValueError as error:
        raise_chemistry_error(error)


def deduplicate_assembly_copies(
    positions: list[list[float]],
    source_atoms: list[int],
    units: list[int],
    copies: list[int],
    eligible: list[bool],
    occupied: list[bool],
) -> tuple[list[int], list[list[int]]]:
    """Remove whole same-source copies with at least 50% overlap within 1 A.

    Lower copy IDs take precedence; rank identity operators first. Mark attached
    multi-atom units ineligible unless their mapping is exact. The result gives
    retained input rows and removed-to-retained source atom pairs. These arrays
    let a caller remap its edges and keep preparation provenance.
    """
    try:
        return _core.deduplicate_assembly_arrays(
            positions, source_atoms, units, copies, eligible, occupied
        )
    except ValueError as error:
        raise_chemistry_error(error)


def validate(
    value: object,
    *,
    profile: Profile | str,
) -> ValidationReport:
    """Validate a document or semantic model against an explicit profile."""
    selected = Profile(profile)
    return _validate_value(value, selected)


def assembly_with_mapping(
    source: Source | CifDocument | MmcifModel,
    assembly_id: str | None = None,
    *,
    registry: Registry | None = None,
) -> tuple[CifDocument, list[tuple[int, str, str]]]:
    """Expand once and return removed-to-retained whole-unit overlap mappings.

    Each triple gives a row in the input atom_site table, the removed label chain,
    and the retained label chain. This does not report the older partial same-atom
    coincidence rule. No second expansion or CIF parse runs to obtain the mapping.
    """
    if assembly_id is not None and not isinstance(assembly_id, str):
        raise TypeError("assembly_id must be a string")
    model = (
        source if isinstance(source, MmcifModel) else read(source, registry=registry)
    )
    try:
        return _core.build_assembly_with_mapping(model, assembly_id)
    except ValueError as error:
        raise_chemistry_error(error)


def _validate_value(value: object, profile: Profile) -> ValidationReport:
    try:
        if isinstance(value, MmcifModel):
            if value.profile != profile.value:
                raise ValueError(
                    f"model profile is {value.profile!r}, not {profile.value!r}"
                )
            fields = _core.validate_mmcif_model(value)
        elif isinstance(value, CifDocument):
            fields = _core.validate_mmcif_document(value, profile.value)
        else:
            raise TypeError(
                "nibbler.mmcif.validate() requires CifDocument or MmcifModel"
            )
    except ValueError as error:
        raise_schema_error(error)
    return validation_report(fields)


def write(
    value: object,
    destination: Destination,
    *,
    profile: Profile | str,
    mode: Literal["canonical", "preserve"] = "canonical",
    validate: Literal["none", "syntax", "dictionary", "profile"] = "profile",
    mirror_local_qa_metric: int | None = None,
    format: Literal["cif", "bcif"] | None = None,
) -> None:
    """Write one complete semantic model with transactional path output."""
    selected = Profile(profile)
    if not isinstance(value, MmcifModel):
        raise TypeError("nibbler.mmcif.write() requires MmcifModel")
    if value.profile != selected.value:
        raise ValueError(f"model profile is {value.profile!r}, not {selected.value!r}")
    if mirror_local_qa_metric is not None:
        if selected is not Profile.MODELCIF:
            raise ValueError("mirror_local_qa_metric requires profile='modelcif'")
        if isinstance(mirror_local_qa_metric, bool) or not isinstance(
            mirror_local_qa_metric, int
        ):
            raise TypeError("mirror_local_qa_metric must be an integer metric ID")
    if mode not in {"canonical", "preserve"}:
        raise ValueError("mode must be 'canonical' or 'preserve'")
    if mode == "preserve" and mirror_local_qa_metric is not None:
        raise ValueError("mode='preserve' cannot mirror a QA metric into B factors")
    if validate not in {"none", "syntax", "dictionary", "profile"}:
        raise ValueError(
            "validate must be 'none', 'syntax', 'dictionary', or 'profile'"
        )
    if validate == "profile":
        _validate_value(value, selected).raise_for_errors()

    generic_validation: Literal["none", "syntax", "dictionary"]
    if validate == "profile":
        generic_validation = "dictionary"
    else:
        generic_validation = validate
    try:
        document = (
            value.source_document()
            if mode == "preserve"
            else _core._model_document(value, mirror_local_qa_metric)
        )
    except ValueError as error:
        raise_write_error(error)
    cif.write(
        document,
        destination,
        mode=mode,
        validate=generic_validation,
        format=format,
    )


def assembly(
    source: Source | CifDocument | MmcifModel,
    assembly_id: str | None = None,
    *,
    registry: Registry | None = None,
) -> CifDocument:
    """Write out every copy of one biological assembly as explicit coordinates.

    ``assembly_id`` selects a ``_pdbx_struct_assembly`` (default: the first). Each
    distinct operator combination of its ``_pdbx_struct_assembly_gen`` rows is one copy;
    copy 1 keeps the chain identifiers and copy *n* appends ``-n`` to the
    ``label_asym_id`` and ``auth_asym_id`` of its chains. An atom on a symmetry axis of
    the assembly, whose copies coincide, is written once. The result holds the entry,
    entity, and component categories, and one row per copy in ``_struct_asym``, the
    sequence schemes, ``_atom_site`` (coordinates transformed, atom ids renumbered),
    ``_atom_site_anisotrop`` (tensors rotated), and ``_struct_conn`` (connections
    within one crystal copy). The crystal cell, symmetry, and assembly definitions are
    left out, so nothing downstream applies the crystal symmetry again.
    """
    if assembly_id is not None and not isinstance(assembly_id, str):
        raise TypeError("assembly_id must be a string")
    model = (
        source
        if isinstance(source, MmcifModel)
        else read(source, profile=Profile.PDBX, registry=registry)
    )
    try:
        return _core.build_assembly_document(model, assembly_id)
    except ValueError as error:
        raise_chemistry_error(error)
