# Native link validation

The link extension retains automatic-link choices from native interpretation
and emits weighted targets in Nibbler. It adds automatic dictionary links,
glycosidic chirality, CCTBX metal defaults, zinc coordination targets,
iron-sulfur cluster coordination and disulfides across crystal images. A
supplied nonstandard RNA component can use the generic pucker modification.
Unknown components still need full typed restraint dictionaries.

The twelve frozen controls in `tests/native_validation_remaining` compare
electron and nuclear modes for one-coordinate zinc, tetrahedral Cys4 zinc,
SF4 with four Cys ligands, NAG C1-O4 linking, a supplied nonstandard RNA
component and a disulfide across a cell translation. All compare atom identity,
processed coordinates, target weights, report selection, Probe atom types,
clash pairs and reported scalar values. The symmetry controls compare the
fractional operators as well as the two directed ASU proxies.

Automatic links preserve monomer targets. Named deletion mods identify absent
leaving atoms, but must not retarget the remaining monomer geometry or delete
observed atoms. The live glycosidic controls caught that distinction. Explicit
supplied links take precedence over automatic discovery.

Nibbler retains metal and glycosidic restraints in the numeric graph. Its
validation report follows mmtbx's origin-0 selection for bonds, angles and
chiralities. `angle_origins` and `chirality_origins` expose that choice. The
explicit `RestraintGraph` evaluator can evaluate every supplied restraint.

## Limits

This is not a complete production replacement yet. General symmetry-linked
covalent and metal chemistry needs weighted image interpretation and reference
tests beyond the disulfide control. Automatic metal partner selection also
needs broader controls for alternate conformers, competing ligands, other
coordination counts and boundary distances. The preparation path includes
metal rules that CCTBX's validation path does not share; an unmatched weighted
bond remains incomplete.

Missing declared connections and undefined target uncertainties remain
incomplete. A structure with a typed dictionary does not automatically have
complete link chemistry. The full 1DWL control has an annotated MET-SD/HEC-FE
connection that the captured CCTBX graph also omits; Nibbler reports that gap.

The 1DWL comparison also exposed a reporting difference. mmtbx's plane wrapper
takes the maximum signed deviation. Nibbler uses the maximum absolute
standardized deviation, as the primitive plane tests specify. One captured
plane reports 8.78267 sigma through mmtbx and 15.63550 through Nibbler. Some
planes therefore produce additional native findings. The implementation must
not claim whole-entry legacy report parity for this case or silently remove
the absolute-deviation check.

The cancelled archive scan remains stopped. These local controls do not
measure archive rejection rates or establish a new archive ETA.
