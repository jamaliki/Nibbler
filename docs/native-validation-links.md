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

## Native validation path

The six agreed checks now use native interpretation, geometry and Probe
contacts. CCTBX remains the independent test oracle; the runtime neither
imports it nor starts a Probe2 process. The test suite blocks reference imports
while checking a complete mixed Cys/His zinc control.

Validation has its own nonbonded metal search. It applies covalent exclusions,
1-4 scaling, conformer compatibility, charged energy radii and the reference's
candidate order. Both MCL hooks read the same initial candidates. Preparation
keeps its separate metal rules. The search uses a spatial index and bounds
crystal shifts before enumerating images, including oblique cells.

The forty frozen controls in `tests/native_validation_mcl` include electron and
nuclear modes for mixed Cys/His counts, ND1/NE2 ordering, competing donors,
alternate conformers, both sides of the 3 A cutoff, FES, F3S, SF4, translated
zinc and annotated image bonds. They compare target weights, operator matrices,
processed atoms, Probe typing, clash pairs and reported values. The effective
zinc table includes an ordered mixed-nitrogen key that the live tests exposed;
its generator records the reference source hash.

Annotated connections now decode operator IDs and cell shifts and require a
matching weighted image. Unknown operators, wrong images and missing targets
remain incomplete. The default CCTBX automatic-link path skips general
symmetry-related covalent candidates. Nibbler does not invent those targets.

## Input limits and report differences

Missing declared connections and undefined target uncertainties remain
incomplete. A structure with a typed dictionary does not automatically have
complete link chemistry. The full 1DWL control has an annotated MET-SD/HEC-FE
connection that the captured CCTBX graph also omits; Nibbler reports that gap.

The 1DWL comparison also exposed a reporting difference. mmtbx's plane wrapper
takes the maximum signed deviation. Nibbler uses the maximum absolute
standardized deviation, as the primitive plane tests specify. One captured
plane reports 8.78267 sigma through mmtbx and 15.63550 through Nibbler. Some
planes therefore produce additional native findings. A recorded direct-primitive
test reproduces that discrepancy and verifies the absolute check. This case
does not have whole-entry legacy report parity; the implementation does not
remove the absolute-deviation check to match the wrapper.

The cancelled archive scan remains stopped. These local controls do not
measure archive rejection rates or establish a new archive ETA.

## Earlier build and timing evidence

The unpatched Nibbler build pins Reduce3
`33a8795bea9b43447c852b5984b9cc4ebc18c608`. It passed 145 Python tests with
zero skips, 106 Nibbler Rust tests and 57 owning backend tests. Nibbler library
Clippy passed; inherited backend warnings remain outside the new modules.

Twenty alternating paired trials included reads, parsing, interpretation,
geometry, Probe contacts and normal report output in warm workers. No agent
builds or tests ran during timing. The laptop remained shared; load ranged from
7.646 to 8.482 on 16 logical CPUs.

| Control | Atoms | Reference median | Native median | Wall speedup | CPU speedup |
| --- | ---: | ---: | ---: | ---: | ---: |
| Protein, nuclear | 642 | 1.46869 s | 6.26319 ms | 234.49x | 159.50x |
| Glycosidic, electron | 57 | 0.783635 s | 0.797958 ms | 982.05x | 272.02x |
| SF4, electron | 48 | 0.754951 s | 0.716834 ms | 1053.17x | 244.45x |

Native wall ranges were 5.935–6.868 ms, 0.709–1.016 ms and 0.584–0.981 ms.
The result file records individual trials and CPU time. RSS comes from the
shared process; it is not an isolated peak measurement for either backend.
See [the recorded results](native-validation-link-results.json).

## Completed native path, 12 October 2026

The unpatched release build pins Reduce3
`9bbb0423fc37c83216ed4e67a06325327bf2bdf7`. It passed 195 Python tests with
zero skips, 109 Nibbler Rust tests and 60 backend Rust tests. Nibbler library
Clippy passed. The full test-target Clippy run found inherited `expect` calls
in the old geometry tests; it does not pass with warnings denied.

Twenty alternating pairs per control used one thread per worker for Rayon,
OpenMP, OpenBLAS and MKL. No builds or tests ran alongside the measurements.
The shared laptop's one-minute load ranged from 6.937 to 8.359 on 16 logical
CPUs. The full 642-atom protein also passed live graph and contact comparisons
in electron and nuclear modes.

| Control | Atoms | Reference median | Native median | Wall speedup | CPU speedup |
| --- | ---: | ---: | ---: | ---: | ---: |
| Protein, nuclear | 642 | 1.47374 s | 6.33602 ms | 232.60x | 232.08x |
| Mixed Cys/His zinc, electron | 61 | 0.75624 s | 0.74075 ms | 1020.91x | 1004.74x |
| FES, electron | 24 | 0.71516 s | 0.36233 ms | 1973.77x | 1932.77x |

The measurements include input reads, parsing, fresh interpretation, restraints,
contacts and normal reports in warm workers. They do not predict archive
throughput. RSS shares a process between backends and does not measure their
separate peaks. Individual trials, worker thread settings, source pins and
binary hashes are in [the completed-path results](native-validation-complete-results.json).
