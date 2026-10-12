# Native model validation

`nibbler.geometry.validate` reads a parsed CIF document, builds chemical
restraints and evaluates geometry and Probe contacts in Rust. It releases the
GIL and caches immutable monomer libraries across calls and threads.

```python
from nibbler.geometry import validate

report = validate(document, chem_data="/absolute/chem_data", nuclear=False)
print(report["checks"])
print(report["failures"])
```

The checks cover bonds, angles, chirality, planes, supplied hydrogen parents
and all-atom clashes. Protein targets include peptide link modifications,
disulfide links and the conformation-dependent backbone library. Probe uses
chemical atom types, surface dots, covalent exclusions, alternate locations,
occupancy, hydrogen bonds, ionic contacts and temporary water hydrogens.

`complete` means that the checks have full native coverage. It does not mean
that the model has no defects. Every finding includes atom IDs, flat indices
and zero-based source atom rows; hydrogen findings also include parent IDs.
The atom IDs encode model, chain, residue number, insertion code, component,
alternate location and atom name as a JSON array string.

Set `include_graph=True` to inspect the independently built weighted graph,
chemical dictionary sources, initial Probe typing, directed contacts and
temporary contact atoms. The normal path does not build this full atom dump.
Validation neither adds real hydrogens nor repairs source coordinates.
Hydrogen preparation remains the separate `nibbler.reduce.run` operation.

## Reference coverage

Unknown components, unvalidated automatic links, ASU symmetry, nucleotide
pucker rules, cross-residue chirality and enol peptide targets remain explicit
fallback cases. Inputs needing CCTBX's symmetric atom-name swaps also require
the reference path. Custom restraint targets or atom-energy definitions in an
input document cannot disappear from coverage. Invalid formal charges and
multiple coordinate blocks return errors; callers must select each model or
assembly view before validation.

An incomplete check or a score near a cutoff must use reference validation.
Do not treat its diagnostic findings as evidence that all checks ran. This
API does not select training entries, change screening criteria or assign a
whole-entry rejection policy.

The reference harness compares native-built graphs against CCTBX 2026.9
targets and weights, then checks atom typing and Probe2 pairs independently.
Separate numerical tests cover the restraint evaluator and charged Probe
physics. These checks do not establish parity across the entire PDB archive.

CCTBX can round prepared hydrogen sites and swap symmetric atom labels during
model processing. Native validation keeps the supplied coordinates. Reference
fixtures that use processed sites include a manifest of every changed site
and the original source hash; the raw-input comparison remains separate.

## Performance checks

`benchmarks/native_validation_timing.py` measures the full input read, parse,
interpretation, restraint, contact and reporting path in warm workers. It
alternates backend order and records every trial, process CPU time, laptop
load and RSS. Run it after the same fixture passes independent parity and
after builds and tests finish. Fresh-process capture timing is separate and
includes imports and extra diagnostic output.

The Rust dependency uses an immutable Reduce3 revision. Its Apache and CCTBX
license notices ship with the wheel. No native speed measurement authorizes a
new scan or changes the training filter.

## Local results, 12 October 2026

Twenty warm paired trials per fixture, with alternating backend order, gave:

| Fixture | CCTBX/Probe2 median | Native median | Wall ratio | CPU ratio |
| --- | ---: | ---: | ---: | ---: |
| 1CRN, 642 atoms, nuclear H | 1.442 s | 6.075 ms | 237.4x | 173.3x |
| Protein, 2160 atoms, electron H | 2.986 s | 18.737 ms | 159.4x | 138.8x |

The reference retains the archive scanner's unchanged-atom guard. Its normal
path builds separate geometry and Probe models; native validation shares its
chemical graph. Both backends read files and emit normal reports. The larger
fixture uses the declared CCTBX-processed sites, with the original input and
every coordinate difference retained in the reference manifest.

Native wall times ranged from 5.758 to 6.295 ms for 1CRN and 18.014 to
19.737 ms for the larger fixture. Process CPU times ranged from 7.024 to
8.888 ms and 20.000 to 23.030 ms. Laptop one-minute load ranged from
29.54 to 40.43 on 16 logical CPUs. No agent build or test ran during these
measurements; other laptop work continued. An earlier 20-pair run measured
243.8x and 156.0x wall speedups, with CPU speedups of 180.6x and 136.6x.
The paired results retain every trial. These fixtures
do not establish archive throughput, large-assembly performance or the rate
of reference fallbacks. Shared-process RSS snapshots cannot compare isolated
backend peak memory.

The final wheel passed all 87 Python tests without skips. The independent
suite accounts for 31 tests, including full native
graph, Probe typing, clash pair and sigma5 scalar comparisons. The two protein
fixtures agree on all four and 134 geometry outliers, respectively. Native API
and numerical tests passed 36 tests, including concurrent calls and 750
seeded CCTBX primitive comparisons. All 107 Rust tests passed with the declared
chemistry tree, and library Clippy and formatting checks passed. All-target
Clippy still reports eight pre-existing `expect_used` lints in the earlier
`tests/geometry.rs` tests.

The 20 preparation tests passed in a separate chemistry tree containing the
missing MBO restraint file from GeoStd `f4d21cc5`. Its SHA256 is
`9fab45f7746c8278219b1ac37365e50f08061585ed5e1b02763aaa571a1d89f2`.
The benchmark tree stayed unchanged; hydrogen preparation code also stayed
unchanged. The implementation pins Reduce3
`4ebfea46adb4707f10441f2dd9f5c3d3cd260aff`.

The committed-build summary and hashes are in
[native-validation-results.json](native-validation-results.json). The full
timing records are in `benchmarks/native-validation-results/`; both runs
remain available. Reference tests check atom-level parity separately from
timing, so an incomplete fast result cannot establish a speedup.
