# Independent native validation references

The captures use CCTBX2026.9 and the declared chemistry tree, with resolved
standard/GeoStd energy and monomer list hashes in the manifest. The legacy
Reduce monomer list has a separate hash. These local assets differ from the
cancelled cluster scan; these captures do not validate that scan's setup.

The compressed JSON references preserve the final processed graph, atom
identities, coordinates, weights, targets, origin IDs, chirality policy flags,
and scalar results. Probe2 captures include initial atom metadata and covalent
neighbors before donor/ionic/phantom processing, its actual internal PHIL,
raw contact output and parsed clash pairs. Capture hooks return the original
objects unchanged and restore their functions before leaving the process.

The comparator maps atom identities rather than assuming parser order. It
checks required graph targets/weights, pseudochirality flags, initial Probe
radii/type flags, rounded pairs and clashscores separately. Unsupported or
missing native checks do not count as parity. ASU proxies and strict chirality
order remain explicit limits. The scalar comparison checks strict sigma>5
findings, excludes pseudochirality, and uses maximum plane score rather than
RMS score. It rejects duplicate findings. Dihedral proxies do not enter
the four-category graph parity claim.

Eight positive/negative chemical controls, prepared five-residue crambin
with and without CDL, nuclear-mode 1CRN, a phantom-water case, and an explicit
CCTBX-sites 2160-atom protein pass independent graph/type/pair and sigma5
scalar checks with Reduce3 `4ebfea46`. The prepared 1CRN input uses nuclear H
distances; the larger protein uses electron distances. The nuclear capture
found 46 N-H target and 246 Probe radius mismatches in the earlier v6 wheel.
The corrected wheel passes that unchanged oracle, including zero clashes.
The microheterogeneity case records Probe adding ten H atoms; it does not
establish native completeness for an input without those atoms. The unknown
component case records failure instead of an empty successful graph.

CCTBX processing rounds the original 2160-atom prepared input's H coordinates
from five decimals to three. Native validation preserves the input values.
The original comparison retains that maximum 0.0005 A difference. A separate
fixture applies the recorded CCTBX sites to the original CIF while retaining
all categories. Its manifest lists every site difference and both hashes.
The common-sites comparison isolates evaluator parity; it is not a native
coordinate repair or a change to deposited data.

Run portable frozen-data tests without a CCTBX installation. Live tests
require an explicit `NIBBLER_REFERENCE_CHEM_DATA` tree and a native-validation
build. The CCTBX unknown-component test also requires the reference toolchain.
An unset tree or missing toolchain causes only those live tests to skip;
the recorded local validation ran with the tree set and no skips.

```sh
NIBBLER_REFERENCE_CHEM_DATA=/absolute/immutable/chem_data \
  python -m pytest -q tests/native_validation_reference_test.py
```

The fresh-process benchmark alternates AB/BA order and includes imports,
setup, validation and JSON output. Reference capture adds instrumentation
overhead, so these timings must remain separate from normal API timings.
No timing alone establishes scientific parity or whole-PDB coverage.
