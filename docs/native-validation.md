# Native model validation

`nibbler.geometry.validate` builds chemical restraints and evaluates bonds,
angles, chirality, planes, supplied hydrogen parents and all-atom Probe contacts
in Rust. It releases the GIL and caches immutable monomer libraries. It does not
import CCTBX or launch Probe2.

```python
from nibbler.geometry import validate

report = validate(document, chem_data="/absolute/chem_data", nuclear=False)
print(report["checks"], report["failures"])
```

`complete` means full coverage of these checks, not absence of defects. Findings
include atom IDs, flat indices and zero-based source rows. IDs encode model,
chain, residue number, insertion code, component, alternate location and atom
name. Set `include_graph=True` for weighted targets, chemical sources, initial
Probe typing and contact details.

Validation does not add real hydrogens or repair deposited coordinates.
Hydrogen preparation remains `nibbler.reduce.run`. Symmetric label normalization
swaps sites and B factors on a disposable copy and records those swaps.

## Chemistry and limits

The native path includes protein CDL and peptide modifications, RNA/DNA pucker
targets, private component/link/modification dictionaries, glycosidic links,
zinc coordination, iron-sulfur clusters and declared crystal-image connections.
Metal validation has a separate candidate search from hydrogen preparation.
Inline definitions do not change the shared library.

Missing typed restraints or annotated link targets remain incomplete.
Save-frame dictionaries, unsupported energy categories, malformed links,
invalid formal charges and multiple coordinate blocks fail explicitly.
Callers must select each training model or assembly view. Scores near a cutoff
also prevent a full-coverage claim. This API does not assign training exclusions.

The default report follows mmtbx origin-0 selection for bonds, angles and
chirality. The numeric graph retains metal and other link restraints; the
explicit graph evaluator can check all supplied restraints. A scanner must
choose that reporting scope.

Native planes use maximum absolute standardized deviation. The mmtbx wrapper
uses maximum signed deviation. The recorded 1DWL ASP B66 control gives 15.6355
versus 8.78267 sigma; the direct CCTBX primitive agrees with native scoring.
Native can therefore report additional plane findings. The annotated
MET-SD/HEC-FE connection in 1DWL lacks a weighted target in the captured reference
graph too. Neither case establishes whole-entry report equivalence.

## Independent test data

The extension suite retains all 84 CCTBX 2026.9/Probe2 controls in a separate
[versioned data asset](https://github.com/jamaliki/Nibbler/releases/tag/native-validation-data-20261012-v1).
Inputs with identical bytes share one stored file. The bundle preserves the
original input, dictionary and reference hashes; no expected value comes from
the implementation under test. Raw timing records and capture generators
remain in its `evidence/` directory, outside the Git source tree.

Fetch explicitly before running the suite. Tests never download data.

```sh
python -m tools.fetch_validation_data
NIBBLER_VALIDATION_CHEM_DATA=/absolute/chem_data \
    python -m pytest -q tests/test_native_validation_extensions.py
```

The fetch command pins the archive size and SHA256, verifies fixture hashes
and rejects unsafe archive paths. Its default destination is
`.cache/native-validation-v1`. Use `--archive` for a local bundle,
`--destination` for another cache and `NIBBLER_VALIDATION_DATA` to point tests
there. With configured chemistry, missing test data fails rather than silently
omitting these controls. Without chemistry, the suite reports explicit skips.

Controls compare atom identities and sites, graph targets and weights,
chirality flags, image operators, Probe typing, clash pairs and scalar findings.
They include modified residues, nucleotides, ligands, metals, alternate
conformers, supplied targets and cutoff cases. They do not establish archive-wide
coverage. MonomerLibrary and GeoStd pins are `1c97ea23` and `f4d21cc5`.

## Performance evidence

Twenty alternating warm pairs used one thread per worker for Rayon, OpenMP,
OpenBLAS and MKL. Measurements included reads, parsing, fresh interpretation,
geometry, contacts and normal reports. No builds or tests ran alongside them;
the shared laptop's one-minute load ranged from 6.937 to 8.359 on 16 logical CPUs.

| Control | Atoms | Reference median | Native median | Wall ratio | CPU ratio |
| --- | ---: | ---: | ---: | ---: | ---: |
| Protein, nuclear | 642 | 1.47374 s | 6.33602 ms | 232.60x | 232.08x |
| Mixed Cys/His zinc, electron | 61 | 0.75624 s | 0.74075 ms | 1020.91x | 1004.74x |
| FES, electron | 24 | 0.71516 s | 0.36233 ms | 1973.77x | 1932.77x |

The full protein passed live comparisons in both hydrogen modes. These results
do not predict archive throughput or include hydrogen preparation.
Shared-process RSS does not measure each backend's isolated peak. All trials,
source pins and binary hashes remain in the external data asset. The runtime
pins Reduce3 `9bbb0423fc37c83216ed4e67a06325327bf2bdf7`; its Apache and CCTBX
notices ship with the wheel. The cancelled scan remains stopped.
