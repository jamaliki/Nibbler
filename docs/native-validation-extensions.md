# Native chemistry extensions

The native API now accepts supplied component, link, modification and
`lib_atom` definitions. Each call uses a private overlay. Overlays share
immutable library tables and file-backed caches; a changed table gets its
own copy. Unknown components still need full, typed restraints. CCD bond
orders and ideal coordinates alone do not provide those restraints.

```python
report = validate(
    document,
    chem_data="/absolute/chem_data",
    restraint_dictionaries=["component-and-link.cif"],
    restraint_links=[("LINK_ID", first_source_atom_row, second_source_atom_row)],
)
```

The zero-based atom rows select the two residue sides in dictionary order.
The API validates the link ID, atom rows, residue match and repeated
applications before building the weighted graph. It applies declared residue
modifications before emitting monomer targets and updates the affected atom
types. A `struct_conn` bond alone does not request the full link dictionary;
CCTBX also distinguishes bond annotations from explicit CIF link applications.
Annotated bonds without weighted native definitions remain incomplete.

Standard A, C, G, U, DA, DC, DG and DT targets now use the effective library
names. Validation preserves the case of monomer aliases because CCTBX
distinguishes classic `Ad` dictionaries from `AD` aliases. Hydrogen preparation
keeps its prior alias rules. The RNA path applies `rna2p`/`rna3p` modifications
with the next phosphate when available.

Label normalization covers ARG, ASP, GLU, PHE, TYR, VAL and LEU, including
paired H/D names. It changes a disposable structure copy and reports each
coordinate/B-factor exchange. A missing member of a swap pair leaves coverage
incomplete. The API does not rewrite deposited structures.

## Independent comparisons

Thirty-two electron/nuclear comparisons pass graph targets and weights, Probe
atom typing, clash pairs/scores and strict sigma5 findings. They cover eight
standard nucleotide bases, an RNA pair with both pucker choices, a deliberately
swapped whole protein, a supplied custom component, typed isolated/linked
enols, a cross-component chirality link and a link with a residue modification.
Two of the controls change the supplied enol bond and angle targets; CCTBX and
the native code preserve those targets when HNO already has a declared type.
The fixtures are controls, rather than a representative archive sample.

The enol tests also found a CCTBX/Probe limitation. CCTBX added HNO geometry
targets, but Probe lost that atom's type and used a generic hydrogen radius.
The reference harness now marks that fallback incomplete. The successful
comparisons supply the same explicit HNO type and restraints to both engines.
Raw fallback captures remain in the local experiment evidence.

The frozen sources, dictionary hashes and packed references are in
`tests/native_validation_extensions/manifest.json`. Live regeneration uses
`benchmarks/native_validation_extension_fixtures.py`,
`native_validation_link_fixture.py` and `native_validation_extension_matrix.py`.
Library sources are MonomerLibrary/monomers
`1c97ea23a7df06373ed01cd7c3156a1f8e9b27e3` and phenix-project/geostd
`f4d21cc5f4dbf84f55cc972fb6b287354300ab2d`; the reference uses CCTBX 2026.9.

Run the regression cases with the declared expanded chemistry tree:

```sh
NIBBLER_VALIDATION_CHEM_DATA=/absolute/chem_data \
    python -m pytest -q tests/test_native_validation_extensions.py
```

## Coverage still requiring reference checks

The later [link extension](native-validation-links.md) covers weighted metal
and glycosidic controls, a supplied nonstandard RNA component and disulfide
image operators. General symmetry-linked chemistry, incomplete dictionaries,
ambiguous conformer links and unsupported dictionary categories remain
incomplete or return
an explicit input error. This extension does not establish full native archive
coverage, change geometry cutoffs, change training inputs or restart the
cancelled scan.

## Build and timing evidence

The extension build and local measurements appear in
`native-validation-extension-results.json`. Full-path timings include reads,
parsing, chemistry, geometry, Probe contacts and normal reporting. They
alternate backend order and retain wall time, process CPU time, laptop load
and RSS for every trial. Detailed parity captures stay outside timed regions.
Shared-process RSS does not measure each backend's isolated peak memory.

The pinned build passed 133 Python tests without skips, all 108 Nibbler Rust
tests and 55 owning Reduce3 tests. Nibbler library Clippy and the changed-file
format/lint checks passed. Reduce3 still has inherited library Clippy warnings
outside the new validation modules.

Twenty paired warm trials per input measured:

| Input | CCTBX/Probe2 median | Native median | Wall ratio | CPU ratio |
| --- | ---: | ---: | ---: | ---: |
| 1CRN, nuclear H, 642 atoms | 1.486 s | 6.349 ms | 234.1x | 161.9x |
| Protein, electron H, 2160 atoms | 3.036 s | 18.949 ms | 160.2x | 136.1x |
| RNA pair, electron H, 66 atoms | 0.771 s | 1.097 ms | 702.9x | 217.8x |

No agent build or test ran during timing; other laptop work continued.
One-minute load ranged from 7.65 to 8.48 on 16 logical CPUs. Native wall-time
ranges were 5.998-7.135 ms, 18.322-19.955 ms and 0.967-1.292 ms, respectively.
The protein input uses the earlier recorded CCTBX-processed sites. These
measurements do not estimate the runtime of an archive scan with unknown
fallback frequency and hydrogen-preparation costs.
