# Reduce3 chem_data test data

Reduce3 and native validation read a cctbx `chem_data` directory: `geostd/`, `mon_lib/`
and `chemical_components/`. The full tree does not belong in this repository, and the
tests need only a small part of it. The [manifest](../tests/chem_data_subset.json) pins
200 files by byte count and SHA-256:

- GeoStd at `phenix-project/geostd` commit `f4d21cc5` and the Monomer Library at
  `MonomerLibrary/monomers` commit `1c97ea23`, the chem_data 2026.9 pins that
  [native validation](native-validation.md) names. Both come from immutable
  `raw.githubusercontent.com` commit URLs.
- CCD definitions from `https://files.rcsb.org/ligands/download/<COMPONENT>.cif`.

The subset holds the list and energy-library files of both libraries and the GeoStd
`rna_dna/` link files. It holds every component of the Rust and Python fixtures, the
native validation data, and the bridge tests (`EJ2`, `NWN`), plus the GeoStd and Monomer
Library files of the synonyms that the list files map those components to, such as `AD`
for `DA`. Each component takes its GeoStd `data_X.cif` with the `_neutron` and `_pH_low`
variants, its Monomer Library `X.cif`, and its CCD file wherever upstream has one.

```console
micromamba run -p .mamba/nibbler-dev python -m tools.fetch_chem_data
REDUCE3_CHEM_DATA=.cache/chem-data-v1 NIBBLER_VALIDATION_CHEM_DATA=.cache/chem-data-v1 \
  NIBBLER_REQUIRE_CHEM_DATA=1 micromamba run -p .mamba/nibbler-dev python -m pytest
```

The ignored cache is `.cache/chem-data-v1`. A warm, verified cache needs no network
access, and `--offline` verifies it without any. `make test` fetches the subset and the
[validation data](native-validation.md), and runs every Rust and Python test against
them. CI restores the cache under a key derived from the manifest.

Tests that need chem_data skip when they cannot find it. `NIBBLER_REQUIRE_CHEM_DATA=1`,
which `make test` and CI set, turns those skips into failures, so a missing or
misplaced tree cannot hide the Reduce3 and validation tests.

RCSB serves only the current revision of a component. When RCSB revises one, its
download stops matching and the fetch fails rather than silently changing a test. Run
`python -m tools.fetch_chem_data --repin` to record the current sizes and digests, then
review the changed components and the test results before committing the manifest.
