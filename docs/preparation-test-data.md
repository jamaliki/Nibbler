# Linked-product and structure-preparation test data

The public coordinate and CCD files do not belong in this repository. The pytest
session fixture downloads each requested entry from RCSB and builds the 112-component
CCD subset only when a product test needs it. The small
[manifest](../tests/preparation_corpus.json) pins byte counts and SHA-256 hashes.
Coordinate hashes cover uncompressed CIF content, so gzip timestamps do not matter.
CCD hashes cover the downloaded component bytes and the assembled subset.

Downloads use RCSB's [file services](https://www.rcsb.org/docs/programmatic-access/file-download-services),
`https://files.rcsb.org/download/<ID>.cif.gz` and
`https://files.rcsb.org/ligands/download/<COMPONENT>.cif`. A changed download fails
verification and needs a reviewed manifest update. It never silently changes a test.

The ignored cache is `.cache/structure-preparation-v1`. Each pytest session verifies
used files once, and later access avoids rereads. A valid assembled CCD avoids all
112 component downloads and reads on later runs. CI reuses the cache under a key
derived from the manifest. Small synthetic Rust tests need no network access.

```console
micromamba run -p .mamba/nibbler-dev python -m tools.fetch_preparation_data
NIBBLER_TEST_DATA_OFFLINE=1 micromamba run -p .mamba/nibbler-dev python -m pytest \
  tests/test_linked_products.py tests/test_conformer_selection.py
```

`NIBBLER_TEST_DATA` selects another cache directory. The fetch command accepts
`--destination` and `--offline`. Missing or changed offline data fails rather than
skipping regressions. No installed Nibbler API fetches these files.

The case selection follows AtomWorks commit
`1b731b247f750e2f7117806a96b5144f94c898b2` and its migration guide. The 31 product
cases include Schiff bases, glycosylation, linked amides, multiple-bond additions,
hydrogen transfer, caps, resonance, phosphorus/boron/selenium charges, and reactions
that must fail. `1VQ7` and `7P1D` contain valid CIF syntax but unsupported products.
The pinned `1TWR` download includes RCSB's metalloprotein remediation, with corrected
`VER` hydrogen templates. Its regression checks that each carbon has one set of
template hydrogens.

`1CBN`, `6LZB`, and `1RCQ` test conformer selection. `6QHP` must fail without an
explicit parent-component substitution policy. `7MUB` and `1XAN` test whole-ligand
assembly overlap. Their deposited metadata has pre-existing profile errors; the
assembly test checks that expansion introduces no new error codes.
