# Changelog

All notable user-facing changes are recorded here. Nibbler follows Semantic Versioning.

## [0.1.0] - Unreleased

### Added

- Strict CIF 1.1 parsing with complete logical-document construction and preserving or
  canonical text output.
- Direct typed Arrow projection with predicates, missing-state preservation, and
  deterministic parallel parsing of large loops.
- Content-detected gzip and BinaryCIF 0.3 input, plus deterministic BinaryCIF output.
- Pinned PDBx 5.416 and ModelCIF 1.4.9 dictionary validation.
- Immutable PDBx and ModelCIF semantic models, explicit chemical-component resolution,
  canonical writers, and optional ModelCIF local-QA mirroring.
- Bounded concurrent scanning, structured diagnostics, and transactional output.
- Bounded property tests, a hash-pinned real-PDB stress corpus, and profile-guided native
  release tooling.
- Portable CPython 3.10-3.14 wheels, isolated artifact qualification, trusted
  publishing, and scheduled PDB correctness, throughput, and peak-RSS regression checks.
- Optional Reduce3 integration (the `reduce3` Cargo feature, enabled in the Python
  package): `nibbler.reduce.run` and `reduce::run` add and optimize hydrogens on a parsed
  document in memory and return a new document that keeps the source block's categories.
- The Python extension allocates with mimalloc, and the CIF text writers borrow values
  they write unchanged instead of copying each one (same output). On a 225-residue
  AlphaFold model, `nibbler.reduce.run` takes 13 ms and `to_canonical` 1.9 ms (from 26 ms
  and 5 ms). Reduce3's vectorized dot scoring takes a single-threaded run from 19.3 ms to
  17.4 ms (same output), and Reduce3 is built without its command-line-only features.
- `nibbler.reduce.run` and `reduce::Params::default()` consider Asn/Gln/His flips by
  default (`add_flip_movers=True`; Reduce2's default is False).
- Reduce3's fixed mode keeps a hydroxyl hydrogen on a planar atom (tyrosine and other
  phenols, enols, carboxylic acids) in that atom's plane and an acid's hydrogen syn to
  its carbonyl oxygen; `nibbler.reduce.run` takes `planar_hydroxyl_preference` and
  `acid_syn_preference` (1.0 each, 0 turns them off; compat mode never applies them).
- Reduce3's fixed mode no longer leaves a histidine with neither ring hydrogen (an
  imidazolate), and removes the thiol hydrogen of a cysteine bound to any metal, not
  just zinc. Against 334 neutron structures, histidine protonation (where the ring is
  oriented as deposited) matches for 62.6%, up from 56.3%, and all 36 metal-bound
  cysteines match, up from none.
- `nibbler.mmcif.assembly` and `pdbx::assembly_document` write one biological assembly
  out as explicit coordinates: every copy of the asymmetric units it uses, chains of copy
  *n* renamed `X-n`, coordinates transformed and anisotropic tensors rotated, atoms on
  symmetry axes written once, without the crystal cell or symmetry. `nibbler.reduce.run` on the result places hydrogens with the
  contacts between copies taken into account.
- Immutable revision-addressed PDB text fixtures, public contribution and security
  policies, checkout installation with `uv`, and original project artwork.

[0.1.0]: https://github.com/jamaliki/Nibbler/releases/tag/v0.1.0
