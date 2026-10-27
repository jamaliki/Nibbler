# Changelog

All notable user-facing changes are recorded here. Nibbler follows Semantic Versioning.

## [0.1.0] - Unreleased

### Added

- `ObservationPacker` takes the caller's ligand token table,
  `{"elements": {symbol: token}, "other": token, "hydrogen": token}`, in place of
  an element list. Nibbler derives no ligand token: before, element `i` of the
  list took `21 + i`, other elements 64, and a bridging metal hydride 65, which
  only held while the list was 43 elements long. The table also replaces the
  molecular schema's `ligand_elements`, which is now ignored. The constructor
  refuses tokens that repeat or fall in the residue range 0-20, hydrogen as an
  element, and lower-case symbols. The 66-token table reproduces the list's
  legacy and molecular packs byte for byte on the test fixtures
  (`tests/test_observation_token_table.py`, recorded at 9b8041e).

- `crop_plans` takes `ligand_limits`, the most ligand tokens and ligand chains a
  crop may keep, with the former (64, 9) as the default. It also returns the
  ligand tokens and chains each accepted crop keeps. `ligand_limits=None` plans
  every crop that keeps a protein chain and also reports the counts of every
  eligible ligand centre's crop. Rejecting the uncapped plans over chosen limits,
  and counting the ligand centres over them as retries, gives that limit's plans
  and weights bit for bit, so callers can choose the limits after planning.

- Seeded `random_clash_aware` conformer selection varies residue alternatives
  even when occupancies differ. It preserves source rows and occupancy values,
  keeps explicitly qualified inter-residue link alternatives together, and
  resolves authored insertion codes. The default remains occupancy-first.

- Native ligand-only SASA through `mmcif.ligand_sasa`, with source-row areas,
  decoded-model reuse, explicit ligand selection, adjustable radii and sampling,
  and a coordinate-array API at `geometry.sasa`. Parsing performs no SASA work
  unless the caller requests this separate operation.

- Reduce3's fixed mode names the hydrogens of Asn, Gln and Arg NH2 groups and of
  base amino groups as IUPAC does (HD21, HE21, HH11 and HH21 cis to CB, CG and NE;
  H61, H41 and H21 cis to N1, N3 and N1), whether or not flips are considered.
  Before, HD21/HD22 and HE21/HE22 swapped between runs with and without flips,
  and every Arg was named the other way. After an Asn, Gln or His flip docks the
  side chain, HA is placed again from the final CB, by up to 0.07 A. Scores, flip
  choices and Reduce2 compatibility mode remain unchanged.

- Reduce3's fixed mode limits conjugated O-H displacement to 0.12 A from both
  local planes and checks the final coordinates independently of contact scores.
  The 6TBR regression retains its Gln contact without a large TYR224 H tilt.
  Ordinary alcohols, thiols and Reduce2 compatibility mode remain unchanged.

- Reduce3 indexes deposited bond endpoints within each unchanged atom hierarchy,
  avoiding repeated full-structure scans while keeping exact match rules and
  candidate order.

- Explicit native linked-product chemistry, clash-aware conformer selection, and
  whole-ligand assembly overlap checks, with source mappings and real-entry regressions.

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
- Reduce3's fixed mode favors a hydroxyl hydrogen in its planar partner's plane
  (tyrosine and other phenols, enols, carboxylic acids), and an acid's hydrogen syn to
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
