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
  document in memory and return a new document.
- Immutable revision-addressed PDB text fixtures, public contribution and security
  policies, checkout installation with `uv`, and original project artwork.

[0.1.0]: https://github.com/jamaliki/Nibbler/releases/tag/v0.1.0
