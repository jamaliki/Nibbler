# Observations

This module packs dense training observations from a PDBx document. kaveh-data
calls it through `nibbler.observations.ObservationPacker`. The bindings live in
`src/python/observations/`. The caller supplies the schema, the ligand token
table, the excluded components and the per-component chemistry. The packer reads
selected atom rows and writes the arrays. It infers no missing coordinate and
derives no ligand token.

The packed bytes are a contract with kaveh-data (`CONTRACT_VERSION` in
`src/lib.rs`). A change that alters any packed array, key or error text needs a
new contract version and a matching kaveh-data change.

## Features

The module builds under the `observations` Cargo feature. The `python` feature
turns it on. Files that read Reduce3 chemistry (`bridges`, `molecular/`,
`parents`) also need `reduce3`. The Rust tests run without libpython:

```sh
make test   # runs `cargo test --features observations,reduce3` among the others
```

Use `--features observations` alone to check that the module builds without
Reduce3. Only the Python layer calls the packer, so builds without `python`
allow dead code in this module.

## Files

| File | Responsibility |
| --- | --- |
| `mod.rs` | Entry points the bindings call: packing, components, chemical sources |
| `schema.rs` | Caller-owned atom slots and allowed residue states |
| `token_table.rs` | Ligand token identities named by the caller |
| `prepared.rs` | One immutable source snapshot for component queries and packing |
| `source.rs` | Selected atom rows and the covalent links training targets use |
| `record.rs` | Owned dense arrays, complete-residue selection and identity alignment |
| `residue.rs` | Observed atoms packed into slots |
| `tokens.rs` | Protein residue tokens and heavy ligand atom tokens |
| `topology.rs` | Deposited links and observed peptide continuity |
| `continuity.rs` | Sequence gaps kept and gaps added at observed polymer breaks |
| `product.rs` | Occurrence-specific product targets applied before token selection |
| `uncertainty.rs` | Source-row diagnostics applied to only the affected targets |
| `crops.rs` | Finite spatial-crop plans over caller-owned token arrays |
| `afdb.rs` | Protein documents rebuilt from dense heavy-atom arrays |
| `afdb_eligibility.rs` | Complete-chain crop eligibility read from A14P bytes |
| `parents.rs` | Hydrogen dictionaries for each retained ligand residue (`reduce3`) |
| `bridges.rs` | Dictionary-backed metal hydrides without a single parent slot (`reduce3`) |
| `molecular/` | Atom84 molecular targets from one source table and one linked product (`reduce3`) |

`*_tests.rs` files and the `tests.rs` files in subdirectories hold the inline
tests for the file they name.

## Determinism

Maps keyed by `FastMap` and `FastSet` (`src/cif/fast_hash.rs`) never decide
output order. Code that emits rows from a map sorts them or walks a `BTreeMap`
or `BTreeSet`. Keep that rule when adding a map here.
