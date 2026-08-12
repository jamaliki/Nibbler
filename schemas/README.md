# Schema locks and compiled dictionaries

[`locks.toml`](locks.toml) pins the exact PDBx/mmCIF 5.416 and ModelCIF 1.4.9
dictionaries used by Nibbler. The large upstream dictionary text is not vendored. Fetch
and verify it into the ignored local cache with:

```console
make schemas
```

Runtime parsing, projection, validation, and writing never fetch schemas from the
network.

Compile the verified DDL2 inputs through Nibbler's own CIF parser with:

```console
make compile-schemas
```

The deterministic [`compiled/`](compiled/) `.nbs` files are committed release inputs.
They contain the definitions needed for projection typing and validation: item and
category identity, DDL2 type patterns, mandatory flags, enumerations, numeric ranges,
keys, and parent-child links. The native extension checks each artifact's format
version and source digest before lazy initialization.

Updating a lock requires an explicit version and digest change, regenerated artifacts,
and conformance review:

```console
make schemas
make compile-schemas
make validate-fixtures
```

Fixture and corpus tests pin the lock metadata so an unreviewed dictionary change fails
the repository checks.
