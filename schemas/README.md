# Schema locks

`locks.toml` pins the exact PDBx/mmCIF and ModelCIF dictionaries targeted by Nibbler.
The dictionaries are large upstream inputs and are not vendored. Fetch and verify them
into the ignored local cache with:

```console
make schemas
```

Updating a lock requires an explicit version change, digest change, regenerated schema
artifacts, and conformance review. Runtime parsing and writing never fetch schemas from
the network.

Compile the verified DDL2 inputs through Nibbler's own CIF parser with:

```console
make compile-schemas
```

The deterministic `compiled/*.nbs` artifacts are committed release inputs. Their
format version, source digest, and included extension dictionaries are recorded in
`locks.toml`.
