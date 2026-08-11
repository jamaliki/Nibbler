# Schema locks

`locks.toml` pins the exact PDBx/mmCIF and ModelCIF dictionaries targeted by Nibbler.
The dictionaries are large upstream inputs and are not vendored. Fetch and verify them
into the ignored local cache with:

```console
make schemas
```

Updating a lock requires an explicit version change, digest change, regenerated schema
artifacts once the compiler exists, and conformance review. Runtime parsing and writing
must never fetch schemas from the network.
