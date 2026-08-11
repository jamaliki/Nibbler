# Nibbler

Nibbler is a native Python toolkit for fast, robust CIF analysis and for writing
correct PDBx/mmCIF and ModelCIF files.

> [!IMPORTANT]
> Nibbler is in contract-first development. The Phase 0 package establishes its API,
> diagnostics, fixtures, and engineering gates. Parsing and writing intentionally fail
> with `FeatureUnavailableError` until the shared Rust implementation is complete.

## Intended API

```python
import nibbler

document = nibbler.chomp("structure.cif")
nibbler.sniff(document, profile="pdbx").raise_for_errors()
nibbler.spit(document, "canonical.cif", profile="pdbx")
```

The themed root facade is backed by conventional APIs:

```python
from nibbler import cif

document = cif.read("structure.cif")
cif.write(document, "canonical.cif")
```

The Python import package is `nibbler`. The distribution is named `nibbler-cif`
because the `nibbler` distribution name is already used by an unrelated project.

## Development

The only host prerequisite is
[micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html).
The repository-local environment does not modify a shared Python installation.

```console
make bootstrap
make develop
make check
```

List the frozen benchmark workloads and verify the fixture corpus:

```console
make corpus
make benchmarks
```

## Project contracts

- [Detailed architecture and behavior](DESIGN.md)
- [Implementation and review standard](ENGINEERING.md)

## License

Nibbler is released under the [MIT License](LICENSE).
