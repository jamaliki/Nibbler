<pre align="center">
          .-----.
         /   o   \
         \       /
          '-----'
    .-----.     .-----.
   /   o   \   /   o   \
   \       /   \       /
    '-----'     '-----'
</pre>

# Nibbler

Nibbler is a typed Python package backed by safe Rust for fast, strict CIF 1.1 analysis
and validated PDBx/mmCIF and ModelCIF output.

Nibbler 0.1.0 is an unreleased alpha. The supported Python contract targets CPython
3.10 through 3.14. Building from source requires Rust 1.88 or newer; the repository's
development environment pins Python 3.12.13 and Rust 1.97.1.

## What it does

- parses text CIF, gzip, and BinaryCIF 0.3 from paths, bytes, or binary streams;
- preserves complete logical documents or projects one category directly into typed,
  segmented Arrow buffers;
- distinguishes CIF unknown (`?`) from not applicable (`.`);
- applies equality, inequality, membership, and missing-state predicates during parse;
- uses pinned PDBx 5.416 or ModelCIF 1.4.9 dictionaries for projection typing and DDL2
  validation;
- exports through the Arrow C Stream protocol without requiring a dataframe library;
- scans many inputs through a bounded, deterministic native worker pool;
- parallelizes qualifying large loops while preserving strict grammar and error order;
- builds immutable PDBx coordinate and ModelCIF prediction models;
- resolves chemical components from embedded definitions, an explicit local CCD cache,
  or a small built-in registry, without network access;
- validates dictionary and semantic-profile invariants with stable diagnostics; and
- writes deterministic CIF or BinaryCIF transactionally.

Parsing is strict: Nibbler does not silently repair malformed syntax, guess chemistry,
or infer a document from a dataframe.

## Install from this checkout

Until 0.1.0 is published, install directly from a checkout with a compatible Python and
Rust toolchain:

```console
uv venv --python 3.12
uv pip install .
```

Activate the environment with `source .venv/bin/activate` on macOS or Linux, or
`.venv\Scripts\activate` on Windows. `python -m pip install .` is an equivalent fallback.
For development, use the repository-local micromamba environment described below.

Nibbler currently provides a Python API, not a command-line program. The commands under
`tools/` and `benchmarks/` are repository-development utilities rather than installed
user interfaces.

## Python API

The root facade has four operations: `chomp`, `feast`, `sniff`, and `dump`.

```python
import io

import nibbler
import pyarrow

# No category returns an immutable logical document.
document = nibbler.chomp("structure.cif.gz", schema="pdbx")
nibbler.sniff(document).raise_for_errors()

# One category returns a projected CifTable.
atoms = nibbler.chomp(
    "structure.bcif",
    category="atom_site",
    columns=["label_comp_id", "Cartn_x", "Cartn_y", "Cartn_z"],
    where={"pdbx_PDB_model_num": 1},
    schema="pdbx",
)
arrow_table = pyarrow.RecordBatchReader.from_stream(
    atoms.with_missing("columns")
).read_all()

# Output format is inferred from a path or selected explicitly for a stream.
nibbler.dump(document, "canonical.cif.gz", validate="dictionary")
nibbler.dump(document, "canonical.bcif", validate="dictionary")
binary_stream = io.BytesIO()
nibbler.dump(document, binary_stream, format="bcif")

# Preserving text output reuses source order and valid original value lexemes.
nibbler.dump(document, "source-order.cif", mode="preserve")

# Semantic models add chemistry and profile validation.
model = nibbler.mmcif.read(document, profile="pdbx")
nibbler.sniff(model, profile="pdbx").raise_for_errors()
nibbler.dump(model, "semantic.cif.gz", profile="pdbx")

registry = nibbler.components.read("components.cif")
model_with_local_ccd = nibbler.mmcif.read(
    "structure.cif", profile="pdbx", registry=registry
)

prediction = nibbler.mmcif.read("prediction.cif", profile="modelcif")
nibbler.dump(prediction, "prediction.cif.gz", profile="modelcif")

# An explicit viewer-oriented output may mirror one local QA metric into B factors.
nibbler.dump(
    prediction,
    "prediction-viewer.cif",
    profile="modelcif",
    mirror_local_qa_metric=1,
)
```

Scan results are bounded and yielded in input order:

```python
result = nibbler.feast(
    files,
    category="atom_site",
    workers=8,
    on_error="collect",
)
for table in result:
    consume(table)
for error in result.errors:
    log(error.code, error.source_name)
```

The searchable conventional APIs are `nibbler.cif.read`, `scan`, `validate`, and
`write`, plus `nibbler.mmcif.read`, `validate`, and `write`. `chomp` and `feast` are
exact aliases of `cif.read` and `cif.scan`.

The import package is `nibbler`. The distribution is named `nibbler-cif` because the
`nibbler` distribution name is used by another project.

The complete source, predicate, missing-state, validation, model-summary, and error
contracts are in the [Python API reference](docs/python-api.md).

## Development

The only host prerequisite is
[micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html). The
environment is repository-local.

```console
make bootstrap
make develop
make robustness
make check
make release-artifacts
```

Fetch and benchmark the hash-pinned PDB stress corpus:

```console
micromamba run -p .mamba/nibbler-dev python -m tools.fetch_pdb_corpus
make pdb-stress
```

Optionally build a correctness-trained, target-specific PGO wheel:

```console
micromamba run -p .mamba/nibbler-dev python -m tools.build_pgo
```

This native optimization workflow writes a wheel and compiler/source/corpus fingerprint
metadata to `dist-pgo/`. Profiles are compiler-, target-, ABI-, source-, and
workload-specific; no PGO flags are stored in `Cargo.toml`, and these wheels are not
part of the portable release matrix.

## Documentation

- [Current architecture and behavior](DESIGN.md)
- [CIF core algorithm diagrams](src/cif/README.md)
- [PDBx semantic algorithm diagrams](src/pdbx/README.md)
- [ModelCIF semantic algorithm diagrams](src/modelcif/README.md)
- [Python binding and scan diagrams](python/README.md)
- [Schema, PGO, and release workflow diagrams](tools/README.md)
- [Python API reference](docs/python-api.md)
- [Engineering and review standard](ENGINEERING.md)
- [Current performance architecture and qualification](docs/performance-report.md)
- [Benchmark commands and measurement contract](benchmarks/README.md)
- [Schema locks and artifacts](schemas/README.md)
- [Historical performance experiment ledger](docs/improvement-beam.md)
- [Release process and trusted publishing](docs/releasing.md)
- [Changelog](CHANGELOG.md)

## License

Nibbler is released under the [MIT License](LICENSE).
