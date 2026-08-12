# Nibbler

Nibbler is a native Python toolkit for fast, robust CIF analysis and for writing
correct PDBx/mmCIF and ModelCIF files.

> [!IMPORTANT]
> Nibbler is in phased development. Phases 1 through 6 implement strict CIF parsing,
> schema-typed columnar projection, Arrow interchange, gzip input, bounded Python batch
> reads, pinned DDL2 validation, transactional writing, chemistry-aware PDBx models,
> prediction-aware ModelCIF models, BinaryCIF, preserving output, and profiled fast
> paths for large coordinate projections.

## Implemented core

The Rust API in `_core::cif` currently provides:

- immutable, reference-counted UTF-8 source buffers;
- a single zero-copy lexer and strict state-machine parser;
- ordered data/global blocks, scalar items, distinct loops, and save frames;
- separate unknown (`?`) and not-applicable (`.`) values;
- stable located parse errors and configurable resource limits;
- projection into Arrow-native UTF-8 buffers without retaining unselected values;
- lossless internal missing states with collapse, companion-column, and extension
  policies at the Arrow boundary;
- content-detected gzip with decompressed-size and expansion-ratio limits;
- content-detected BinaryCIF 0.3 decoding and deterministic encoding behind the same
  document, projection, semantic-model, and Python APIs;
- GIL-free single-source reads and deterministic bounded native file parallelism;
- deterministic compiled artifacts for pinned PDBx 5.416 and ModelCIF 1.4.9 schemas;
- lazy, lock-verified dictionary validation for known items, DDL2 types, mandatory
  definitions, enumerations, ranges, keys, and parent-child links;
- schema-guided integer and floating-point Arrow builders without heuristic inference;
- transactional canonical or preserving path output and content-detected gzip round
  trips, with valid source value lexemes retained in preserve mode;
- immutable PDBx entities, asym units, polymer sequences, non-polymers, water,
  branched glycans, components, atom sites, and explicit connections;
- deterministic embedded, caller-selected local CCD, minimal-registry, or unresolved
  component resolution without network access;
- strict PDBx semantic validation and canonical profile ordering, with unresolved
  chemistry blocked before output;
- immutable ModelCIF targets, prediction models and groups, versioned software,
  protocol data flow, templates, alignments, QA metrics, and associated archives;
- strict ModelCIF reference and completeness validation, deterministic profile ordering,
  and opt-in local-confidence mirroring without conflating confidence with B-factor
  semantics;
- external qualification of source and Nibbler-written chemistry goldens with Gemmi
  and the pinned PDBe validator, plus Python-ModelCIF qualification for predictions;
- deterministic CIF 1.1 text formatting and canonical logical round trips; and
- property tests plus lexer, parser-round-trip, and formatter fuzz targets.

## Python API

```python
import io

import nibbler

document = nibbler.chomp("structure.cif")

atoms = nibbler.chomp(
    "structure.cif.gz",
    category="atom_site",
    columns=["label_comp_id", "Cartn_x", "Cartn_y", "Cartn_z"],
    where={"pdbx_PDB_model_num": 1},
    schema="pdbx",
)
arrow_table = atoms.to_pyarrow(missing="columns")

checked = nibbler.chomp("structure.cif", schema="pdbx")
nibbler.sniff(checked).raise_for_errors()
nibbler.dump(checked, "canonical.cif.gz", validate="dictionary")
nibbler.dump(checked, "canonical.bcif", validate="dictionary")

# Input format and gzip are detected from content. Binary output is inferred from a
# .bcif or .bcif.gz path, or selected explicitly for a binary stream.
binary_atoms = nibbler.chomp("structure.bcif", category="atom_site")
binary_stream = io.BytesIO()
nibbler.dump(checked, binary_stream, format="bcif")
nibbler.dump(checked, "source-order.cif", mode="preserve")

model = nibbler.mmcif.read("structure.cif")
nibbler.sniff(model, profile="pdbx").raise_for_errors()
nibbler.dump(model, "semantic.cif.gz", profile="pdbx")

registry = nibbler.components.Registry.from_ccd_cache("components.cif")
model_with_local_ccd = nibbler.mmcif.read("structure.cif", registry=registry)

prediction = nibbler.mmcif.read("prediction.cif", profile="modelcif")
nibbler.sniff(prediction, profile="modelcif").raise_for_errors()
nibbler.dump(prediction, "prediction.cif.gz", profile="modelcif")

# Optional compatibility view for legacy structure viewers. The named QA metric is
# retained in ModelCIF; only output B-factor values are mirrored.
nibbler.dump(
    prediction,
    "prediction-viewer.cif",
    profile="modelcif",
    mirror_local_qa_metric=1,
)

result = nibbler.feast(files, category="atom_site", workers=8, on_error="collect")
for batch in result.batches:
    consume(batch)
diagnostics = result.errors
```

The themed root facade is backed by conventional APIs:

```python
from nibbler import cif

document = cif.read("structure.cif")
report = cif.validate(document, schema="pdbx")
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

Run the correctness-qualified benchmark smoke test and verify the fixture corpus:

```console
make corpus
make benchmarks
python -m tools.fetch_pdb_corpus
python -m benchmarks.pdb_stress --warmups 1 --samples 3
```

The benchmark target installs an optimized release extension before measuring it.

Release maintainers can build a native, corpus-trained PGO wheel after fetching the
pinned PDB corpus:

```console
micromamba run -p .mamba/nibbler-dev python tools/build_pgo.py
```

The workflow publishes a wheel plus compiler/source/corpus fingerprint metadata under
`dist-pgo/`. Profiles are target- and compiler-specific, so no persistent PGO flags are
stored in `Cargo.toml`.

## Project contracts

- [Detailed architecture and behavior](DESIGN.md)
- [Implementation and review standard](ENGINEERING.md)

## License

Nibbler is released under the [MIT License](LICENSE).
