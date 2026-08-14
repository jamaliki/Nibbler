# Python API reference

Nibbler's import package is `nibbler`. The four short root names are a facade over
conventional modules:

| Root API | Conventional API | Purpose |
| --- | --- | --- |
| `nibbler.chomp` | `nibbler.cif.read` | read one document or project one category |
| `nibbler.feast` | `nibbler.cif.scan` | read many sources through a bounded worker pool |
| `nibbler.sniff` | `nibbler.cif.validate` or `nibbler.mmcif.validate` | validate a document or semantic model |
| `nibbler.dump` | `nibbler.cif.write` or `nibbler.mmcif.write` | write a complete document or semantic model |

`chomp` and `feast` are exact aliases. `sniff` and `dump` dispatch between generic CIF
and semantic APIs from `schema=`, `profile=`, and the supplied Nibbler object.

## Sources and formats

`cif.read`, `cif.scan`, `mmcif.read`, and `components.read` accept:

- a `str` or string-valued `os.PathLike` filesystem path;
- `bytes`, `bytearray`, or `memoryview`; or
- a binary file object whose `read()` method returns bytes-like data.

A binary file object is read in full. Text CIF, gzip, and BinaryCIF 0.3 are detected
from content, not from a filename suffix. Gzip is decompressed before text or BinaryCIF
detection and is subject to the input, decompressed-size, and expansion-ratio limits in
[DESIGN.md](../DESIGN.md#12-resource-and-safety-contract).

## Read one source

```python
nibbler.cif.read(
    source,
    *,
    category=None,
    columns=None,
    where=None,
    schema=None,
)
```

With no `category`, `read` returns an immutable `CifDocument`. `columns` and `where`
are invalid in this mode. Passing `schema="pdbx"` or `schema="modelcif"` attaches that
dictionary selector to the document for later validation and writing; it does not by
itself validate every dictionary constraint. `schema="mmcif"` is accepted as an alias
for `"pdbx"`.

With `category`, `read` returns a `CifTable`:

```python
atoms = nibbler.cif.read(
    "structure.cif.gz",
    category="atom_site",
    columns=["id", "label_comp_id", "Cartn_x"],
    where={"pdbx_PDB_model_num": 1},
    schema="pdbx",
)
```

Category and item matching is ASCII case-insensitive. A category may be written with or
without its leading underscore. Item selectors may be bare names such as `Cartn_x` or
full tags from the selected category. Explicit output columns retain caller order.
Omitting `columns` retains the category's source item order. Predicates may refer to an
item that is not returned. A category absent from the source produces a zero-row table;
a matching category occurrence missing a requested or predicate item is an error.

Without a schema, projected columns are text. A schema-guided projection rejects
unknown items and builds dictionary-declared text, `int64`, or `float64` columns.
Present numeric values must parse and be finite.

### Predicates

`where` maps item names to native, text-based conditions:

```python
from nibbler import MissingKind

where = {
    "label_comp_id": {"ALA", "GLY"},  # membership
    "label_asym_id": ("!=", "B"),  # inequality
    "pdbx_PDB_model_num": 1,  # equality after str(1)
    "occupancy": None,  # either CIF missing state
    "Cartn_x": MissingKind.UNKNOWN,  # exactly ?
    "Cartn_y": ("is", MissingKind.NOT_APPLICABLE),  # exactly .
}
```

The tuple operators are `==`/`eq`, `!=`/`ne`, `in`, and `is`/`missing`. A non-string
collection is shorthand for membership; every other value is equality. Present-value
comparisons use the decoded CIF text exactly and do not run Python callbacks in the row
loop. Missing values match only a missing-state predicate.

## Documents and tables

`CifDocument` exposes:

- `block_count`: number of data and global blocks;
- `schema`: attached dictionary selector, or `None`;
- `to_canonical()`: deterministic CIF 1.1 text;
- `to_preserving()`: deterministic text using source order and valid original value
  lexemes; and
- `to_binary()`: deterministic BinaryCIF 0.3 bytes when the document is representable.

`CifTable` exposes `category`, `columns`, `len(table)`, `with_missing(policy)`, and the
Arrow C Stream protocol. PyArrow and Polars consume the object directly; neither is a
Nibbler runtime dependency.

```python
import pyarrow

arrow = pyarrow.RecordBatchReader.from_stream(atoms.with_missing("columns")).read_all()
```

The Arrow missing policies are:

| Policy | Representation |
| --- | --- |
| `collapse` | both CIF missing states become Arrow null; counts remain in field metadata |
| `columns` | each value column gets a non-null `uint8` `<name>__missing_kind` companion |
| `extension` | each value is a `struct<value, kind>` annotated as `nibbler.cif_missing` |

Missing-kind codes are `0` for present, `1` for unknown (`?`), and `2` for not
applicable (`.`). `collapse` is the default.

## Scan many sources

```python
nibbler.cif.scan(
    sources,
    *,
    category=None,
    columns=None,
    where=None,
    schema=None,
    workers=None,
    on_error="raise",
)
```

`scan` lazily admits at most one source per worker and yields results in input order.
The default worker count is `min(32, os.cpu_count() or 1)`; explicit counts must be
between 1 and 256. With `on_error="raise"`, the first failure in input order is raised
as `BatchError` and outstanding work is cancelled. With `on_error="collect"`, successful
results are yielded and failures accumulate in `result.errors`, also in input order.
There is no silent-skip mode.

Projected scan tables include `_nibbler_source`, `_nibbler_block`, and
`_nibbler_frame` provenance columns when exported through Arrow. `BatchDiagnostics`
supports sequence access and optional `to_pyarrow()` and `to_polars()` conversions. A
caller that stops iteration early should call `result.close()` to stop admission and
join native workers.

## Validation

`nibbler.cif.validate(document, schema=None)` accepts only `CifDocument`. With no
explicit or attached schema, successful parsing is reported as `cif-1.1-syntax`
coverage. With `pdbx` or `modelcif`, it applies the pinned DDL2 dictionary.

`nibbler.mmcif.validate(value, profile=...)` accepts `CifDocument` or `MmcifModel` and
requires `profile="pdbx"` or `profile="modelcif"`. A document receives dictionary and
semantic validation; an existing model receives semantic validation. A model cannot be
validated under a profile different from the one used to build it.

`ValidationReport` is immutable and contains `diagnostics`, `schema`,
`dictionary_version`, and explicit `coverage`. `is_valid` is true when no error
diagnostic exists; `errors` returns only errors, and `raise_for_errors()` raises
`ValidationError` carrying those diagnostics.

## Semantic models and components

```python
model = nibbler.mmcif.read(
    source_or_document,
    profile="pdbx",
    registry=None,
)
```

`MmcifModel` is immutable and source-backed. Common PDBx summaries are `profile`,
`entry_id`, `entity_count`, `asym_unit_count`, `component_count`, `atom_site_count`,
`connection_count`, `entity_kinds`, and `component_ids`. ModelCIF adds
`prediction_model_count`, `target_entity_count`, `template_count`, `qa_metric_count`,
`qa_value_count`, `software_names`, `qa_metric_names`, and `qa_metric_modes`.
`source_document()` returns the shared logical source document with its profile schema
attached.

`nibbler.components.read(source_or_document)` builds an immutable local component
`Registry`; `len(registry)` returns its component count. Supplying the registry to
`mmcif.read` never enables network access. Embedded definitions take priority, but a
conflicting embedded and local definition is an error rather than a silent override.

## Writing

Generic writing accepts only `CifDocument`:

```python
nibbler.cif.write(
    document,
    destination,
    *,
    mode="canonical",
    validate="syntax",
    format=None,
)
```

Semantic writing accepts only `MmcifModel` and requires its matching profile:

```python
nibbler.mmcif.write(
    model,
    destination,
    *,
    profile="pdbx",
    mode="canonical",
    validate="profile",
    mirror_local_qa_metric=None,
    format=None,
)
```

`destination` is a filesystem path or writable binary stream. Path suffixes select
BinaryCIF for `.bcif`/`.bcif.gz` and gzip for `.gz`; streams default to text CIF and
need `format="bcif"` for BinaryCIF. Streams are never implicitly gzip-compressed.

`mode="preserve"` is text-only and reuses source order and valid original value lexemes,
not comments or whitespace. Canonical semantic output orders profile categories and may
materialize resolved component definitions. ModelCIF may explicitly mirror one local QA
metric into B factors with `mirror_local_qa_metric=<integer metric ID>`; preserving mode
and non-local metrics reject that option.

Path writes use a temporary sibling file, flush and `fsync` it, optionally reparse and
validate it, then atomically replace the destination. A failed operation leaves an
existing destination unchanged. BinaryCIF cannot represent global blocks or save frames,
and its writer rejects values outside its lossless supported encodings.

## Errors

All public failures derive from `NibblerError`:

- `ParseError`: input, decompression, UTF-8, resource, or CIF syntax failure;
- `ProjectionError`: malformed projection or incompatible category occurrence;
- `SchemaError`: schema selector, artifact, dictionary item, or typed-value failure;
- `ChemistryError`: semantic-model construction or component-resolution failure;
- `WriteError`: serialization or destination failure;
- `BatchError`: one source failure from `scan`, with `source_index`; and
- `ValidationError`: one or more error diagnostics from a validation report.

Parse-derived errors expose a stable `code`, message, source name, byte span, and
one-based line and display column when available. Code and structured fields are the
machine contract; human-readable messages may add detail.
