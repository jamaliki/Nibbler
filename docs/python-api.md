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
a BinaryCIF category with no rows counts as absent, as it does in the decoded document.
A matching category occurrence missing a requested or predicate item is an error.

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

## Biological assemblies

`mmcif.assembly` writes out every copy of one biological assembly as explicit
coordinates, so a consumer of the coordinates sees the contacts between copies:

```python
document = nibbler.mmcif.assembly(source_document_or_model, assembly_id=None)
result = nibbler.reduce.run(document)  # hydrogens see the other copies
```

It accepts a `CifDocument`, an `MmcifModel`, or any source `cif.read` accepts (built
with `registry=`, as `mmcif.read` does), and returns a new `CifDocument` with one block.
`assembly_id` selects a `_pdbx_struct_assembly`; the default is the first. Each distinct
operator combination of the assembly's `_pdbx_struct_assembly_gen` rows is one copy,
numbered from 1 in order of first appearance. An operator expression is a list (`1,2`,
`1-4`), a parenthesized list, or a product of lists (`(1-60)(61)`) whose combinations
apply the rightmost operator first. Copy 1 keeps the chain identifiers; copy *n* appends
`-n` to the `label_asym_id` and `auth_asym_id` of its chains (`A` becomes `A-2`). A copy of
an atom within 0.2 Å of an earlier copy of the same atom (an atom on a symmetry axis of
the assembly) is written once, as gemmi's `transform_to_assembly` does. A chain copy all
of whose atoms are written that way (an ion on an axis) is left out, and connections to
it name the earlier copy instead.

The result keeps the entry, entity, and component categories and `_atom_type`
unchanged, except that `_entity_poly.pdbx_strand_id` lists every copy's chains. It has
one row per copy in `_struct_asym`, the sequence schemes, `_atom_site` (coordinates
transformed and written with three decimals, atom ids renumbered, model by model),
`_atom_site_anisotrop` (tensors rotated; uncertainties of a rotated tensor become `?`),
and `_struct_conn` (only connections within one crystal copy, symmetry `1_555`). Every
other category is left out, including the crystal cell and symmetry, so nothing
downstream applies the crystal symmetry again, and the assembly definitions, which
described the source. The result builds a PDBx model and validates against the profile.

Failures raise `ChemistryError` with `PDBX_ASSEMBLY_ABSENT` (no assembly is defined),
`PDBX_ASSEMBLY_UNKNOWN`, `PDBX_ASSEMBLY_EXPRESSION`, `PDBX_ASSEMBLY_OPERATOR` (an
operator missing from `_pdbx_struct_oper_list`), `PDBX_ASSEMBLY_ASYM`,
`PDBX_ASSEMBLY_ID_COLLISION` (a renamed chain would take an existing name),
`PDBX_ASSEMBLY_TOO_LARGE` (over 50 million atom sites), or a `PDBX_ITEM_*` code for a
missing or untyped operator or coordinate value.

## Hydrogens with Reduce3

Builds with the `reduce3` Cargo feature, which the Python package enables, can run
[Reduce3](https://github.com/jamaliki/reduce3) on a model without leaving memory:

```python
result = nibbler.reduce.run(
    source_document_or_model,
    approach="add",
    add_flip_movers=True,
    compat=False,
    chem_data=None,
)
result.document  # a new CifDocument
result.report  # Reduce3's description text
```

`run` accepts a `CifDocument`, an `MmcifModel` (its source document is used), or any
source `cif.read` accepts. Reduce3 reads the first data block with an `_atom_site` loop
directly from the parsed values, adds and places hydrogens, optimizes rotatable and
flippable groups, and builds the result as a new `CifDocument`. No text is written or
parsed again, and the GIL is released while it runs.

The result block keeps the source block code and the source document's schema selector.
By default it also keeps every category of the source block, in order: `_atom_site` is
rebuilt with the source's items and label identifiers (new hydrogens take their
residue's), atom ids are renumbered, `_atom_site_anisotrop` follows the new ids, and
`_atom_type` gains the elements it lacks. Everything else, including `_struct_conn`, the
entities and the sequence schemes, is shared with the source unchanged. With
`compat=True` the block has the layout Reduce2 writes instead: cell, space group,
`_struct_asym`, `_chem_comp`, `_atom_site`, and `_atom_site_anisotrop`, with regenerated
label identifiers. Either way the entries equal a parse of the file the `reduce3`
program writes in that mode.

The keyword options are Reduce2's parameters with Reduce2's defaults, except that
`add_flip_movers` is True (Reduce2: False), so Asn, Gln and His flips are considered
unless it is set to False: `approach` (`"add"`, `"remove"`, or `"optimize"`),
`add_flip_movers`, `n_terminal_charge`, `keep_existing_h`, `exclude_water`,
`use_neutron_distances`, `preference_magnitude`, `non_flip_preference`,
`planar_hydroxyl_preference`, `acid_syn_preference`, `skip_bond_fix_up`,
`set_flip_states`, `model_id`, `alt_id`,
`bonded_neighbor_depth`, `stop_on_any_missing_hydrogen`, `ignore_missing_restraints`,
`verbosity`, and `probe`, a mapping of Probe scoring parameters such as
`{"probe_radius": 0.25}`. `compat=True` reproduces Reduce2 exactly, including its known
defects; the default corrects them and goes on where Reduce2 gives up. Residues that
neither monomer library describes get restraints built from their chemical component
definition, including the entries Reduce2's RDKit step rejects, with GeoStd-style atom
types (so their donors and acceptors take part in scoring) and physiological protonation
of their acids; a hydroxyl hydrogen on a planar atom (a phenol such as tyrosine, an enol,
a carboxylic acid) prefers that atom's plane by `planar_hydroxyl_preference` and an
acid's hydrogen prefers syn to its carbonyl oxygen by `acid_syn_preference` (both 1.0
Probe score units by default, 0 turns them off; Reduce2 has neither); atoms of unknown element
(`UNX`) are kept unchanged; and a residue with no definition anywhere keeps its input
hydrogens and is reported in the report instead of stopping the run
(`stop_on_any_missing_hydrogen=True` stops it).

Reduce3 needs the cctbx `chem_data` monomer library. `chem_data` names its directory;
otherwise it is found through `REDUCE3_CHEM_DATA`, `CHEM_DATA`, the parent of
`MMTBX_CCP4_MONOMER_LIB` or `CLIBD_MON`, or the active conda environment. It is loaded
once per directory and process. Failures raise `ChemistryError` with a `REDUCE_*` code:
`REDUCE_MONOMER_LIBRARY_NOT_FOUND`, `REDUCE_MONOMER_LIBRARY_INVALID`,
`REDUCE_INVALID_MODEL`, `REDUCE_FAILED` (for example, missing restraints in compat mode),
or
`REDUCE_RESULT_TOO_LARGE`. `nibbler.reduce.available()` reports whether the build
includes Reduce3.

## Errors

All public failures derive from `NibblerError`:

- `ParseError`: input, decompression, UTF-8, resource, or CIF syntax failure;
- `ProjectionError`: malformed projection or incompatible category occurrence;
- `SchemaError`: schema selector, artifact, dictionary item, or typed-value failure;
- `ChemistryError`: semantic-model construction, component-resolution, or Reduce3
  failure;
- `WriteError`: serialization or destination failure;
- `BatchError`: one source failure from `scan`, with `source_index`; and
- `ValidationError`: one or more error diagnostics from a validation report.

Parse-derived errors expose a stable `code`, message, source name, byte span, and
one-based line and display column when available. Code and structured fields are the
machine contract; human-readable messages may add detail.
