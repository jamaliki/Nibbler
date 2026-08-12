# Nibbler engineering standard

- Status: current
- Updated: 2026-08-12
- Applies to: Rust, Python, tests, benchmarks, schema artifacts, and build tooling

This is the implementation and review standard for [DESIGN.md](DESIGN.md). Correctness,
performance, and readability are independent release gates; success in one does not
excuse regression in another.

## 1. Reviewability

A reviewer familiar with Rust or Python must be able to determine:

1. which invariant a change preserves or introduces;
2. which bytes, rows, or schema objects it reads and owns;
3. how missing values and failures propagate;
4. why a non-obvious optimization exists; and
5. which test or benchmark detects a regression.

Prefer the smallest implementation that satisfies the measured requirement. Line count
is useful pressure, not proof of simplicity: explicit invariants and clear ownership are
worth lines; forwarding layers, duplicate semantics, and speculative abstractions are
not.

## 2. Architecture

### One production path

- Text CIF has one lexer and parser.
- `DocumentSink` and `TableSink` share that grammar owner.
- Canonical and preserving text output share syntax and value-formatting rules.
- BinaryCIF converges on the same document and table models.
- Python wrappers delegate to Rust and must not implement a second parser, validator,
  chemistry model, or writer.
- Reference algorithms may exist in tests, never as alternate production paths.

The dependency direction is:

```text
input -> lexer -> parser -> document/table sinks -> document or table
BinaryCIF codec -------------------------------> document or table

compiled dictionary -> projection typing and dictionary validation
document + optional local registry -> PDBx semantics -> ModelCIF semantics
document or semantic model -> writer -> destination

Rust API -> PyO3 conversion -> typed Python facade
```

Lower layers do not import chemistry, profile, or Python concerns. Side effects are
confined to input, compression, output, FFI, corpus fetching, and build workflows.
Parsing, validation, semantic construction, and writing perform no hidden filesystem or
network access.

### Boundaries that earn their place

- `SourceBuffer` owns immutable source text shared by spans.
- `ParseSink` separates grammar from document retention and projection.
- `CifDocument` and `CifTable` are the common text/BinaryCIF boundary.
- compiled dictionaries isolate DDL2 loading from validation and projection.
- PDBx and ModelCIF types name domain records and protect invariants.
- the PyO3 layer translates native values and failures once.

A new trait, generic framework, crate, registry, service layer, or callback system needs
at least two real implementations or a demonstrated ownership/performance boundary.

## 3. Names, values, and ownership

Use CIF and PDBx terms precisely: block, frame, category, tag, item, entity,
`label_asym_id`, `auth_seq_id`, component, and atom site. Do not replace distinct domain
identifiers with ambiguous names such as `id`, `obj`, or `data` when the distinction
matters. Include units or representation in names such as `byte_offset`, `row_count`,
and `source_span`.

The themed names `chomp`, `feast`, `sniff`, and `dump` belong only at the Python root.
Native modules, errors, configuration, and diagnostics use conventional domain names.

Unknown (`?`), not applicable (`.`), absent, and present values remain distinct. They
must not be encoded as empty strings, sentinel numbers, or NaN. Finite numeric values
and checked conversions are required at trust boundaries. Chemistry resolution carries
an explicit `ComponentResolution`; residue-name sets are not a general chemistry type
system.

Values are immutable after construction. Temporary mutation is local to an owning
builder. Shared ownership is explicit through `Arc` or immutable Python wrappers; no
global mutable model or cache is allowed. The one process-wide parallel-worker counter
is an atomic resource lease, not application state.

## 4. Rust rules

### Control flow and modules

- Prefer `match`, early returns, and small state machines over callbacks and hidden
  control flow.
- A function does one conceptual job. Fifty lines is a review prompt, not a mechanical
  limit; a longer linear parser transition may be clearer than fragmented helpers.
- A module has one reason to change. Five hundred lines is a decomposition prompt, not
  a reason to create one-function files.
- Extract shared code at real duplication or a strong domain boundary, not anticipated
  reuse.
- Keep hot-path helpers local when measurement shows that abstraction changes code
  generation.
- Macros are limited to small, obvious repetition or generated schema data.

### Failure policy

Library behavior must not panic because of input bytes, Python arguments, integer
overflow, allocation sizes, or destination failures.

- `unwrap`, `expect`, `panic!`, `todo!`, and `unimplemented!` are denied on reachable
  library paths.
- Direct indexing requires a locally visible bound or checked alternative.
- `debug_assert!` may document a programmer invariant, never validate input.
- Tests and fuzz targets may use assertions and setup `expect` calls.
- Errors carry stable codes and structured context. String matching must not drive
  native control flow.

The crate sets `unsafe_code = "forbid"`. Production code must remain safe Rust.

### Documentation and comments

Public Rust APIs document ownership, missing states, failure modes, and useful examples.
Comments explain grammar rules, proof obligations, or counterintuitive performance
constraints; they do not paraphrase the next line. A performance-sensitive comment
names the guardrail or benchmark that justifies the shape.

## 5. Python rules

- The package is a typed facade, not a second implementation.
- Per-row Python loops are forbidden in native parse, validation, semantic-build, and
  write paths.
- Public functions and classes have complete annotations; `py.typed` and `_core.pyi`
  ship with the package.
- pandas, Polars, and PyArrow remain optional imports at the method that needs them.
- Public dispatch is explicit. No monkey patching, metaclass registry, or dynamic
  attribute protocol defines behavior.
- Validate external Python inputs once at the facade, then pass normalized plain values
  across PyO3.
- Translate native error tuples centrally into the public exception hierarchy.
- An obvious convenience wrapper should fit on one screen; complex behavior belongs in
  a named domain function with direct tests.

Public APIs may change before 1.0 when the result is materially simpler. Such a change
updates exports, stubs, callers, examples, tests, and documentation in the same patch.
Do not keep aliases, forwarding modules, or deprecated wrappers unless compatibility is
an explicit release requirement.

## 6. Dependencies

Every runtime dependency must satisfy a current requirement, remove more code or risk
than it adds, have an acceptable license and maintenance posture, and preserve the
supported toolchain. Default features are disabled unless required.

Current native dependencies are:

| Dependency | Purpose |
| --- | --- |
| `arrow-array`, `arrow-schema` | typed column buffers and Arrow C Stream export |
| `flate2` with `zlib-rs` | portable gzip decoding |
| `memchr` | byte scanning in lexer and loop-boundary kernels |
| `pyo3` | optional Python extension boundary |
| `regex` | compiled DDL2 type-pattern validation |
| `rmp-serde`, `serde`, `serde_bytes` | BinaryCIF MessagePack model |
| `proptest` (development) | property tests |

There is one Rust crate and one Python package. Cargo and Python lock inputs must remain
reproducible. Runtime dependencies may not introduce network access.

## 7. Performance work

Performance changes follow a scientific loop:

1. name the end-to-end workload and correctness guardrails;
2. record a reproducible comparison on pinned inputs;
3. profile before changing architecture;
4. implement the smallest credible mechanism;
5. compare interleaved samples and peak RSS;
6. run semantic, determinism, and adversarial tests; and
7. keep the change only when the representative result is neutral or better.

The unit of success is a usable `CifDocument`, `CifTable`, Arrow table, or semantic
model, not a token-count microbenchmark. Report logical decompressed throughput for
gzip and source-byte throughput for CIF/BinaryCIF. Record hardware, compiler, build
profile, corpus digest, warmups, samples, latency distribution, throughput, RSS, and
output digest.

Optimizations must preserve strict grammar, stable errors, resource-limit ordering,
source order, worker-count determinism, and small-file latency. Specialized probes may
decline work and fall back to the grammar owner; they must never redefine validity.

Current mechanisms and measurements live in
[docs/performance-report.md](docs/performance-report.md). Experimental history is
isolated in [docs/improvement-beam.md](docs/improvement-beam.md), not mixed into current
architecture documents.

## 8. Tests

Tests are organized by invariant rather than implementation detail:

- Rust unit and integration tests cover grammar, errors, limits, projection,
  dictionaries, BinaryCIF, semantic models, and writers.
- Python tests cover the public facade, optional interchange, transactional I/O,
  deterministic scanning, error conversion, and typing contracts.
- manifest tests pin every fixture by size and SHA-256.
- interoperability tests qualify outputs with installed external readers/validators.
- fuzz targets exercise the lexer, parse/write/parse equivalence, and formatter.
- the PDB stress corpus exercises small structures, chemistry, multiple models, a
  ribosome, and a multi-million-atom assembly in CIF, gzip, and BinaryCIF.

Every bug fix adds the smallest failing case at the layer that owns the invariant. A
performance path has a serial or general-path equivalence test. Tests do not depend on
network access; fetch commands populate ignored caches separately.

## 9. Required checks

The repository-local micromamba environment is the only supported development runtime.
The `Makefile` is authoritative:

```console
make bootstrap          # create .mamba/nibbler-dev
make develop            # build the development extension
make format             # Rust and Python formatting
make lint               # rustfmt, Clippy, Ruff
make typecheck          # strict mypy
make test               # Rust and Python tests
make docs               # rustdoc with warnings denied
make corpus             # verify pinned local fixtures
make validate-fixtures  # external fixture qualification
make check              # all gates plus benchmark smoke test
```

Direct equivalents use the same environment and flags:

```console
micromamba run -p .mamba/nibbler-dev cargo fmt --all -- --check
micromamba run -p .mamba/nibbler-dev cargo clippy --all-targets \
  --no-default-features -- -D warnings
micromamba run -p .mamba/nibbler-dev cargo test --no-default-features
micromamba run -p .mamba/nibbler-dev pytest
micromamba run -p .mamba/nibbler-dev ruff format --check .
micromamba run -p .mamba/nibbler-dev ruff check .
micromamba run -p .mamba/nibbler-dev mypy --strict python benchmarks tools
micromamba run -p .mamba/nibbler-dev env RUSTDOCFLAGS=-D warnings \
  cargo doc --no-deps --no-default-features
micromamba run -p .mamba/nibbler-dev python -m tools.verify_corpus
```

Release benchmarks require `maturin develop --release`. Native PGO wheels are built
per target and Python ABI with:

```console
micromamba run -p .mamba/nibbler-dev python -m tools.build_pgo
```

## 10. Change policy and definition of done

A change is complete only when:

- the design has one clear owner per invariant;
- obsolete code and compatibility layers are removed;
- public types, stubs, docstrings, examples, and design documents agree;
- formatting, lint, typing, tests, rustdoc, corpus verification, and affected
  interoperability checks pass;
- performance-sensitive changes include representative A/B evidence; and
- the diff has no generated caches, build artifacts, fetched corpora, or secrets.

Reviews check correctness and also ask whether a smaller, flatter implementation would
make local reasoning easier. The chosen design must be the smallest sufficient design;
historical alternatives belong only in the experiment ledger.
