# Nibbler Engineering Standard

- Status: Draft 0.1
- Date: 2026-08-11
- Applies to: Rust core, Python package, tests, benchmarks, generators, and build tooling

This is the implementation and review companion to [DESIGN.md](DESIGN.md). The words
MUST, MUST NOT, SHOULD, SHOULD NOT, and MAY are normative.

Nibbler does not trade readability for speed. Correctness, performance, and
maintainability are independent release gates: passing two does not compensate for
failing the third.

## 1. Reviewability contract

A change is reviewable when a contributor familiar with Rust or Python, but not with
its author, can determine:

1. which invariant it preserves or introduces;
2. which bytes, rows, or schema objects it reads and owns;
3. how errors and missing values propagate;
4. why any non-obvious optimization is necessary; and
5. which test would fail if the behavior regressed.

If those answers require reconstructing hidden control flow, macro expansion, global
state, or benchmark folklore, the change is not ready to merge.

The smallest implementation that meets the measured requirement is preferred. Line
count alone is not simplicity: duplicated semantics, implicit state, and unnecessary
abstractions are larger costs than a few explicit lines.

## 2. Architectural discipline

### 2.1 One production path

- There MUST be one tokenizer and parser for all production reads.
- Canonical and preserving output modes MUST share the same value formatter and syntax
  rules.
- Python convenience APIs MUST delegate to conventional submodule APIs rather than
  reimplement behavior.
- Reference implementations MAY exist in tests for differential testing. They MUST NOT
  become a second production path.
- A specialized kernel MAY replace an implementation detail only when it preserves the
  same contract and wins a representative benchmark.

### 2.2 Dependency direction

The intended dependency direction is:

```text
source -> lexer -> parser -> sinks -> document/tables
                                  -> diagnostics

dictionary -> schema validation -> mmCIF/ModelCIF semantics
document/tables + semantics -> writer -> destination
Rust public API -> PyO3 adapter -> Python facade
```

Lower layers MUST NOT import profile-specific chemistry or Python concerns. Domain
layers MUST consume parser output through narrow typed interfaces. Cyclic module or
package dependencies are forbidden.

Nibbler SHOULD begin as one Rust library crate plus one Python package. A new crate,
plugin system, service locator, registry, or framework requires evidence that the
existing structure is causing a concrete problem.

### 2.3 Functional boundaries

- Pure transformations are preferred for schema compilation, category planning,
  validation, and formatting decisions.
- I/O, allocation-heavy materialization, compression, and FFI MUST be kept behind
  narrow boundaries.
- Values SHOULD be immutable after construction. Mutation used for streaming table
  assembly MUST remain local to the owning builder.
- There MUST be no hidden filesystem or network access in parsing, validation, model
  construction, or writing.
- Global mutable state is forbidden. Read-only compiled tables MAY be initialized once
  through a thread-safe standard-library primitive.

## 3. Names and types

### 3.1 Domain language

Internal names MUST use CIF and PDBx terminology: `data_block`, `category`, `tag`,
`label_asym_id`, `auth_seq_id`, `entity`, and `component`. Do not shorten distinct
identifiers to ambiguous names such as `chain`, `residue`, `id`, `obj`, or `data` when
the distinction matters.

Names SHOULD expose units or representation where confusion is possible, for example
`byte_offset`, `row_count`, `compressed_size`, and `source_span`.

The character-themed API belongs only at the root facade:

- `nibbler.chomp(...)` may provide a friendly parse entry point.
- Conventional APIs such as `nibbler.cif.read(...)`, `nibbler.cif.write(...)`, and
  `nibbler.cif.validate(...)` remain the searchable, composable foundation.
- Internal types, error names, module names, configuration, and diagnostics MUST NOT
  use jokes or character lore.

### 3.2 Make invalid states difficult to represent

- Semantically distinct identifiers SHOULD use distinct types when accidental mixing
  would be plausible.
- CIF unknown (`?`), not-applicable (`.`), absent, and present values MUST remain
  distinct. They MUST NOT be encoded as magic strings, negative numbers, empty strings,
  or floating-point NaNs.
- Boolean parameters are forbidden when the call site would not reveal what `true`
  means. Use a named enum.
- Input sizes and byte offsets MUST use checked conversions and checked arithmetic.
- Chemistry classifications MUST carry provenance and confidence as specified in
  [DESIGN.md](DESIGN.md); a residue-name set is not a type system.

## 4. Rust implementation rules

### 4.1 Control flow

- Prefer explicit `match`, early returns, and small state machines over callback webs,
  deeply nested branches, or encoded bit tricks.
- Functions SHOULD do one conceptual job. Around 50 lines is a review prompt, not a
  mechanical limit; a longer linear parser transition can be clearer than fragmented
  helpers.
- Modules SHOULD have a single reason to change. A module approaching 500 lines should
  be reviewed for separable responsibilities, but MUST NOT be split into meaningless
  one-function files.
- Abstractions are introduced after concrete duplication or a demonstrated boundary,
  not in anticipation of reuse.
- A generic or trait MUST have at least two real implementations, define an important
  boundary, or deliver a measured static-dispatch benefit.
- Macros are reserved for generated schema tables and small, obvious repetition.
  Public declarative macros and procedural macros are out of scope initially.

### 4.2 Panic policy

Library code MUST NOT panic because of file contents, Python arguments, destination
failures, integer overflow, or allocation-size calculations.

- `unwrap`, `expect`, `panic!`, `todo!`, and `unimplemented!` are forbidden on reachable
  library input paths.
- Direct indexing requires a locally obvious bound or a checked alternative.
- Internal invariants SHOULD use typed construction. `debug_assert!` MAY document an
  invariant when violating it indicates a programmer defect, never malformed input.
- Test code MAY use `unwrap` or `expect` when it makes the asserted setup clearer.
- Fuzz targets MUST treat every panic, abort, excessive allocation, and nontermination
  as a defect.

Errors MUST be structured. Each diagnostic carries a stable code, severity, message,
source span when available, and causal context. Errors are translated at the PyO3
boundary once; string matching MUST NOT drive control flow.

### 4.3 Unsafe code

The initial workspace MUST compile with `unsafe_code = "forbid"`. If profiling later
proves that `unsafe` is required for FFI or a material hot-path improvement, relaxing
that rule requires a focused design review and all of the following:

1. a benchmark showing the representative benefit;
2. a safe implementation retained in the benchmark or differential test;
3. isolation in the smallest practical module and function;
4. a `SAFETY:` comment stating every caller and memory invariant;
5. targeted property, fuzz, sanitizer, and Miri coverage where applicable; and
6. an entry in the pull request's complexity ledger.

An unsafe optimization that is merely faster in a microbenchmark, or whose invariants
cannot be stated locally, MUST NOT merge.

### 4.4 Comments and documentation

Comments explain grammar decisions, invariants, lifetimes, or counterintuitive
performance findings. They do not translate the next line into English.

Public types and functions MUST document missing-value behavior, ownership or borrowing,
failure modes, and a minimal example. Non-obvious parser states and writer ordering
rules MUST link to the governing specification section or dictionary concept.

## 5. Python implementation rules

- The Python package is a typed facade over the native core, not a second parser or
  semantic engine.
- Per-row Python loops are forbidden in native read, validation, and write paths.
- Python MUST NOT duplicate component classification, identifier reconciliation, value
  quoting, or category-ordering logic.
- Public Python APIs MUST have complete type annotations and ship `py.typed` plus useful
  generated or maintained stubs for native symbols.
- Imports MUST remain light. pandas, Polars, PyArrow, NumPy, and compression libraries
  are optional integrations unless the core contract explicitly requires them.
- Dynamic attribute dispatch, metaclass registries, and monkey patching MUST NOT be used
  to define the public API.
- A convenience wrapper SHOULD normally fit on one screen. Complex behavior belongs in
  a named Rust or Python domain function with direct tests.

## 6. Dependency policy

Every runtime dependency must:

1. solve a current requirement;
2. remove more risk or maintained code than it introduces;
3. have an acceptable license and maintenance posture;
4. expose a narrower, clearer boundary than a local implementation; and
5. be tested with default features disabled unless those features are required.

The pull request adding a dependency MUST record its purpose, enabled features,
alternatives considered, and effect on wheel size and clean build time. Transitive
dependency count is a design signal, not a vanity metric.

No dataframe engine, query engine, async runtime, logging framework, parser generator,
or general plugin framework belongs in the core merely for convenience. Standard
library solutions are preferred when they remain direct and correct.

Dependencies used only for tests, fuzzing, generation, or benchmarks MUST remain out of
the runtime graph. Generated dictionary artifacts MUST record the generator version and
input hash and MUST be reproducible without network access.

## 7. Performance without obscurity

Performance work follows this loop:

1. name the user-visible workload and input distribution;
2. record a reproducible baseline;
3. identify a measured bottleneck;
4. state a mechanical hypothesis;
5. implement the smallest change that tests it;
6. run correctness, allocation, memory, and wall-time comparisons; and
7. keep the change only if the benefit is material for the named workload.

Benchmark evidence MUST include corpus hashes, command, hardware, software versions,
sample count, and variance or confidence information. Microbenchmarks explain a
mechanism; the Section 17 end-to-end workloads decide promotion.

The following require benchmark evidence and an explanatory comment or design note:

- custom SIMD or branchless parsing;
- manual buffering that replaces a standard primitive;
- caching with eviction or invalidation behavior;
- duplicated representations;
- specialized allocators;
- unsafe code; and
- platform-specific kernels.

Do not retain speculative fast paths. If an optimization does not materially improve
its target workload, delete it. Machine-specific code MUST sit behind a small safe
interface and have a portable implementation tested for semantic equality.

## 8. Test design

- Tests assert public behavior, format invariants, or important internal state-machine
  transitions, not private call sequences.
- Every parser or writer bug fix begins with the smallest durable regression fixture.
- Table-driven tests are preferred for grammar and quoting boundaries.
- Property tests cover token boundaries, missing states, integer/float formatting,
  projection equivalence, and parse-write-parse invariants.
- Fuzzing covers tokenization, parsing, validation, and parse-write-parse composition.
- Differential tests compare trusted implementations where contracts overlap, but
  Nibbler's specified behavior remains authoritative.
- Golden files are appropriate for canonical output and diagnostics. Reviewers MUST be
  able to inspect the semantic change; blindly regenerated snapshots cannot merge.
- Randomized tests MUST record seeds and minimize failures into committed fixtures.
- Tests MUST NOT depend on the network, wall-clock ordering, locale, thread schedule, or
  mutable global state.

Ligands, ions, waters, branched carbohydrates, modified residues, alternate locations,
negative sequence identifiers, insertion codes, and unknown components are first-class
test dimensions, not an edge-case bucket.

## 9. Tooling and continuous integration

The bootstrap implementation SHOULD establish these gates before optimization begins:

```text
cargo fmt --check
cargo clippy --workspace --all-targets --all-features -- -D warnings
cargo test --workspace --all-features
RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps
ruff format --check .
ruff check .
mypy --strict python
pytest
```

Exact directories MAY change with the scaffold, but equivalent checks MUST remain.
Compiler and linter exceptions require a local explanation; workspace-wide suppression
is forbidden for a local issue.

Pull-request CI MUST include supported-platform unit tests and Python wheel smoke tests.
Scheduled CI SHOULD run the malformed corpus, fuzzers, Miri, sanitizers, dependency
audits, and the reproducible benchmark suite. Expensive scheduled checks do not replace
targeted regression tests on pull requests.

Performance regression thresholds MUST be based on a stable baseline and noise study.
A noisy benchmark must be repaired, not converted into a meaningless hard gate.

## 10. Change and review policy

Each change MUST be focused. Unrelated cleanup, generated artifacts, vendored data, API
changes, and performance work should be separate unless splitting would make the
correctness argument harder.

A pull request description records:

- the problem and affected invariant;
- the chosen approach and the simpler alternatives rejected;
- new public API or compatibility effects;
- correctness tests and validators run;
- benchmark results when hot code changes;
- new dependencies, unsafe code, caches, generated code, or platform specialization;
  and
- remaining risk.

The last-but-one item is the **complexity ledger**. It is empty for an ordinary change.
Each entry must state why the complexity is earned, where its invariant is tested, who
owns it, and what evidence would justify deleting it.

Reviewers MUST reject a change when:

- behavior is duplicated across parser, writer, Python, or profile layers;
- malformed input can panic or cause an uncontrolled allocation;
- chemistry is inferred silently from an incomplete name list;
- a new abstraction has no present use;
- an optimization lacks representative evidence;
- a dependency substitutes a framework for a small explicit boundary;
- diagnostics discard location or causal information;
- tests encode implementation accidents rather than the contract; or
- the author cannot state the invariant locally.

Before 1.0, obsolete internal and provisional API paths SHOULD be deleted rather than
wrapped in compatibility layers. Stable public compatibility is deliberate and tested;
dead code, commented-out implementations, and indefinite deprecations are forbidden.

## 11. Definition of done

A production change is done only when:

1. its contract and names are clear without oral context;
2. malformed and boundary inputs have explicit behavior;
3. tests cover the success, missing-value, and failure paths;
4. format, lint, type, documentation, and relevant sanitizer checks pass;
5. affected CIF and profile fixtures revalidate;
6. performance claims are reproducible and correctness-equivalent;
7. public documentation and examples match the implementation; and
8. the complexity ledger is empty or every entry is justified.

The desired result is intentionally boring code on a fast path: explicit state,
predictable ownership, narrow interfaces, and enough benchmark evidence that no reader
has to trust cleverness.
