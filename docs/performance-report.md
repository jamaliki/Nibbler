# Performance architecture and qualification

- Status: current implementation
- Qualified: 2026-08-12
- Reference host: 16-core Apple M4 Max
- Corpus: hash-pinned PDB CIF, gzip, and BinaryCIF inputs in
  [`benchmarks/pdb_corpus.toml`](../benchmarks/pdb_corpus.toml)

This document describes the performance design present in the repository and its
qualified checkpoint. Historical experiments are recorded separately in
[improvement-beam.md](improvement-beam.md).

## 1. Performance contract

Nibbler measures end-to-end construction of a usable projected table or complete
logical document. A timing is publishable only when:

- strict CIF behavior and stable error locations are unchanged;
- projected Arrow tables or canonical document bytes are exactly equal;
- resource limits and the earliest error are stable across worker counts;
- concurrency is bounded;
- small-file latency and peak RSS remain acceptable; and
- the release extension, input digests, compiler, hardware, warmups, and samples are
  recorded.

MB/s uses source bytes for CIF and BinaryCIF and logical decompressed CIF bytes for
gzip. See [../benchmarks/README.md](../benchmarks/README.md) for commands and report
fields.

## 2. Qualified checkpoint

Medians below use the corpus-trained, fingerprinted PGO wheel on the reference host.
Large projection ranges are two independent runs with two warmups and nine samples;
full documents use one warmup and three samples in isolated processes. They are
workload-specific measurements, not portable thresholds.

Measured artifact and toolchain:

- package: `nibbler-cif 0.1.0`;
- wheel: `nibbler_cif-0.1.0-cp312-cp312-macosx_11_0_arm64.whl`, SHA-256
  `222b0f3f2fb44075c83bbee30081be13df0a364cda81c2bb9407595418f9605d`;
- native extension SHA-256:
  `fb27ea5be2034ad9cbf68c28bdaf1af2788de617fd4632ea6e7348cf17b3d9b0`;
- Python/platform: CPython 3.12.13 on `macOS-26.5.2-arm64-arm-64bit`;
- compiler: `rustc 1.97.1` commit `8bab26f4f68e0e26f0bb7960be334d5b520ea452`,
  LLVM 22.1.6, target `aarch64-apple-darwin`; and
- PGO input fingerprint:
  `865099cece8e32f4a8796045eebdd8959badfe91b16fd4dcefd9ac27d46ff963`.

| Workload | PDB | Throughput |
| --- | --- | ---: |
| text projection | `6qnr` | 1,393-1,486 MB/s |
| text projection | `3j3q` | 1,771-1,869 MB/s |
| text full document | `6qnr` | 950 MB/s |
| text full document | `3j3q` | 1,422 MB/s |
| BinaryCIF projection | `6qnr` | 1,774-1,917 MB/s |
| BinaryCIF projection | `3j3q` | 1,755-1,787 MB/s |
| BinaryCIF full document | `6qnr` | 716 MB/s |
| BinaryCIF full document | `3j3q` | 839 MB/s |
| gzip text projection | `6qnr` | 646-653 MB/s logical CIF |
| gzip text projection | `3j3q` | 750-766 MB/s logical CIF |

Small `1crn` latency is 0.31-0.32 ms for text projection, 0.30 ms for a text document,
0.70-0.72 ms for BinaryCIF projection, and 0.99 ms for a BinaryCIF document.

Peak RSS for full documents is:

| Format | `6qnr` | `3j3q` |
| --- | ---: | ---: |
| text CIF | 264.5 MiB | 1,438.2 MiB |
| BinaryCIF | 185.3 MiB | 878.2 MiB |

Arrow import of the already projected table takes 0.6-0.8 ms for `6qnr` and 5.0-5.2 ms
for `3j3q`.

## 3. Text representation

`SourceBuffer` owns one immutable UTF-8 allocation behind `Arc`. Raw and decompressed
owned input moves into that allocation rather than being copied into a second source
string. A `CifDocument` owns its block slice behind another `Arc`, so semantic models
can retain the source document without duplicating it.

Text loop values are 8-byte source descriptors. The source bytes determine text bounds,
quote style, and missing kind lazily. Parallel document construction stores descriptors
in 65,536-cell segments that become final document storage directly; it does not join
all cells into a second contiguous vector.

BinaryCIF keeps decoded typed/dictionary columns. It does not expand every cell into a
row-major enum.

Semantic construction uses the same borrowed category-occurrence view as dictionary
validation. A decoder resolves its fixed set of columns once per occurrence and reads
rows lazily, so the `_atom_site` path creates neither a map nor an intermediate row
object per atom.

## 4. Lexer

Ordinary unquoted values use one scalar pass for delimiter discovery and forbidden
character validation. Quoted values, semicolon text, comments, control words, UTF-8,
and token limits stay in the same strict production lexer.

The parser is the only owner of block, frame, tag, and loop grammar. Performance kernels
may find candidate loop ranges or materialize values, but they cannot accept input that
the parser would reject.

## 5. Bounded intra-file parallelism

When at least 8 MiB remains from a loop's first value, the loop may request
`ceil(remaining bytes / 2 MiB)` workers, capped by available parallelism. Document mode
first proves that the loop itself reaches the threshold. A process-wide atomic lease
makes concurrent intra-file loop kernels share the host's available parallelism. The
outer multi-file scan pool is bounded independently. Fewer than two loop workers means
serial loop execution.

The large-loop path is:

1. **Eligibility.** Document mode uses an allocation-free byte probe to prove that the
   current loop reaches the 8 MiB threshold. Ambiguous or malformed input declines the
   fast path.
2. **Safe boundaries.** Workers summarize line-leading semicolons; the coordinator
   composes the exact two-state text-field state machine and produces ranges that never
   split a semicolon value.
3. **Validation.** Every range runs the production lexer, counts values, records the
   last value/control token, and reports a located error.
4. **Document retention.** Validation workers store compact source cells directly in
   final source-ordered segments.
5. **Projection.** Workers validate and build typed selected columns in one pass while
   assuming local column phase zero. Global prefix counts prove every chunk boundary is
   row-aligned before output or typed errors are committed. If the proof fails, Nibbler
   discards speculative output and runs the general aligned projection path.
6. **Commit.** The coordinator selects the earliest error before the loop control token,
   checks row/value limits, and exposes chunks in source order without copying payload
   buffers.

Custom limits that could change observable error ordering force serial loop processing.
The safe-boundary index and speculative commit rules have adversarial serial/parallel
equivalence tests for comments, quotes, multiline text, row wrapping, control-looking
values, incomplete rows, predicates, and lexical/resource errors.

## 6. Projection and Arrow

Projection plans normalize requested items and compile predicate column indexes once
per category occurrence. Unselected values are validated but not retained. A compiled
dictionary selects text, `int64`, or `float64` builders; schema-less values remain text.

Columns retain parallel Arrow chunks rather than concatenating them. Missing-kind counts
are computed once per column, and Arrow import reuses the segmented buffers through the
C Stream interface. The public missing policy is applied only at export.

## 7. Gzip and input

`flate2` uses the portable pure-Rust `zlib-rs` backend. Decompression enforces input,
output, and expansion-ratio limits. Its initial output reserve is bounded by:

- five times compressed size;
- the configured decompressed limit;
- the configured expansion-ratio limit; and
- 256 MiB.

The reserve reduces reallocations without trusting the compressed stream. Text,
BinaryCIF, and gzip are detected by content, and all feed the same document/projection
contracts.

## 8. Optional profile-guided native build

`tools/build_pgo.py` builds one host-native optimization artifact. It is separate from
the portable release workflow. The tool:

1. rejects inherited profile flags and non-host targets;
2. builds an instrumented wheel;
3. verifies the pinned corpus by size and SHA-256;
4. trains text, BinaryCIF, gzip, projection, and full-document paths while checking
   output counts;
5. merges data with the active Rust toolchain's `llvm-profdata`;
6. rebuilds the optimized wheel; and
7. writes the ABI-specific wheel plus compiler, target, source, and corpus fingerprint
   metadata to `dist-pgo/`.

```console
micromamba run -p .mamba/nibbler-dev python -m tools.build_pgo
```

Profiles are never checked into Cargo configuration because they are coupled to the
compiler, target, ABI, source, and training workload. Portable release wheels are built
without PGO on the full supported platform and CPython matrix.

## 9. Hardware interpretation

On the reference host, warm file reads measure 12-24 GB/s and a payload copy measures
about 169 GB/s, while strict text parsing reaches 1-2 GB/s. Modeled document traffic is
only about 4 GB/s. The parser is therefore compute- and control-flow-bound, not DRAM- or
storage-bandwidth-bound.

End-to-end average utilization varies with serial metadata and output ownership, from
roughly 2.7 to 7.4 cores on the large corpus cases. The current design is near the
practical throughput limit of its strict generic lexer and logical-output contract, not
the machine's literal memory-copy limit.

## 10. Qualification commands

The current implementation is gated by:

```console
make check
micromamba run -p .mamba/nibbler-dev python -m tools.build_pgo
PGO_WHEEL=$(find dist-pgo -name '*.whl' -print -quit)
micromamba run -p .mamba/nibbler-dev python -m tools.qualify_release "$PGO_WHEEL"
python -m venv .cache/pgo-benchmark-env
.cache/pgo-benchmark-env/bin/python -m pip install "$PGO_WHEEL" pyarrow
.cache/pgo-benchmark-env/bin/python -m benchmarks.pdb_stress \
  --ids 1crn 6qnr 3j3q --formats all --warmups 2 --samples 9 --json
# Run the projection command twice.
.cache/pgo-benchmark-env/bin/python -m benchmarks.pdb_stress \
  --ids 1crn 6qnr 3j3q --formats both --full-document --warmups 1 --samples 3 --json
```

Qualification covers the Rust and Python suites, Clippy, rustdoc, Ruff, strict mypy,
manifest verification, equal typed projections across all three input encodings, and
canonical parse/write/parse equality for complete documents.
