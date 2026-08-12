# Nibbler performance search report

## Hardware-limit cycle

The second cycle eliminated the generic second projection scan on row-aligned chunks,
parallelized exact safe-boundary discovery, made full-document retention single-pass,
reduced each text cell descriptor from 12 to 8 bytes, replaced a serial generic-lexer
eligibility pass with a conservative allocation-free byte probe, and removed the final
source-sized input copy. It used fresh current-tree baselines because the installed
development extension was not byte-identical to the earlier published PGO artifact.

## Scope and guardrails

The Phase 6 optimization search targeted end-to-end construction of usable projected
tables and complete documents, not isolated token counters. Tests used the hash-pinned
PDB corpus in `.cache/pdb-stress` on a 16-core Apple M4 Max. Strict CIF semantics,
deterministic source order and errors, configured resource limits, bounded concurrency,
small-file latency, and `unsafe_code = "forbid"` were hard gates.

Baseline and final JSON artifacts are under `.cache/performance/`. Representative
commands are:

```console
micromamba run -p .mamba/nibbler-dev maturin develop --release
micromamba run -p .mamba/nibbler-dev python -m benchmarks.pdb_stress \
  --formats all --warmups 2 --samples 7 --json
micromamba run -p .mamba/nibbler-dev python -m benchmarks.pdb_stress \
  --ids 6qnr 3j3q --formats cif --full-document --warmups 1 --samples 3 --json
micromamba run -p .mamba/nibbler-dev python tools/build_pgo.py
```

## Results

Large-file medians compare the initial Phase 6 release extension with the final
corpus-trained native PGO wheel. MB/s uses source bytes for CIF/BinaryCIF and logical
decompressed CIF bytes for gzip.

| Workload | PDB | Baseline | Final | Change |
|---|---|---:|---:|---:|
| text projection | 6qnr | 321 MB/s | 1,603 MB/s | 4.99x |
| text projection | 3j3q | 278 MB/s | 2,092 MB/s | 7.53x |
| text full document | 6qnr | 272 MB/s | 969 MB/s | 3.56x |
| text full document | 3j3q | 243 MB/s | 1,512 MB/s | 6.22x |
| BinaryCIF projection | 6qnr | 1,269 MB/s | 2,068 MB/s | 1.63x |
| BinaryCIF projection | 3j3q | 1,303 MB/s | 1,971 MB/s | 1.51x |
| BinaryCIF full document | 6qnr | 188 MB/s | 811 MB/s | 4.32x |
| BinaryCIF full document | 3j3q | 206 MB/s | 874 MB/s | 4.24x |
| gzip text projection | 6qnr | - | 694 MB/s | - |
| gzip text projection | 3j3q | - | 808 MB/s | - |

Final full-document peak RSS was 405 MB for 6qnr text, 1.90 GB for 3j3q text,
279 MB for 6qnr BinaryCIF, and 1.12 GB for 3j3q BinaryCIF. Relative to the initial
BinaryCIF full-document representation, the large cases use 55-77% less peak memory.
Small 1crn text projection/full latency is 0.30/0.28 ms and BinaryCIF is
0.67/0.98 ms.

Against the final same-source and same-compiler plain wheel, PGO improves 6qnr/3j3q
text projection by 16.8/18.1% and text full construction by 13.2/13.7%.

Arrow export originally repeated whole-column missing-state counts for every parallel
chunk. Hoisting those counts reduced Python Arrow import from 9.7 to 0.7 ms on 6qnr and
from 73.8 to 4.9 ms on 3j3q.

## Promoted mechanisms

### Compact documents

Text loops store 8-byte source token spans with one shared source allocation. Quote style,
missing kind, and content bounds are derived lazily from the immutable source. Parallel
text documents own segmented cell chunks directly in the validation pass, avoiding both
a second scan and a concatenation copy. BinaryCIF
documents retain decoded typed/dictionary columns rather than expanding row-major value
enums. Public APIs expose borrowed values and rows.

### Fused strict lexer

Ordinary unquoted values use one scalar pass for delimiter discovery and forbidden
character validation. Real 6qnr tokens average 2.95 bytes and 99.5% are shorter than
eight bytes, so scalar fusion beat speculative long-token SIMD. Isolated gains were
8.7-10.2% on the large text files and neutral on the small file.

### Bounded intra-file parallelism

The serial parser remains the sole block/frame/tag grammar owner. Loops at least 8 MiB
use up to `min(available cores, ceil(loop tail / 2 MiB))` workers under a process-wide
atomic budget:

1. summarize column-one semicolons in parallel and compose the exact two-state text-field
   DFA to find safe boundaries;
2. tokenize, validate, classify, and count each range with the production lexer;
3. for projection, speculatively build each chunk at local column phase zero, then commit
   only when global prefix counts prove every chunk is row-aligned; otherwise use the
   generic second-pass fallback;
4. for documents, retain 8-byte source spans during that same validation pass;
5. expose chunks in source order without copying their payload buffers.

Document mode first proves that a loop itself reaches 8 MiB, preventing small metadata
loops near the start of a large file from repeatedly scanning the tail. This eligibility
probe recognizes whitespace, comments, quotes, semicolon text, bare values, and control
tokens directly without allocating token objects. Ambiguous or malformed input falls
back to the canonical lexer, so the probe never owns parse semantics. Against the
adjacent generic-lexer gate it reduced 6qnr/3j3q full-document latency by 23.2/12.1%.
Custom value/row limits that could become observable force the serial path. Workers
report located parser errors, and the coordinator selects the earliest source error
before the first loop control token.

### Input and release pipeline

`flate2` uses the portable pure-Rust `zlib-rs` backend. Gzip output gets a bounded reserve
limited by compressed size, decompressed limit, expansion-ratio limit, and 256 MiB. This
reduced large gzip end-to-end time by 17.6-22.6% before intra-file parallelism and lowered
one-shot 3j3q gzip RSS by 16.7%.

Raw and decompressed UTF-8 input is moved into the immutable `SourceBuffer` rather than
copied from an owned `Vec<u8>` into a second reference-counted string. The outer
`Arc<SourceInner>` already provides shared ownership. An adjacent copy-vs-consume A/B
improved 6qnr projection/full by 4.0/1.5% and 3j3q by 9.4/5.6%, while reducing large-text
peak RSS by 7-22%.

`tools/build_pgo.py` builds an instrumented wheel, trains it against verified text,
BinaryCIF, gzip, projection, and document workloads, merges profiles using the active
rustc's `llvm-profdata`, rebuilds, and publishes a wheel plus fingerprint metadata.
Profiles are deliberately not stored in Cargo configuration because they are coupled to
the compiler, source, target, Python ABI, and workload.

## Rejected or parked experiments

- Parallel full-document vectors followed by concatenation: rejected after measuring
  45-87 MB/s and 0.74-4.4 GB RSS; exact chunk ownership replaced it.
- `target-cpu=native`, thin LTO, and static native flags: no stable gain.
- one codegen unit/fat LTO: helped BinaryCIF by up to 6.5% but regressed text about 15%.
- mmap input: an owned-buffer source copy was removed without mapping or `unsafe`; mapping
  itself still has only the remaining file-read fraction as its ceiling and would add
  platform-specific ownership complexity.
- SIMD for ordinary PDB tokens: corpus tokens are too short; a hybrid is only justified
  if adversarial long unquoted-token throughput becomes a maintained guardrail.
- Reusing the serial 8 MiB eligibility prefix: 9-16% slower because it serialized cell
  materialization and reduced parallel work.
- Removing the eligibility probe: about 8x slower on large PDB files because many small
  metadata loops repeatedly launched worker waves against a large remaining tail.
- A specialized batch loop scanner: improved quote-free 3j3q slightly but regressed
  quote-heavy 6qnr full construction by 5.4%; rejected rather than corpus-dispatched.

## Validation

- the complete Rust unit/integration suite plus one doctest, including a >8 MiB
  forced-serial versus parallel adversarial loop with multiline text, quoted controls,
  predicates, row spans, and `stop_` handoff.
- `cargo clippy --release --all-targets -- -D warnings`.
- 65 Python tests against the final PGO wheel.
- strict Ruff and mypy checks for the PGO and stress tools.
- all five pinned PDB structures in raw CIF, BinaryCIF, and gzip: exact manifest sizes and
  SHA-256 hashes, expected atom counts, and equal projected Arrow tables.
- all five full raw CIF, gzip CIF, and BinaryCIF documents serialized canonically,
  reparsed, and reproduced identical canonical bytes; the largest output was 225 MB.

## Residual ceiling

The final PGO build reaches 1.60-2.09 GB/s text projection and 0.97-1.51 GB/s full
construction on the two large corpus files. Final average utilization ranges from 2.7
to 7.4 cores end to end; the parallel atom-loop work is diluted by serial metadata,
input, and output ownership. Pre-B15 plain-build profiling attributed 57% of projection
and 74% of document worker samples to the strict lexer, boundary indexing about 1%, and
cell stores at most 14%; B15 then removed the 24.5% main-thread eligibility-lexer share.
Modeled document traffic is roughly 4 GB/s versus 169 GB/s measured `memcpy`, so the
residual is branch/control-flow compute, not DRAM bandwidth. Warm file reads are
12-24 GB/s. Nibbler is therefore near the practical limit of this exact generic-lexer
architecture, not the machine's literal byte-copy limit. Another large gain requires a
universally faster strict lexer; a quote-sensitive specialization was measured and
rejected because it regressed a representative PDB file.
