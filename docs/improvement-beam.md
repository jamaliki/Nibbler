# Historical performance experiment ledger

> This file is an engineering history, not the current design. It is the single place
> where Nibbler records superseded baselines, promoted experiments, failed approaches,
> and restart conditions. See [performance-report.md](performance-report.md) for the
> implementation and qualification that are current.

## Hardware-limit cycle: final PGO checkpoint, 2026-08-12

This cycle started from a fresh exact-tree plain wheel at 930/1,020 MB/s text projection
and 499/661 MB/s full-document construction for 6qnr/3j3q. Paired measurements on the
promoted plain tree reach 1,373/1,772 MB/s projection and 856/1,330 MB/s full construction.
The fingerprinted PGO wheel reaches 1,603/2,092 MB/s projection and 969/1,512 MB/s
full construction. Against a same-source, same-compiler plain wheel, PGO adds 17-18%
large text projection and about 13-14% large text full-document throughput.

| Rank | Beam | Type | Status | Hypothesis | Decisive experiment |
|---:|---|---|---|---|---|
| 1 | B5 single-scan projection | structural/high-risk | promoted | speculative chunk outputs are committed only after global prefixes prove row alignment | +39-41% versus the cycle plain baseline |
| 2 | B6 parallel boundary index | exploit | promoted | exact two-state semicolon summaries compose across byte blocks | index -94-96%; exhaustive equality |
| 3 | B7 hardware profile | diagnostic/near-miss | complete | current loss is compute, memory, allocation, or synchronization | lexer 57-74%; boundary index about 1% |
| 4 | B8 direct retained documents | structural | promoted | the validation pass can own compact cell chunks without a rescan or concatenation | full document +20-53% before B10 |
| 5 | B10 8-byte source cells | representation | promoted | token spans already encode content and quote style | full -2-6% latency; RSS -19-23% |
| 6 | B15 lightweight eligibility probe | serial-path | promoted | allocation-free classification can replace a generic 8 MiB lexer pass without owning semantics | full -23.2/-12.1% adjacent A/B |
| 7 | B16 owned UTF-8 source | ownership | promoted | the loader can move its validated buffer into the existing shared source | projection -4.0/-9.4%; full -1.5/-5.6% adjacent A/B |

Success requires a paired end-to-end gain outside run noise, not only a faster isolated
counter. The final combination must retain the one grammar implementation, exact errors,
resource limits, small-file latency, bounded RSS, and all corpus equality gates.

## Prior cycle checkpoint: 2026-08-12

The search started from Phase 6 text projection at 321/278 MB/s and text full-document
construction at 272/243 MB/s for 6qnr/3j3q. Its PGO checkpoint reached 1,174/1,258 MB/s
projection and 588/787 MB/s full construction before the hardware-limit cycle above.

| Rank | Beam | Status | Evidence | Decision |
|---:|---|---|---|---|
| 1 | B3 deterministic intra-file parallelism | promoted | 3.65-4.53x text projection; 2.16-3.24x text full versus initial baseline | bounded production path |
| 2 | B1 compact document storage | promoted | BinaryCIF full 3.47-4.18x; RSS down 55-77% | shared spans and retained columns |
| 3 | B4 build/input pipeline | promoted | PGO +13-17% on final parallel text; gzip backend/preallocation -17.6-22.6% before parallelism | native artifact workflow |
| 4 | B2 fused strict scanner | promoted | +8.7-10.2% large text, small file neutral | one-pass scalar unquoted scan |

## Combination result

B1+B2+B3+B4 passed the complete Rust/Python suite and all five raw CIF, BinaryCIF, and
gzip corpus cases. PGO improved rather than disrupted the parallel path. Chunked Arrow
output avoided a final buffer concatenation; hoisted missing-state metadata removed an
accidental O(chunks x rows) export pass.

## Failed and parked branches

| Candidate | Status | Evidence | Restart condition |
|---|---|---|---|
| concatenate parallel document cell vectors | killed | 45-87 MB/s, up to 4.4 GB RSS | none; replaced by owned chunks |
| long-token SIMD | parked | 99.5% of real tokens shorter than 8 bytes | representative corpus changes |
| static LTO/CGU/native flags | parked | mixed/noisy; text regression up to 15% | new compiler or target |
| mmap source | parked | owned-buffer copy removed; only file-read fraction remains | a safe portable mapping proves a further end-to-end gain |
| retained serial 8 MiB prefix | killed | moved work onto the critical path; full latency +9-16% | none |
| remove document eligibility probe | killed | metadata loops repeatedly spawned workers; large files about 8x slower | none |
| rolling projection column counter | killed | large projection about 10% slower than optimized remainder | none |
| specialized batch loop scanner | killed | 3j3q up to 4% faster but 6qnr full 5.4% slower | a universal exact scanner with cross-corpus gain |
| stable-ABI Python wheels | parked | repeated large-PDB runs were about 2-3% slower than interpreter-specific wheels | PyO3 or compiler changes eliminate the measured loss |

## Remaining beams

1. A future lexer redesign must beat both quote-free 3j3q and quote-heavy 6qnr; the
   profile gives a realistic 15-35% opportunity, not another multi-x gain.
2. Revisit file mapping only with a safe ownership design and an end-to-end win; the
   redundant in-memory source copy is already gone.

Promotion continues to require corpus equality, stable errors and resource limits,
small-file latency, bounded RSS, clippy/tests, and an end-to-end gain outside run noise.
