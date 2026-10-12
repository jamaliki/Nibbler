# Structure preparation qualification

The local CPU measurements used an Apple M4 Max with 16 logical CPUs and 128 GiB
RAM, macOS 26.7, CPython 3.12.13, and Rust 1.97.1. Both Nibbler extensions used
`maturin build/develop --release` with the same lockfile and default package features.
The baseline is Nibbler `63b9e19bca04d4d53af2b5793af9b3ba043b36b1`;
the baseline kaveh-data selector is `78ca71629b353e5c47cdca6496cf7d93031cd6e0`.
The measurements cover the preparation changes on those baseline checkouts.
The Nibbler PR rebases the preparation diff onto `main`
`3a72ea6`; the saved timings describe the measured build, not a new run after that
rebase. The tests below ran again on the PR branches. The kaveh-data PR uses
`codex/pdb-geometry-audit` at `1684a99` as its base, excluding later local
geometry-analysis commits.

The measured native extension SHA-256 is
`25901d5815e5f872e421aabb6dc5ac9317bf5aecc68321020832198b87991fff`.
The Cargo lockfile SHA-256 is
`7cb2dc6cd75145225e7a94ac873e1af0f8fedf5482747f10794ad93fe4cda299`.
Each [worker report](../benchmarks/preparation-results/) records input hashes,
output hashes, all samples, p50/p90, and peak RSS. `summary.json` pools samples
from matching workers. Document and semantic-read output hashes match the baseline.

Saved reports refer to the former local fixture paths and CCD subset. The public
CIF blobs no longer belong in the source tree. The test downloader pins fresh RCSB
component bytes and coordinate content in `tests/preparation_corpus.json`.
Thirty-five coordinate files match the old content. RCSB's newer `1TWR` revision
keeps the observed atoms and authored links but corrects embedded `VER` hydrogen
templates. Tests qualify those downloads separately; the saved product timings did
not rerun with the new CCD subset or the remediated `1TWR` input.

## Latency

Workers loaded input bytes before timing and ran three warmups and 31 samples.
The first comparisons ran in baseline/current/current/baseline order. Final
selection workers reran after the hydrogen mask and crystal-image contact checks.
Read, assembly, and product measurements precede those contact-only changes;
those workloads do not call the conformer contact code. No other CPU test ran
during timing. Disk reads, registry loading, and result hashing/export lie outside
the measured operations.

Times below are pooled p50 / p90 in milliseconds.

| Operation and input | Baseline | Current |
| --- | ---: | ---: |
| Document read, 1CBN | 0.451 / 0.493 | 0.444 / 0.505 |
| Document read, 6QNR | 80.879 / 83.165 | 80.014 / 81.945 |
| Semantic read, 1CBN | 0.935 / 1.027 | 0.980 / 1.053 |
| Semantic read, 6DP5 | 3.273 / 3.369 | 3.368 / 3.523 |
| Semantic read, 6QNR | 239.347 / 258.351 | 236.891 / 253.684 |
| kaveh-data selection, 1CBN | 1.996 / 2.114 | 1.698 / 1.797 |
| kaveh-data selection, 6LZB | 11.104 / 11.540 | 7.374 / 7.521 |
| kaveh-data selection, 1RCQ | 11.401 / 11.869 | 7.640 / 7.912 |
| kaveh-data selection, 1IVO | 21.373 / 21.977 | 8.118 / 8.297 |
| kaveh-data selection, 6QNR | 1126.548 / 1145.984 | 912.363 / 953.589 |
| Assembly expansion, 7MUB | 34.069 / 34.435 | 34.258 / 34.751 |
| Assembly expansion, 1XAN | 7.549 / 7.674 | 7.891 / 8.037 |
| Assembly expansion, 6QNR | 1030.250 / 1043.956 | 1027.837 / 1031.630 |

Replacing the Gemmi/category-copy/reparse selector reduces measured selection
latency by 15% on 1CBN, 34% on 6LZB, 33% on 1RCQ, 62% on 1IVO, and 19% on 6QNR.
This includes native decoding, clash selection, and the filtered document view.
Whole-ligand overlap adds 0.19 ms on 7MUB and 0.34 ms on 1XAN. Ordinary semantic
reads increase by 0.046 ms on 1CBN and 0.095 ms on 6DP5; the large read stays near
baseline. These are local observations, not portable speed claims.

Native product preparation succeeds on 29 guide cases, with p50 from 0.129 to
35.075 ms and a median per-case p50 of 2.703 ms. The largest successful graph,
4NDZ, contains 69,355 expected atoms and takes 35.075 / 41.141 ms. Expected
rejections cost 145.079 / 166.192 ms for 1VQ7 and 2.714 / 6.542 ms for 7P1D.
The rejected 1VQ7 case still pays for graph construction and resolution before
validation fails. JSON graph export is outside these timings.

## Memory

Whole-worker peak RSS includes setup, retained input documents, allocator caches,
and untimed output hashing/export. It is not a unique decoded-shard footprint.

| Worker workload | Baseline peak range (MiB) | Current peak range (MiB) |
| --- | ---: | ---: |
| Document reads | 405.5-407.6 | 434.3-436.4 |
| Semantic reads | 917.0-947.9 | 980.7-1022.8 |
| kaveh-data selection | 1816.9-1933.3 | 1483.5-1746.1 |
| Assembly expansion | 532.0-533.1 | 534.0-535.5 |

The ordinary-read workers show higher peak RSS despite similar latency. This gate
does not prove memory neutrality for those workers. Selection uses less peak RSS
in this run; assembly differs by under 1%. The product worker peaks at 641.3 MiB,
including Python JSON export. No change bounds a resident training cache or adds
training-time disk reads.

## Correctness and limits

- Rust tests pass without Reduce3 (118) and with Reduce3 (133), including 18 new
  product, conformer, and overlap tests.
- The 57 preparation regressions and 12 download/cache tests pass. They cover
  observed coordinates, charges,
  hydrogen counts, pi additions, resonance, missing leaving atoms, alternate link
  states, deliberate failures, and complete special-position ligands.
  A fresh-cache run downloads the corpus during pytest and passes all 69 tests
  in 7.42 s. An offline rerun passes in 2.98 s. These whole-test times include
  setup and assertions; they are separate from the operation benchmarks above.
- kaveh-data's PR suite passes 256 tests, with three optional corpus checks
  skipped. The measured checkout included 13 additional geometry-analysis tests
  and passed all 269. The optional translated two-copy 1T4W assembly/Reduce3 test
  passed on that checkout.
- The full Nibbler Python run has 209 passes and 47 optional skips. Its one failure,
  `test_acid_hydrogens_prefer_syn_unless_turned_off`, also fails on the baseline
  with the selected chemistry tree (`HZ2` is absent). The optional parallel real
  kaveh-data preparation gate cannot pass with that partial tree because HEM, PO4,
  TPO, SEP, and MN restraints are absent.
- Modified Python APIs and benchmark/generator files pass strict mypy and Ruff.
  Library Clippy passes; the Reduce3-enabled run needs the pre-existing
  `too_many_arguments` warning suppressed. Full-repo mypy and all-target Clippy
  still report pre-existing issues outside these changes. Rust format checks also
  fail in unchanged `src/validation/mod.rs` and `src/validation/report.rs`.

The conda Arrow 25.0.0 mimalloc pool and Nibbler's existing mimalloc allocator
corrupt each other's buffers on this host. The baseline reproduces that conflict.
Using Arrow's system pool fixes a small reproducer but does not qualify the full
conda test run. Full Python qualification uses the pip Arrow 25.0.1 build and
`ARROW_DEFAULT_MEMORY_POOL=system`. No allocator code changed in this work.

These gates do not cover a full PDB archive rebuild, feature export from the product
graph, a deployed pipeline, or GPU training. Nibbler's product graph remains
separate from kaveh-data's CCD-derived ligand atom features.

## Reproduction

```console
micromamba run -p .mamba/nibbler-dev maturin develop --release
micromamba run -p .mamba/nibbler-dev python -m tools.fetch_preparation_data
micromamba run -p .mamba/nibbler-dev python benchmarks/preparation_cost.py \
  --workload product --ccd .cache/structure-preparation-v1/components.cif \
  --file .cache/structure-preparation-v1/4ndz.cif.gz --warmups 3 --samples 31
```

For baseline comparisons, build a wheel from the pinned Git archive, extract it
outside the checkout, and select it through `PYTHONPATH`. Run the current benchmark
script against both extensions. `pipeline-select` also needs the matching
kaveh-data source tree, Gemmi, NumPy, SciPy, and the qualified Arrow build.
