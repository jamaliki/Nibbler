# Benchmark contract

Phase 0 froze workloads and comparison engines before optimization. Phase 2 adds an
executable, correctness-qualified `projected-atom-site` adapter for Nibbler and every
baseline. Each engine runs in a fresh process so peak RSS and imports are isolated.
The runner refuses to publish results from an unoptimized Nibbler extension.

```console
make develop-release
python -m benchmarks.run --list
python -m benchmarks.run --file tests/fixtures/chemistry/ligand_ion_water.cif --json
python -m benchmarks.run --synthetic-rows 100000 --warmups 1 --samples 3
```

`make benchmarks` installs the release extension automatically.

The report records the input digest, engine version, command, hardware, warmup/sample
counts, p50/p95/p99 latency, decompressed throughput, peak resident memory, and a digest
of projected rows. A semantic mismatch fails the command rather than publishing timing.
End-to-end workloads, rather than isolated token microbenchmarks, decide whether an
optimization is promoted.

## Real PDB stress corpus

Phase 6 adds five hash-pinned RCSB structures spanning a small protein, ligand/metal
chemistry, a multi-model NMR ensemble, a ribosome, and a 2.44-million-atom assembly.
The files are downloaded into the ignored `.cache/pdb-stress` directory; tests never
depend on network access.

```console
python -m tools.fetch_pdb_corpus
python -m benchmarks.pdb_stress --warmups 1 --samples 3
python -m benchmarks.pdb_stress --formats all --warmups 1 --samples 3
python -m benchmarks.pdb_stress --full-document --warmups 0 --samples 1
```

The staged runner measures warm file reads, native `_atom_site` projection, Arrow import,
and optional isolated full-document materialization separately. `--formats all` adds the
pinned gzip inputs and reports their throughput against logical decompressed CIF bytes.
CIF, BinaryCIF, and gzip typed projections must compare equal before timings are
reported. Manifest sizes, SHA-256 digests, and atom counts are enforced rather than only
printed. Workload rationale, hardware, software, and sample counts are included in the
manifest or JSON report.

On the 2026-08-12 16-core Apple M4 Max development machine, the final corpus-trained PGO
wheel projected the two large text CIF files at 1.60-2.09 GB/s, BinaryCIF at
1.97-2.07 GB/s, and gzip at 0.69-0.81 GB/s of logical CIF. Arrow import took 0.7 ms for
6qnr and 4.9 ms for 3j3q. Full text materialization reached 0.97-1.51 GB/s. These are
workload-specific observations, not portable thresholds.
