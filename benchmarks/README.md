# Benchmark contract

Phase 0 freezes workloads and comparison engines before optimization. It does not
publish meaningless timings for an unimplemented parser.

```console
python -m benchmarks.run --list
```

Timed adapters will be added with the shared Rust parser. Each result must record the
input digest, engine and version, command, hardware, warmup count, sample count,
distribution statistics, peak resident memory, and correctness result. End-to-end
workloads, rather than isolated token microbenchmarks, decide whether an optimization is
promoted.
