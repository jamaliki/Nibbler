# Contributing to Nibbler

Nibbler welcomes focused bug reports and changes that improve correctness, robustness,
performance, or clarity. Please use GitHub Discussions for open-ended design proposals
and GitHub Issues for reproducible defects or bounded feature requests.

## Development setup

The supported development environment is repository-local:

```console
make bootstrap
make check
```

`make bootstrap` creates `.mamba/nibbler-dev` from `environment-dev.yml`. The complete
gate formats and lints Rust and Python, type-checks Python, builds the extension, runs
unit, integration, robustness, documentation, corpus, schema, fixture, release-metadata,
and benchmark checks.

For a minimal installation test from a checkout:

```console
uv venv --python 3.12
uv pip install .
```

This requires Rust 1.88 or newer.

## Changes

- Keep changes narrowly scoped and preserve strict CIF behavior.
- Add tests at the closest stable public boundary.
- Do not make parser recovery implicit or add runtime network access.
- Treat malformed input, resource limits, and deterministic error ordering as part of
  correctness.
- Include before/after measurements for performance claims. Use release builds, pinned
  fixtures, exact-output checks, and the process described in
  [`benchmarks/README.md`](benchmarks/README.md).
- Update user-facing documentation and `CHANGELOG.md` with behavioral changes.

Run `make check` before opening a pull request. CI uses Rust 1.88, so code accepted by a
newer local compiler must still pass that minimum-version gate.

## Reporting security issues

Do not open a public issue for a suspected vulnerability. Follow
[`SECURITY.md`](SECURITY.md) instead.

By participating, you agree to follow the [`CODE_OF_CONDUCT.md`](CODE_OF_CONDUCT.md).
