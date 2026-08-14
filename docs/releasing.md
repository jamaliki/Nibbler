# Releasing Nibbler

Portable release artifacts are built once, installed and qualified on their target
platform, then assembled before any publication job receives credentials. A manual
workflow run produces a release candidate without publishing it. Only an exact version
tag can publish.

The release matrix contains 25 interpreter-specific wheels plus one source distribution:

| Platform | Architecture | CPython |
| --- | --- | --- |
| manylinux 2.17 | x86_64, AArch64 | 3.10-3.14 |
| macOS 11+ | arm64, x86_64 | 3.10-3.14 |
| Windows | x86_64 | 3.10-3.14 |

The Rust minimum supported version is 1.88. Release wheels are portable builds; the
optional host-specific PGO artifact is not uploaded.

## One-time PyPI configuration

Before pushing the first release tag, create a PyPI
[pending trusted publisher](https://docs.pypi.org/trusted-publishers/creating-a-project-through-oidc/)
for:

- owner: `jamaliki`
- repository: `Nibbler`
- workflow: `release.yml`
- environment: `pypi`

Create the matching GitHub
[`pypi` environment](https://docs.github.com/en/actions/how-tos/deploy/configure-and-manage-deployments/manage-environments).
If the repository plan supports environment reviewers, require one. Otherwise, tag
creation and trusted-publisher configuration are the manual approval boundary. No
long-lived PyPI token is used. The publish job alone receives `id-token: write`; build
and qualification jobs have read-only repository access.

## Release candidate

1. Update `CHANGELOG.md` and synchronize the versions in `pyproject.toml` and
   `Cargo.toml`. Keep the changelog entry marked `Unreleased` while qualifying the
   candidate.
2. Run `make check` and `make release-artifacts` locally. `dist/` must be absent or empty.
3. Open and merge the release-candidate change into `main` after continuous integration
   passes.
4. Run the `Release` workflow manually from that exact commit. Download
   `release-packages` and `release-metadata`, verify `SHA256SUMS`, confirm 25 wheels and
   one source distribution, and install representative wheels.
5. Review the scheduled `PDB regression` report for correctness, throughput, and peak
   RSS.

Manual runs never upload to PyPI or create a GitHub release.

## Publish

Replace `Unreleased` in the version entry with the release date, rerun
`make release-metadata`, merge that change into `main`, and create an annotated tag on
that exact commit. The tag must match the package version:

```console
git tag -a v0.1.0 -m "Nibbler 0.1.0"
git push origin v0.1.0
```

The tag-triggered workflow rebuilds and requalifies every artifact, publishes the
complete set to PyPI through the configured environment, and creates the GitHub release
with the same archives and checksums. A metadata, build, qualification, or publication
failure leaves no partial GitHub release.

After publication, install `nibbler-cif==<version>` into a fresh environment, verify the
native and Python versions, run a minimal text and BinaryCIF read, and compare the GitHub
release checksums with `SHA256SUMS`.

Profile-guided wheels remain a separately fingerprinted native workflow. They are not mixed
into the portable release matrix until dedicated runners can train and qualify them on every
published target.
