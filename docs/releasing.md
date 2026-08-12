# Releasing Nibbler

Release artifacts are built once, installed and qualified on their target platform, then
assembled before any publication job receives credentials. A manual workflow run produces
a release-candidate artifact without publishing it. Only an exact version tag can publish.

## One-time PyPI configuration

Create a PyPI trusted publisher for:

- owner: `jamaliki`
- repository: `Nibbler`
- workflow: `release.yml`
- environment: `pypi`

Create the matching GitHub `pypi` environment and require manual approval. No long-lived
PyPI token is used. The publish job alone receives `id-token: write`; build and qualification
jobs have read-only repository access.

## Release candidate

1. Update `CHANGELOG.md` and synchronize the versions in `pyproject.toml` and
   `Cargo.toml`. Keep the changelog entry marked `Unreleased` while qualifying the
   candidate.
2. Run `make check` and `make release-artifacts` locally. `dist/` must be absent or empty.
3. Push the candidate commit to `main`.
4. Run the `Release` workflow manually. Download `release-packages` and
   `release-metadata`, verify `SHA256SUMS`, and install representative wheels.
5. Review the scheduled `PDB regression` report for correctness, throughput, and peak RSS.

Manual runs never upload to PyPI or create a GitHub release.

## Publish

Replace `Unreleased` in the version entry with the release date, rerun
`make release-metadata`, commit, and push. Then create and push an annotated tag matching
the package version exactly:

```console
git tag -a v0.1.0 -m "Nibbler 0.1.0"
git push origin v0.1.0
```

The tag-triggered workflow rebuilds and requalifies every artifact, publishes the complete
set to PyPI through the protected environment, and creates the GitHub release with the same
archives and checksums. A metadata, build, qualification, or publication failure leaves no
partial GitHub release.

Profile-guided wheels remain a separately fingerprinted native workflow. They are not mixed
into the portable release matrix until dedicated runners can train and qualify them on every
published target.
