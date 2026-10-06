DEV_PREFIX := $(CURDIR)/.mamba/nibbler-dev
RUN := micromamba run -p $(DEV_PREFIX)

.PHONY: bootstrap develop develop-release format lint typecheck test robustness docs check corpus benchmarks pdb-corpus pdb-stress schemas compile-schemas validate-fixtures release-metadata release-artifacts clean

RELEASE_DIR ?= dist

bootstrap:
	micromamba create -y -p $(DEV_PREFIX) -f environment-dev.yml

develop:
	$(RUN) maturin develop

develop-release:
	$(RUN) maturin develop --release

format:
	$(RUN) cargo fmt --all
	$(RUN) ruff format .

lint:
	$(RUN) cargo fmt --all --check
	$(RUN) cargo clippy --all-targets --no-default-features -- -D warnings
	$(RUN) cargo clippy --all-targets --features robustness -- -D warnings
	$(RUN) cargo clippy --all-targets --features python,reduce3 -- -D warnings
	$(RUN) ruff format --check .
	$(RUN) ruff check .

typecheck:
	$(RUN) mypy --strict python benchmarks tools

test: develop
	$(RUN) cargo test --no-default-features
	$(RUN) cargo test --features reduce3
	$(RUN) pytest

robustness:
	$(RUN) cargo test --release --features robustness --test robustness

docs:
	$(RUN) env RUSTDOCFLAGS="-D warnings" cargo doc --no-deps --no-default-features
	$(RUN) env RUSTDOCFLAGS="-D warnings" cargo doc --no-deps --features reduce3

corpus:
	$(RUN) python -m tools.verify_corpus

benchmarks: develop-release
	$(RUN) python -m benchmarks.run --file tests/fixtures/chemistry/ligand_ion_water.cif --warmups 1 --samples 1

pdb-corpus:
	$(RUN) python -m tools.fetch_pdb_corpus

pdb-stress: develop-release
	$(RUN) python -m benchmarks.pdb_stress --warmups 1 --samples 3

schemas:
	$(RUN) python -m tools.fetch_schemas

compile-schemas: schemas
	$(RUN) python -m tools.compile_schemas

validate-fixtures: develop schemas
	$(RUN) python -m tools.validate_fixtures

release-metadata:
	$(RUN) python -m tools.validate_release --metadata-only

release-artifacts: release-metadata
	@test ! -e "$(RELEASE_DIR)" || test -z "$$(find "$(RELEASE_DIR)" -mindepth 1 -maxdepth 1 -print -quit)" || (echo "$(RELEASE_DIR) must be empty"; exit 1)
	mkdir -p "$(RELEASE_DIR)"
	$(RUN) maturin build --release --locked --compatibility pypi --sdist --features extension-module --out "$(RELEASE_DIR)"
	$(RUN) twine check --strict "$(RELEASE_DIR)"/*
	$(RUN) python -m tools.validate_release --dist "$(RELEASE_DIR)" --expected-wheel-count 1
	$(RUN) python -m tools.qualify_release "$(RELEASE_DIR)"

check: lint typecheck test robustness docs corpus validate-fixtures release-metadata
	$(MAKE) benchmarks

clean:
	$(RUN) cargo clean
	rm -rf .mypy_cache .pytest_cache .ruff_cache
