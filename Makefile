DEV_PREFIX := $(CURDIR)/.mamba/nibbler-dev
RUN := micromamba run -p $(DEV_PREFIX)

.PHONY: bootstrap develop develop-release format lint typecheck test docs check corpus benchmarks pdb-corpus pdb-stress schemas compile-schemas validate-fixtures clean

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
	$(RUN) ruff format --check .
	$(RUN) ruff check .

typecheck:
	$(RUN) mypy --strict python benchmarks tools

test: develop
	$(RUN) cargo test --no-default-features
	$(RUN) pytest

docs:
	$(RUN) env RUSTDOCFLAGS="-D warnings" cargo doc --no-deps --no-default-features

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

check: lint typecheck test docs corpus validate-fixtures
	$(MAKE) benchmarks

clean:
	$(RUN) cargo clean
	rm -rf .mypy_cache .pytest_cache .ruff_cache
