DEV_PREFIX := $(CURDIR)/.mamba/nibbler-dev
RUN := micromamba run -p $(DEV_PREFIX)

.PHONY: bootstrap develop format lint typecheck test docs check corpus benchmarks schemas validate-fixtures clean

bootstrap:
	micromamba create -y -p $(DEV_PREFIX) -f environment-dev.yml

develop:
	$(RUN) maturin develop

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

benchmarks:
	$(RUN) python -m benchmarks.run --list

schemas:
	$(RUN) python -m tools.fetch_schemas

validate-fixtures: schemas
	$(RUN) python -m tools.validate_fixtures

check: lint typecheck test docs corpus benchmarks

clean:
	$(RUN) cargo clean
	rm -rf .mypy_cache .pytest_cache .ruff_cache
