"""Smoke tests for the private Rust extension."""

from __future__ import annotations

import nibbler


def test_native_metadata_matches_python_package() -> None:
    assert nibbler.contract_version() == 1
    assert nibbler.native_version() == nibbler.__version__
