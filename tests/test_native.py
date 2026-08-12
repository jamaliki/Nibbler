"""Smoke tests for the private Rust extension."""

from __future__ import annotations

import nibbler


def test_native_metadata_matches_python_package() -> None:
    assert nibbler.contract_version() == 1
    assert nibbler.native_version() == nibbler.__version__


def test_public_objects_are_the_native_classes() -> None:
    assert nibbler.CifDocument is nibbler._core.CifDocument
    assert nibbler.CifTable is nibbler._core.CifTable
    assert nibbler.MmcifModel is nibbler._core.MmcifModel
    assert nibbler.components.Registry is nibbler._core.Registry
    assert nibbler.CifDocument.__module__ == "nibbler"
    assert nibbler.components.Registry.__module__ == "nibbler.components"
    assert not hasattr(nibbler.MmcifModel, "to_document")
