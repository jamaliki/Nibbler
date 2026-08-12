"""Tests for the conventional and themed Python API facades."""

from __future__ import annotations

import pytest

import nibbler
from nibbler import cif


def test_parse_facade_is_the_conventional_function() -> None:
    assert nibbler.chomp is cif.read
    assert nibbler.feast is cif.scan


def test_profile_dump_requires_a_semantic_model() -> None:
    with pytest.raises(TypeError, match="requires MmcifModel"):
        nibbler.dump(object(), "out.cif", profile="modelcif")


def test_profile_and_schema_are_mutually_exclusive() -> None:
    with pytest.raises(TypeError, match="mutually exclusive"):
        nibbler.sniff(object(), profile="pdbx", schema="pdbx")


def test_profile_validation_rejects_unknown_profiles() -> None:
    with pytest.raises(ValueError):
        nibbler.sniff(object(), profile="deposition")
