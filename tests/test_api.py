"""Tests for the conventional and themed Phase 0 APIs."""

from __future__ import annotations

from collections.abc import Callable
from typing import NoReturn

import pytest

import nibbler
from nibbler import FeatureUnavailableError, Profile, cif


def test_parse_facade_is_the_conventional_function() -> None:
    assert nibbler.chomp is cif.read
    assert nibbler.feast is cif.scan


@pytest.mark.parametrize(
    ("call", "operation", "phase"),
    [
        (lambda: nibbler.chomp(b"data_example\n"), "nibbler.cif.read", 2),
        (lambda: nibbler.feast([]), "nibbler.cif.scan", 2),
        (lambda: nibbler.sniff(object()), "nibbler.cif.validate", 3),
        (
            lambda: nibbler.sniff(object(), profile=Profile.PDBX),
            "nibbler.mmcif.validate[pdbx]",
            4,
        ),
        (
            lambda: nibbler.spit(object(), "out.cif"),
            "nibbler.cif.write",
            3,
        ),
        (
            lambda: nibbler.spit(object(), "out.cif", profile="modelcif"),
            "nibbler.mmcif.write[modelcif]",
            4,
        ),
    ],
)
def test_unimplemented_operations_fail_explicitly(
    call: Callable[[], NoReturn], operation: str, phase: int
) -> None:
    with pytest.raises(FeatureUnavailableError) as raised:
        call()
    assert raised.value.operation == operation
    assert raised.value.required_phase == phase


def test_profile_and_schema_are_mutually_exclusive() -> None:
    with pytest.raises(TypeError, match="mutually exclusive"):
        nibbler.sniff(object(), profile="pdbx", schema="pdbx")


def test_profile_validation_rejects_unknown_profiles() -> None:
    with pytest.raises(ValueError):
        nibbler.sniff(object(), profile="deposition")
