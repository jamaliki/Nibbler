"""Explicit, immutable chemical-component registries."""

from __future__ import annotations

from . import _core, cif
from ._core import CifDocument
from ._core import Registry as Registry
from ._input import Source
from ._native import raise_chemistry_error


def read(source: Source | CifDocument) -> Registry:
    """Load an immutable registry from a caller-selected local CCD artifact."""
    document = source if isinstance(source, CifDocument) else cif.read(source)
    try:
        return _core.build_component_registry(document)
    except ValueError as error:
        raise_chemistry_error(error)
