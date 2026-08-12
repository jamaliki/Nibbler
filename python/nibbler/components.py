"""Explicit, immutable chemical-component registries."""

from __future__ import annotations

from . import _core, cif
from ._input import Source
from ._objects import CifDocument
from .errors import ChemistryError


class Registry:
    """Chemical component definitions loaded from a caller-selected local CIF cache."""

    __slots__ = ("_native",)

    def __init__(self, native: _core._ComponentRegistry) -> None:
        self._native = native

    @classmethod
    def from_ccd_cache(cls, source: Source | CifDocument) -> Registry:
        """Load an existing immutable CCD CIF artifact without network access."""
        document = source if isinstance(source, CifDocument) else cif.read(source)
        if not isinstance(document, CifDocument):  # pragma: no cover - defensive
            raise TypeError("a CCD cache must parse as a full CIF document")
        try:
            return cls(_core.build_component_registry(document._native))
        except ValueError as error:
            if len(error.args) != 3:
                raise
            code, message, context = error.args
            raise ChemistryError(
                code=str(code),
                message=str(message),
                context=tuple(str(value) for value in context),
            ) from None

    def __len__(self) -> int:
        return len(self._native)

    def __repr__(self) -> str:
        return f"Registry(components={len(self)})"
