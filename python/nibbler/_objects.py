"""Immutable Python objects and optional dataframe interchange."""

from __future__ import annotations

from importlib import import_module
from typing import Literal, Protocol, TypeAlias, cast

from . import _core
from ._native import raise_write_error

MissingPolicy: TypeAlias = Literal["collapse", "columns", "extension"]


class _ArrowTable(Protocol):
    def to_pandas(self) -> object: ...


class CifDocument:
    """An immutable, order-preserving generic CIF document."""

    __slots__ = ("_native", "_schema")

    def __init__(self, native: _core._CifDocument, schema: str | None = None) -> None:
        self._native = native
        self._schema = schema

    @property
    def block_count(self) -> int:
        """Return the number of data and global blocks."""
        return self._native.block_count

    @property
    def schema(self) -> str | None:
        """Return the dictionary selector requested at read time, if any."""
        return self._schema

    def to_canonical(self) -> str:
        """Serialize the document using Nibbler's canonical CIF policy."""
        return self._native.to_canonical()

    def to_preserving(self) -> str:
        """Serialize with source ordering and valid original value lexemes."""
        return self._native.to_preserving()

    def to_binary(self) -> bytes:
        """Serialize the document as BinaryCIF 0.3 MessagePack."""
        return self._native.to_binary()

    def __repr__(self) -> str:
        return f"CifDocument(block_count={self.block_count})"


class CifTable:
    """A projected CIF category with lossless internal missing states."""

    __slots__ = ("_native",)

    def __init__(self, native: _core._CifTable) -> None:
        self._native = native

    @property
    def category(self) -> str:
        """Return the category name without a leading underscore."""
        return self._native.category

    @property
    def columns(self) -> tuple[str, ...]:
        """Return projected item names in output order."""
        return tuple(self._native.columns)

    def __len__(self) -> int:
        return len(self._native)

    def __arrow_c_stream__(self, requested_schema: object = None) -> object:
        """Export this table through the zero-copy Arrow C Stream protocol."""
        return self._native.__arrow_c_stream__(requested_schema)

    def to_pyarrow(self, *, missing: MissingPolicy = "collapse") -> object:
        """Return a PyArrow table without making PyArrow a runtime dependency."""
        try:
            pa = import_module("pyarrow")
        except ImportError as error:
            raise ImportError(
                "CifTable.to_pyarrow() requires the optional 'pyarrow' package"
            ) from error
        stream = self._native._with_missing(missing)
        return pa.RecordBatchReader.from_stream(stream).read_all()

    def to_polars(self, *, missing: MissingPolicy = "collapse") -> object:
        """Return a Polars dataframe without making Polars a runtime dependency."""
        try:
            pl = import_module("polars")
        except ImportError as error:
            raise ImportError(
                "CifTable.to_polars() requires the optional 'polars' package"
            ) from error
        return pl.DataFrame(self._native._with_missing(missing))

    def to_pandas(self, *, missing: MissingPolicy = "collapse") -> object:
        """Return a pandas dataframe through the PyArrow interchange path."""
        try:
            import_module("pandas")
        except ImportError as error:
            raise ImportError(
                "CifTable.to_pandas() requires the optional 'pandas' package"
            ) from error
        arrow_table = cast(_ArrowTable, self.to_pyarrow(missing=missing))
        return arrow_table.to_pandas()

    def __repr__(self) -> str:
        return (
            f"CifTable(category={self.category!r}, rows={len(self)}, "
            f"columns={len(self.columns)})"
        )


class MmcifModel:
    """An immutable, source-backed macromolecular coordinate model."""

    __slots__ = ("_native",)

    def __init__(self, native: _core._MmcifModel) -> None:
        self._native = native

    @property
    def profile(self) -> Literal["pdbx", "modelcif"]:
        """Return the semantic profile used to construct this model."""
        return self._native.profile

    @property
    def entry_id(self) -> str:
        """Return the PDBx entry identifier."""
        return self._native.entry_id

    @property
    def entity_count(self) -> int:
        """Return the number of semantic entities."""
        return self._native.entity_count

    @property
    def asym_unit_count(self) -> int:
        """Return the number of entity instances."""
        return self._native.asym_unit_count

    @property
    def component_count(self) -> int:
        """Return the number of resolved or explicitly unresolved components."""
        return self._native.component_count

    @property
    def atom_site_count(self) -> int:
        """Return the number of coordinate atom sites."""
        return self._native.atom_site_count

    @property
    def connection_count(self) -> int:
        """Return the number of explicit inter-site connections."""
        return self._native.connection_count

    @property
    def entity_kinds(self) -> tuple[str, ...]:
        """Return entity kinds in source order."""
        return tuple(self._native.entity_kinds)

    @property
    def component_ids(self) -> tuple[str, ...]:
        """Return component identifiers in stable order."""
        return tuple(self._native.component_ids)

    @property
    def prediction_model_count(self) -> int:
        """Return the number of deposited prediction models."""
        return self._native.prediction_model_count

    @property
    def target_entity_count(self) -> int:
        """Return the number of declared prediction targets."""
        return self._native.target_entity_count

    @property
    def template_count(self) -> int:
        """Return the number of structural templates."""
        return self._native.template_count

    @property
    def qa_metric_count(self) -> int:
        """Return the number of quality-metric definitions."""
        return self._native.qa_metric_count

    @property
    def qa_value_count(self) -> int:
        """Return the number of global, local, and pairwise QA values."""
        return self._native.qa_value_count

    @property
    def software_names(self) -> tuple[str, ...]:
        """Return prediction software names in source order."""
        return tuple(self._native.software_names)

    @property
    def qa_metric_names(self) -> tuple[str, ...]:
        """Return quality-metric names in source order."""
        return tuple(self._native.qa_metric_names)

    @property
    def qa_metric_modes(self) -> tuple[str, ...]:
        """Return quality-metric modes in source order."""
        return tuple(self._native.qa_metric_modes)

    def to_document(self, *, mirror_local_qa_metric: int | None = None) -> CifDocument:
        """Return a generic document in canonical profile order.

        ``mirror_local_qa_metric`` explicitly copies one local ModelCIF QA metric
        into output B factors for viewer compatibility. The QA records remain intact.
        """
        if self.profile == "pdbx" and mirror_local_qa_metric is not None:
            raise ValueError("mirror_local_qa_metric requires a ModelCIF model")
        try:
            document = self._native.to_document(mirror_local_qa_metric)
        except ValueError as error:
            raise_write_error(error)
        return CifDocument(document, self.profile)

    def __repr__(self) -> str:
        return (
            f"MmcifModel(profile={self.profile!r}, entry_id={self.entry_id!r}, "
            f"entities={self.entity_count}, atoms={self.atom_site_count})"
        )
