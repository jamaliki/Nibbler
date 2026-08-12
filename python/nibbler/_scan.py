"""Bounded deterministic orchestration over the native worker pool."""

from __future__ import annotations

from collections.abc import Iterable, Iterator, Sequence
from importlib import import_module
from typing import Literal, TypeAlias, overload

from . import _core
from ._input import (
    PredicateSpec,
    Source,
    _batch_error,
    _normalize_source,
    _raise_native_error,
)
from ._objects import CifDocument, CifTable
from .errors import BatchError

ReadResult: TypeAlias = CifDocument | CifTable


class BatchDiagnostics(Sequence[BatchError]):
    """Collected batch failures in deterministic source order."""

    __slots__ = ("_errors",)

    def __init__(self, errors: Sequence[BatchError]) -> None:
        self._errors = tuple(errors)

    @overload
    def __getitem__(self, index: int) -> BatchError: ...

    @overload
    def __getitem__(self, index: slice) -> tuple[BatchError, ...]: ...

    def __getitem__(self, index: int | slice) -> BatchError | tuple[BatchError, ...]:
        return self._errors[index]

    def __len__(self) -> int:
        return len(self._errors)

    def to_pyarrow(self) -> object:
        """Return diagnostics as a PyArrow table."""
        try:
            pa = import_module("pyarrow")
        except ImportError as error:
            raise ImportError(
                "BatchDiagnostics.to_pyarrow() requires the optional 'pyarrow' package"
            ) from error
        return pa.table(self._column_mapping())

    def to_polars(self) -> object:
        """Return diagnostics as a Polars dataframe."""
        try:
            pl = import_module("polars")
        except ImportError as error:
            raise ImportError(
                "BatchDiagnostics.to_polars() requires the optional 'polars' package"
            ) from error
        return pl.DataFrame(self._column_mapping())

    def _column_mapping(self) -> dict[str, list[object]]:
        return {
            "source_index": [error.source_index for error in self._errors],
            "source": [error.source_name for error in self._errors],
            "code": [error.code for error in self._errors],
            "message": [error.message for error in self._errors],
            "line": [error.line for error in self._errors],
            "column": [error.column for error in self._errors],
            "byte_start": [error.byte_start for error in self._errors],
            "byte_end": [error.byte_end for error in self._errors],
        }


class ScanResult(Iterator[ReadResult]):
    """A bounded, deterministic iterator over native parse results."""

    __slots__ = (
        "_closed",
        "_errors",
        "_exhausted",
        "_inflight",
        "_native",
        "_on_error",
        "_schema",
        "_sources",
        "_workers",
    )

    def __init__(
        self,
        sources: Iterable[Source],
        *,
        category: str | None,
        columns: Sequence[str] | None,
        predicates: list[PredicateSpec],
        schema: str | None,
        workers: int,
        on_error: Literal["raise", "collect"],
    ) -> None:
        self._sources = iter(sources)
        self._workers = workers
        self._on_error = on_error
        self._schema = schema
        self._errors: list[BatchError] = []
        self._inflight = 0
        self._exhausted = False
        self._closed = False
        try:
            self._native = _core._NativeScan(
                workers,
                category,
                None if columns is None else list(columns),
                predicates,
                schema,
            )
        except ValueError as error:
            _raise_native_error(error)
        for _ in range(workers):
            if not self._submit_next():
                break

    @property
    def batches(self) -> ScanResult:
        """Return the single-pass batch iterator."""
        return self

    @property
    def errors(self) -> BatchDiagnostics:
        """Return an immutable snapshot of failures observed so far."""
        return BatchDiagnostics(self._errors)

    def __iter__(self) -> ScanResult:
        return self

    def __next__(self) -> ReadResult:
        while self._inflight:
            item = self._native.next()
            if item is None:
                self.close()
                raise RuntimeError("native scan ended before all submitted sources")
            source_index, native, fields = item
            self._inflight -= 1
            self._submit_next()
            if fields is None:
                if isinstance(native, _core._CifDocument):
                    return CifDocument(native, self._schema)
                if isinstance(native, _core._CifTable):
                    return CifTable(native)
                self.close()
                raise RuntimeError("native scan returned an unknown batch type")

            error = _batch_error(source_index, fields)
            self._errors.append(error)
            if self._on_error == "raise":
                self._native.cancel()
                self._closed = True
                raise error

        self.close()
        raise StopIteration

    def close(self) -> None:
        """Stop admitting input and join native workers without holding the GIL."""
        if not self._closed:
            self._native.close()
            self._closed = True

    def _submit_next(self) -> bool:
        if self._exhausted:
            return False
        try:
            source = next(self._sources)
        except StopIteration:
            self._exhausted = True
            return False
        source_name, content = _normalize_source(source)
        if content is None:
            self._native.submit_file(source_name)
        else:
            self._native.submit_bytes(source_name, content)
        self._inflight += 1
        return True
