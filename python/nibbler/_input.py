"""Source and predicate normalization for the thin Python facade."""

from __future__ import annotations

import os
from collections.abc import Collection, Mapping
from os import PathLike
from typing import BinaryIO, TypeAlias

from .contracts import MissingKind

Source: TypeAlias = str | PathLike[str] | bytes | bytearray | memoryview | BinaryIO
PredicateSpec: TypeAlias = tuple[str, str, list[str]]


def _normalize_source(source: Source) -> tuple[str, bytes | None]:
    if isinstance(source, str) or isinstance(source, os.PathLike):
        file = os.fspath(source)
        if not isinstance(file, str):
            raise TypeError("filesystem sources must resolve to a string path")
        return file, None
    if isinstance(source, bytes):
        return "<memory>", source
    if isinstance(source, (bytearray, memoryview)):
        return "<memory>", bytes(source)

    read_method = getattr(source, "read", None)
    if read_method is None:
        raise TypeError(
            "source must be a path, bytes-like object, or binary file object"
        )
    content = read_method()
    if not isinstance(content, (bytes, bytearray, memoryview)):
        raise TypeError("binary file source.read() must return bytes-like data")
    name = getattr(source, "name", "<binary-stream>")
    source_name = os.fspath(name) if isinstance(name, (str, os.PathLike)) else str(name)
    if not isinstance(source_name, str):
        source_name = "<binary-stream>"
    return source_name, bytes(content)


def _normalize_predicates(
    where: Mapping[str, object] | None,
) -> list[PredicateSpec]:
    if where is None:
        return []
    predicates: list[PredicateSpec] = []
    for column, condition in where.items():
        if not isinstance(column, str):
            raise TypeError("where keys must be column names")
        operator, operands = _normalize_condition(condition)
        predicates.append((column, operator, operands))
    return predicates


def _normalize_condition(condition: object) -> tuple[str, list[str]]:
    if isinstance(condition, MissingKind):
        return "missing", [condition.value]
    if condition is None:
        return "missing", []
    if isinstance(condition, tuple) and len(condition) == 2:
        operator, operand = condition
        if operator in {"==", "eq"}:
            return "eq", [_predicate_text(operand)]
        if operator in {"!=", "ne"}:
            return "ne", [_predicate_text(operand)]
        if operator == "in":
            return "in", _membership_values(operand)
        if operator in {"is", "missing"}:
            if operand is None:
                return "missing", []
            if isinstance(operand, MissingKind):
                return "missing", [operand.value]
        raise ValueError(f"unsupported native predicate {condition!r}")
    if isinstance(condition, Collection) and not isinstance(
        condition, (str, bytes, bytearray, memoryview)
    ):
        return "in", sorted(_predicate_text(value) for value in condition)
    return "eq", [_predicate_text(condition)]


def _membership_values(value: object) -> list[str]:
    if not isinstance(value, Collection) or isinstance(
        value, (str, bytes, bytearray, memoryview)
    ):
        raise TypeError("an 'in' predicate requires a non-string collection")
    return sorted(_predicate_text(item) for item in value)


def _predicate_text(value: object) -> str:
    if value is None or isinstance(value, MissingKind):
        raise TypeError("missing values require an 'is' or 'missing' predicate")
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)
