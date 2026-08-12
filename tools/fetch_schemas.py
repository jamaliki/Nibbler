"""Fetch exact schema inputs declared in the repository lock file."""

from __future__ import annotations

import hashlib
import shutil
import sys
import tempfile
import urllib.request
from pathlib import Path
from typing import TypedDict, cast

if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib

ROOT = Path(__file__).resolve().parents[1]
LOCK_FILE = ROOT / "schemas" / "locks.toml"
DEFAULT_CACHE = ROOT / ".cache" / "schemas"


class SchemaLock(TypedDict):
    """One serialized schema lock."""

    name: str
    dictionary_name: str
    version: str
    file: str
    url: str
    sha256: str
    retrieved_on: str
    artifact: str
    extensions: list[str]


class SchemaLocks(TypedDict):
    """Serialized schema-lock document."""

    format_version: int
    artifact_format_version: int
    schema: list[SchemaLock]


def load_schema_locks() -> SchemaLocks:
    """Load the schema lock file."""
    raw = tomllib.loads(LOCK_FILE.read_text(encoding="utf-8"))
    return cast(SchemaLocks, raw)


def file_digest(file: Path) -> str:
    """Calculate a file's SHA-256 digest without loading it all at once."""
    digest = hashlib.sha256()
    with file.open("rb") as input_file:
        for chunk in iter(lambda: input_file.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fetch_schema(lock: SchemaLock, cache_dir: Path = DEFAULT_CACHE) -> Path:
    """Return a verified cached dictionary, fetching it when necessary."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    destination = cache_dir / lock["file"]
    if destination.is_file() and file_digest(destination) == lock["sha256"]:
        return destination

    request = urllib.request.Request(
        lock["url"], headers={"User-Agent": "nibbler-schema-bootstrap/0.1"}
    )
    temporary_file = tempfile.NamedTemporaryFile(
        mode="wb", dir=cache_dir, prefix=f".{lock['file']}.", delete=False
    )
    temporary_name = Path(temporary_file.name)
    try:
        with temporary_file, urllib.request.urlopen(request, timeout=120) as response:
            shutil.copyfileobj(response, temporary_file)
        actual_digest = file_digest(temporary_name)
        if actual_digest != lock["sha256"]:
            raise RuntimeError(
                f"digest mismatch for {lock['name']}: "
                f"expected {lock['sha256']}, found {actual_digest}"
            )
        temporary_name.replace(destination)
    except BaseException:
        temporary_name.unlink(missing_ok=True)
        raise
    return destination


def fetch_all(cache_dir: Path = DEFAULT_CACHE) -> dict[str, Path]:
    """Fetch and verify every locked schema."""
    locks = load_schema_locks()
    if locks["format_version"] != 1:
        raise RuntimeError(
            f"unsupported schema lock version: {locks['format_version']}"
        )
    return {lock["name"]: fetch_schema(lock, cache_dir) for lock in locks["schema"]}


def main() -> int:
    """Fetch all schemas and report verified cache locations."""
    for name, file in fetch_all().items():
        print(f"{name}: {file}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
