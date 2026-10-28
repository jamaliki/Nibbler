"""Fetch and verify the pinned Reduce3 chem_data subset that the tests read.

Reduce3 reads a ``chem_data`` directory with ``geostd/``, ``mon_lib/`` and
``chemical_components/``. The tests need only the components of their fixtures and of
the native validation data, so ``tests/chem_data_subset.json`` pins those files by size
and SHA-256: GeoStd and the Monomer Library at upstream commits, and CCD definitions
from RCSB. RCSB serves only the current revision of a component, so a revised component
stops matching; ``--repin`` downloads every file again and records the new sizes and
digests.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from collections.abc import Sequence
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path, PurePosixPath

from benchmarks.corpus import CorpusFile, verify_content
from tools.fetch_pdb_corpus import atomic_write, download

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "tests" / "chem_data_subset.json"
DEFAULT_DESTINATION = ROOT / ".cache" / "chem-data-v1"
SCHEMA = "nibbler.chem-data-subset.v1"
SOURCES = {
    "chemical_components": "https://files.rcsb.org/ligands/download/",
    "geostd": "https://raw.githubusercontent.com/phenix-project/geostd/",
    "mon_lib": "https://raw.githubusercontent.com/MonomerLibrary/monomers/",
}
PATH = re.compile(r"(chemical_components|geostd|mon_lib)(/[A-Za-z0-9_.-]+)+")
COMPONENT = re.compile(r"chemical_components/[0-9a-z]/data_([A-Z0-9]{1,5})\.cif")


class Subset:
    """The pinned files of the subset and the URL of each."""

    def __init__(self, manifest: Path = DEFAULT_MANIFEST) -> None:
        raw = json.loads(manifest.read_text(encoding="utf-8"))
        if not isinstance(raw, dict) or raw.get("schema") != SCHEMA:
            raise ValueError("Unsupported chem_data subset schema")
        sources = raw.get("sources")
        if not isinstance(sources, dict) or set(sources) != set(SOURCES):
            raise ValueError("The chem_data subset must name each source once")
        for name, prefix in SOURCES.items():
            url = sources[name]
            if not isinstance(url, str) or not url.startswith(prefix):
                raise ValueError(f"Unexpected {name} source: {url!r}")
            if name != "chemical_components" and not re.fullmatch(
                re.escape(prefix) + r"[0-9a-f]{40}/", url
            ):
                raise ValueError(f"The {name} source must name a commit: {url!r}")
        files = raw.get("files")
        if not isinstance(files, dict) or not files:
            raise ValueError("The chem_data subset needs nonempty file records")
        self.sources: dict[str, str] = sources
        self.files: dict[str, CorpusFile] = {}
        for path, record in files.items():
            self.url(path)
            if not isinstance(record, dict):
                raise ValueError(f"{path}: file record must be an object")
            size, digest = record.get("bytes"), record.get("sha256")
            if (
                isinstance(size, bool)
                or not isinstance(size, int)
                or size < 1
                or not isinstance(digest, str)
                or re.fullmatch(r"[0-9a-f]{64}", digest) is None
            ):
                raise ValueError(f"{path}: file record needs a byte count and SHA-256")
            self.files[path] = CorpusFile(size, digest)

    def url(self, path: object) -> str:
        """Check a manifest path before using it in a cache path or a download URL."""
        if (
            not isinstance(path, str)
            or PATH.fullmatch(path) is None
            or any(part in {".", ".."} for part in PurePosixPath(path).parts)
        ):
            raise ValueError(f"Invalid chem_data path: {path!r}")
        source, rest = path.split("/", 1)
        if source != "chemical_components":
            return self.sources[source] + rest
        component = COMPONENT.fullmatch(path)
        if component is None or path.split("/")[1] != component[1][0].lower():
            raise ValueError(f"Invalid chemical component path: {path!r}")
        return self.sources[source] + f"{component[1]}.cif"


def fetch(subset: Subset, destination: Path, *, offline: bool = False) -> int:
    """Download each file that the cache lacks or holds changed; return the count."""

    def one(path: str) -> bool:
        file, spec = destination / path, subset.files[path]
        if file.is_file():
            content = file.read_bytes()
            if (
                len(content) == spec.bytes
                and hashlib.sha256(content).hexdigest() == spec.sha256
            ):
                return False
        if offline:
            raise FileNotFoundError(f"Missing or changed cached chem_data: {file}")
        atomic_write(file, verify_content(download(subset.url(path)), spec, path))
        return True

    with ThreadPoolExecutor(max_workers=8) as pool:
        return sum(pool.map(one, sorted(subset.files)))


def repin(subset: Subset, manifest: Path) -> list[str]:
    """Record the current size and digest of every file; return the changed paths."""
    paths = sorted(subset.files)
    with ThreadPoolExecutor(max_workers=8) as pool:
        contents = list(pool.map(download, map(subset.url, paths)))
    raw = json.loads(manifest.read_text(encoding="utf-8"))
    changed = []
    for path, content in zip(paths, contents, strict=True):
        record = {"bytes": len(content), "sha256": hashlib.sha256(content).hexdigest()}
        if raw["files"][path] != record:
            raw["files"][path] = record
            changed.append(path)
    atomic_write(manifest, (json.dumps(raw, indent=2) + "\n").encode())
    return changed


def main(argv: Sequence[str] | None = None) -> None:
    """Fetch the subset, verify a cache without network I/O, or re-pin the manifest."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--destination", type=Path, default=DEFAULT_DESTINATION)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--offline", action="store_true")
    mode.add_argument("--repin", action="store_true")
    args = parser.parse_args(argv)
    subset = Subset(args.manifest)
    if args.repin:
        for path in repin(subset, args.manifest):
            print(f"re-pinned {path}")
        return
    count = fetch(subset, args.destination, offline=args.offline)
    print(
        f"Verified {len(subset.files)} chem_data files in {args.destination} "
        f"({count} downloaded)"
    )


if __name__ == "__main__":
    main()
