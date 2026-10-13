"""Fetch and verify public structure-preparation test data from RCSB."""

from __future__ import annotations

import argparse
import gzip
import json
import re
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from benchmarks.corpus import CorpusFile, verify_content
from tools.fetch_pdb_corpus import atomic_write, download

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "tests" / "preparation_corpus.json"
DEFAULT_DESTINATION = ROOT / ".cache" / "structure-preparation-v1"


def file_spec(value: object) -> CorpusFile:
    """Read a pinned byte count and SHA-256 from the manifest."""
    if not isinstance(value, dict):
        raise ValueError("Corpus file record must be an object")
    size, digest = value.get("bytes"), value.get("sha256")
    if (
        isinstance(size, bool)
        or not isinstance(size, int)
        or size < 1
        or not isinstance(digest, str)
        or re.fullmatch(r"[0-9a-f]{64}", digest) is None
    ):
        raise ValueError("Corpus file record needs a byte count and SHA-256")
    return CorpusFile(size, digest)


def records(value: object, pattern: str) -> dict[str, CorpusFile]:
    """Check IDs before using them in cache paths or download URLs."""
    if not isinstance(value, dict) or not value:
        raise ValueError("Corpus manifest needs nonempty file records")
    result = {}
    for identifier, spec in value.items():
        if not isinstance(identifier, str) or re.fullmatch(pattern, identifier) is None:
            raise ValueError(f"Invalid corpus ID: {identifier!r}")
        result[identifier] = file_spec(spec)
    return result


class PreparationCorpus:
    """Download only requested data and verify each used cache file once."""

    def __init__(
        self,
        destination: Path = DEFAULT_DESTINATION,
        *,
        offline: bool = False,
        manifest: Path = DEFAULT_MANIFEST,
    ) -> None:
        raw = json.loads(manifest.read_text(encoding="utf-8"))
        if (
            not isinstance(raw, dict)
            or raw.get("schema") != "nibbler.preparation-corpus.v1"
        ):
            raise ValueError("Unsupported preparation corpus schema")
        self.entries = records(raw.get("entries"), r"[0-9a-z]{4}")
        self.component_records = records(raw.get("components"), r"[A-Z0-9]{1,10}")
        self.combined = file_spec(raw.get("combined_components"))
        self.destination = destination
        self.offline = offline
        self._verified: set[str] = set()

    def _fetch(
        self, file: Path, spec: CorpusFile, url: str, *, compressed: bool = False
    ) -> bytes:
        cached = file.is_file()
        if not cached and self.offline:
            raise FileNotFoundError(f"Missing cached test data: {file}")
        content = file.read_bytes() if cached else download(url)
        decoded = gzip.decompress(content) if compressed else content
        verify_content(decoded, spec, str(file))
        if not cached:
            atomic_write(file, content)
        return decoded

    def entry(self, pdb_id: str) -> Path:
        """Fetch one coordinate file; hash its uncompressed CIF content."""
        spec = self.entries[pdb_id]
        file = self.destination / f"{pdb_id}.cif.gz"
        if pdb_id not in self._verified:
            self._fetch(
                file,
                spec,
                f"https://files.rcsb.org/download/{pdb_id.upper()}.cif.gz",
                compressed=True,
            )
            self._verified.add(pdb_id)
        return file

    def _component(self, identifier: str) -> bytes:
        return self._fetch(
            self.destination / "ccd" / f"{identifier}.cif",
            self.component_records[identifier],
            f"https://files.rcsb.org/ligands/download/{identifier}.cif",
        )

    def components(self) -> Path:
        """Assemble the small CCD subset once, with at most four downloads at a time."""
        file = self.destination / "components.cif"
        if "components" not in self._verified:
            if file.is_file():
                verify_content(file.read_bytes(), self.combined, str(file))
            else:
                with ThreadPoolExecutor(max_workers=4) as pool:
                    content = b"\n".join(
                        pool.map(self._component, sorted(self.component_records))
                    )
                verify_content(content, self.combined, str(file))
                atomic_write(file, content)
            self._verified.add("components")
        return file


def main() -> None:
    """Fetch the full corpus, or verify an existing cache without network I/O."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--destination", type=Path, default=DEFAULT_DESTINATION)
    parser.add_argument("--offline", action="store_true")
    args = parser.parse_args()
    corpus = PreparationCorpus(args.destination, offline=args.offline)
    for pdb_id in corpus.entries:
        corpus.entry(pdb_id)
    corpus.components()
    print(f"Verified preparation test data in {args.destination}")


if __name__ == "__main__":
    main()
