"""Fetch the pinned native validation controls outside the source tree."""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import urllib.request
import zipfile
from pathlib import Path, PurePosixPath

from tools.fetch_pdb_corpus import atomic_write

DEFAULT_DESTINATION = Path(".cache/native-validation-v1")
URL = (
    "https://github.com/jamaliki/Nibbler/releases/download/"
    "native-validation-data-20261012-v1/nibbler-native-validation-data-v1.zip"
)
SHA256 = "75238538098acd661091ae28353383be0861b6538477d240a4bcdd3db0e3e5f3"
SIZE = 3064261


def install(content: bytes, destination: Path) -> None:
    """Verify the bundle before writing files, then publish its manifest last."""
    if len(content) != SIZE or hashlib.sha256(content).hexdigest() != SHA256:
        raise ValueError("Validation data does not match the pinned size and SHA256")
    with zipfile.ZipFile(io.BytesIO(content)) as archive:
        entries = archive.infolist()
        names = [PurePosixPath(entry.filename) for entry in entries]
        if (
            len(set(names)) != len(names)
            or any(
                name.is_absolute()
                or ".." in name.parts
                or "\\" in entry.filename
                or entry.is_dir()
                or (entry.external_attr >> 16) & 0o170000 == 0o120000
                for name, entry in zip(names, entries, strict=True)
            )
            or sum(entry.file_size for entry in entries) > 64 * 1024 * 1024
        ):
            raise ValueError("Validation data contains invalid archive members")
        records = {
            str(name): archive.read(entry)
            for name, entry in zip(names, entries, strict=True)
        }
    manifest = json.loads(records.pop("manifest.json"))
    if manifest["schema"] != "nibbler.validation-data.v1":
        raise ValueError("Validation data has an unsupported schema")
    for name in manifest["files"].values():
        if hashlib.sha256(records[name]).hexdigest() != Path(name).name.split(".")[0]:
            raise ValueError("Validation data contains a changed fixture")
    for name, data in records.items():
        atomic_write(destination / name, data)
    atomic_write(destination / "manifest.json", json.dumps(manifest).encode())


def main() -> None:
    """Download once, or install an already downloaded bundle without network I/O."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--destination", type=Path, default=DEFAULT_DESTINATION)
    parser.add_argument("--archive", type=Path)
    args = parser.parse_args()
    if args.archive:
        content = args.archive.read_bytes()
    else:
        request = urllib.request.Request(URL, headers={"User-Agent": "Nibbler/0.1"})
        with urllib.request.urlopen(request, timeout=120) as response:
            content = response.read(SIZE + 1)
    install(content, args.destination)
    print(f"Verified native validation data in {args.destination}")


if __name__ == "__main__":
    main()
