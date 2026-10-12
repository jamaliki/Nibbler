"""Keep validation data downloads pinned and archive writes contained."""

import hashlib
import io
import json
import zipfile

import pytest

from tools import fetch_validation_data as bundle


def archive(monkeypatch, extra=None):
    content = b"data_control\n"
    name = f"files/{hashlib.sha256(content).hexdigest()}.cif"
    records = {
        name: content,
        "manifest.json": json.dumps(
            {"schema": "nibbler.validation-data.v1", "files": {"control": name}}
        ).encode(),
        **(extra or {}),
    }
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as output:
        for key, value in records.items():
            output.writestr(key, value)
    data = buffer.getvalue()
    monkeypatch.setattr(bundle, "SIZE", len(data))
    monkeypatch.setattr(bundle, "SHA256", hashlib.sha256(data).hexdigest())
    return data, name


def test_valid_bundle_keeps_source_bytes_and_publishes_manifest(monkeypatch, tmp_path):
    data, name = archive(monkeypatch)
    bundle.install(data, tmp_path)
    assert (tmp_path / name).read_bytes() == b"data_control\n"
    assert (
        json.loads((tmp_path / "manifest.json").read_text())["files"]["control"] == name
    )


@pytest.mark.parametrize("change", [b"", b"extra", b"different bytes"])
def test_changed_bundle_cannot_write_files(change, monkeypatch, tmp_path):
    data, _ = archive(monkeypatch)
    with pytest.raises(ValueError, match="pinned size and SHA256"):
        bundle.install(data[:-1] + change, tmp_path)
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize("name", ["../outside", "/absolute", "folder\\outside"])
def test_archive_paths_cannot_escape_destination(name, monkeypatch, tmp_path):
    data, _ = archive(monkeypatch, {name: b"bad"})
    with pytest.raises(ValueError, match="invalid archive members"):
        bundle.install(data, tmp_path)
    assert list(tmp_path.iterdir()) == []


def test_changed_fixture_cannot_publish_partial_bundle(monkeypatch, tmp_path):
    _, name = archive(monkeypatch)
    data, _ = archive(monkeypatch, {name: b"changed"})
    with pytest.raises(ValueError, match="changed fixture"):
        bundle.install(data, tmp_path)
    assert list(tmp_path.iterdir()) == []


def test_symlink_members_cannot_write_files(monkeypatch, tmp_path):
    member = zipfile.ZipInfo("link")
    member.create_system = 3
    member.external_attr = 0o120777 << 16
    data, _ = archive(monkeypatch, {member: b"../outside"})
    with pytest.raises(ValueError, match="invalid archive members"):
        bundle.install(data, tmp_path)
    assert list(tmp_path.iterdir()) == []


def test_normalized_duplicate_names_cannot_write_files(monkeypatch, tmp_path):
    data, _ = archive(monkeypatch, {"./manifest.json": b"duplicate"})
    with pytest.raises(ValueError, match="invalid archive members"):
        bundle.install(data, tmp_path)
    assert list(tmp_path.iterdir()) == []
