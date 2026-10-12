"""Public test data downloads stay lazy, hash-checked, and outside Git."""

import gzip
import hashlib
import json

import pytest

from benchmarks.linked_products import CASES
from tools import fetch_preparation_data as data


def record(content):
    return {"bytes": len(content), "sha256": hashlib.sha256(content).hexdigest()}


@pytest.fixture
def corpus(tmp_path, monkeypatch):
    entry = b"data_entry\n_entry.id 1abc\n"
    first, second = b"data_A\n_chem_comp.id A\n", b"data_Z\n_chem_comp.id Z\n"
    manifest = {
        "schema": "nibbler.preparation-corpus.v1",
        "entries": {"1abc": record(entry)},
        "components": {"Z": record(second), "A": record(first)},
        "combined_components": record(first + b"\n" + second),
    }
    file = tmp_path / "manifest.json"
    file.write_text(json.dumps(manifest))
    responses = {
        "https://files.rcsb.org/download/1ABC.cif.gz": gzip.compress(entry),
        "https://files.rcsb.org/ligands/download/A.cif": first,
        "https://files.rcsb.org/ligands/download/Z.cif": second,
    }
    requests = []

    def download(url):
        requests.append(url)
        return responses[url]

    monkeypatch.setattr(data, "download", download)
    return data.PreparationCorpus(tmp_path / "cache", manifest=file), requests


def test_constructor_does_not_fetch_or_create_a_cache(corpus):
    cache, requests = corpus
    assert not requests and not cache.destination.exists()


def test_entry_download_is_verified_once_per_session(corpus, monkeypatch):
    cache, requests = corpus
    file = cache.entry("1abc")
    assert gzip.decompress(file.read_bytes()) == b"data_entry\n_entry.id 1abc\n"

    def unexpected_verification(*args):
        pytest.fail("Repeated access should not reread and hash the same entry")

    monkeypatch.setattr(data, "verify_content", unexpected_verification)
    assert cache.entry("1abc") == file
    assert len(requests) == 1


def test_warm_cache_works_offline(corpus, tmp_path):
    cache, requests = corpus
    file = cache.entry("1abc")
    components = cache.components()
    repeat = data.PreparationCorpus(
        cache.destination, offline=True, manifest=tmp_path / "manifest.json"
    )
    assert repeat.entry("1abc") == file
    assert repeat.components() == components
    assert (
        components.read_bytes()
        == b"data_A\n_chem_comp.id A\n\ndata_Z\n_chem_comp.id Z\n"
    )
    assert len(requests) == 3


def test_gzip_metadata_does_not_change_the_cif_fingerprint(corpus, monkeypatch):
    cache, _ = corpus
    monkeypatch.setattr(
        data,
        "download",
        lambda url: gzip.compress(b"data_entry\n_entry.id 1abc\n", mtime=123),
    )
    assert cache.entry("1abc").is_file()


def test_changed_download_does_not_enter_the_cache(corpus, monkeypatch):
    cache, _ = corpus
    monkeypatch.setattr(data, "download", lambda url: gzip.compress(b"different"))
    with pytest.raises(RuntimeError, match=r"expected .* bytes"):
        cache.entry("1abc")
    assert not cache.destination.exists()


def test_corrupt_cached_content_fails_without_redownloading(corpus, tmp_path):
    cache, requests = corpus
    cache.entry("1abc").write_bytes(gzip.compress(b"data_entry\n_entry.id 1abd\n"))
    repeat = data.PreparationCorpus(
        cache.destination, manifest=tmp_path / "manifest.json"
    )
    with pytest.raises(RuntimeError, match="expected sha256"):
        repeat.entry("1abc")
    assert len(requests) == 1


def test_missing_offline_data_fails_without_network(corpus):
    cache, requests = corpus
    cache.offline = True
    with pytest.raises(FileNotFoundError, match="Missing cached test data"):
        cache.entry("1abc")
    assert not requests


def test_unknown_entry_cannot_construct_a_download_url(corpus):
    cache, requests = corpus
    with pytest.raises(KeyError):
        cache.entry("../outside")
    assert not requests


def test_bad_component_does_not_publish_a_combined_ccd(corpus, monkeypatch):
    cache, _ = corpus
    monkeypatch.setattr(data, "download", lambda url: b"changed CCD")
    with pytest.raises(RuntimeError, match=r"expected .* bytes"):
        cache.components()
    assert not (cache.destination / "components.cif").exists()


@pytest.mark.parametrize("kind", ["entries", "components"])
def test_manifest_ids_cannot_escape_cache_paths(corpus, tmp_path, kind):
    file = tmp_path / "manifest.json"
    manifest = json.loads(file.read_text())
    manifest[kind] = {"../outside": record(b"data")}
    file.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="Invalid corpus ID"):
        data.PreparationCorpus(tmp_path / "cache", manifest=file)


def test_pinned_manifest_covers_product_and_structure_cases():
    cache = data.PreparationCorpus()
    assert set(CASES) <= cache.entries.keys()
    assert {"1cbn", "6lzb", "1rcq", "6qhp", "7mub", "1xan"} <= cache.entries.keys()
