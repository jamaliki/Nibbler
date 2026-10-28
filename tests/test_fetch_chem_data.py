"""The pinned chem_data subset downloads from fixed sources and is hash-checked."""

import hashlib
import json

import pytest

from tools import fetch_chem_data as data

GEOSTD = "https://raw.githubusercontent.com/phenix-project/geostd/" + "a" * 40 + "/"
MON_LIB = "https://raw.githubusercontent.com/MonomerLibrary/monomers/" + "b" * 40 + "/"
CCD = "https://files.rcsb.org/ligands/download/"


def record(content):
    return {"bytes": len(content), "sha256": hashlib.sha256(content).hexdigest()}


CONTENT = {
    "geostd/list/geostd_list.cif": b"data_geostd_list\n",
    "mon_lib/a/ALA.cif": b"data_comp_ALA\n",
    "chemical_components/f/data_FES.cif": b"data_FES\n",
}


def write_manifest(file, files=None, **sources):
    manifest = {
        "schema": data.SCHEMA,
        "sources": {
            "chemical_components": CCD,
            "geostd": GEOSTD,
            "mon_lib": MON_LIB,
            **sources,
        },
        "files": files
        if files is not None
        else {path: record(content) for path, content in CONTENT.items()},
    }
    file.write_text(json.dumps(manifest))
    return file


@pytest.fixture
def subset(tmp_path, monkeypatch):
    responses = {
        GEOSTD + "list/geostd_list.cif": CONTENT["geostd/list/geostd_list.cif"],
        MON_LIB + "a/ALA.cif": CONTENT["mon_lib/a/ALA.cif"],
        CCD + "FES.cif": CONTENT["chemical_components/f/data_FES.cif"],
    }
    requests = []

    def download(url):
        requests.append(url)
        return responses[url]

    monkeypatch.setattr(data, "download", download)
    manifest = write_manifest(tmp_path / "manifest.json")
    return data.Subset(manifest), manifest, responses, requests


def test_each_file_comes_from_its_pinned_source(subset, tmp_path):
    files, _, responses, requests = subset
    assert data.fetch(files, tmp_path / "cache") == 3
    assert sorted(requests) == sorted(responses)
    for path, content in CONTENT.items():
        assert (tmp_path / "cache" / path).read_bytes() == content


def test_a_verified_cache_needs_no_network(subset, tmp_path):
    files, _, _, requests = subset
    data.fetch(files, tmp_path / "cache")
    requests.clear()
    assert data.fetch(files, tmp_path / "cache", offline=True) == 0
    assert not requests


def test_offline_fetch_fails_on_a_changed_cache(subset, tmp_path):
    files, _, _, _ = subset
    data.fetch(files, tmp_path / "cache")
    (tmp_path / "cache" / "mon_lib/a/ALA.cif").write_bytes(b"changed\n")
    with pytest.raises(FileNotFoundError, match=r"ALA\.cif"):
        data.fetch(files, tmp_path / "cache", offline=True)
    assert data.fetch(files, tmp_path / "cache") == 1


def test_a_changed_download_is_never_written(subset, tmp_path):
    files, _, responses, _ = subset
    responses[CCD + "FES.cif"] = b"data_FES\nrevised\n"
    with pytest.raises(RuntimeError, match=r"data_FES\.cif"):
        data.fetch(files, tmp_path / "cache")
    assert not (tmp_path / "cache/chemical_components/f/data_FES.cif").exists()


def test_repin_records_the_current_downloads(subset):
    files, manifest, responses, _ = subset
    revised = b"data_FES\nrevised\n"
    responses[CCD + "FES.cif"] = revised
    assert data.repin(files, manifest) == ["chemical_components/f/data_FES.cif"]
    pinned = json.loads(manifest.read_text())["files"]
    assert pinned["chemical_components/f/data_FES.cif"] == record(revised)
    assert data.Subset(manifest).files.keys() == files.files.keys()


@pytest.mark.parametrize(
    "path",
    [
        "geostd/../outside.cif",
        "/geostd/a/data_ALA.cif",
        "other/a/ALA.cif",
        "chemical_components/a/data_FES.cif",
        "chemical_components/f/FES.cif",
    ],
)
def test_paths_outside_the_layout_are_rejected(tmp_path, path):
    manifest = write_manifest(tmp_path / "manifest.json", {path: record(b"x")})
    with pytest.raises(ValueError, match="path"):
        data.Subset(manifest)


def test_dictionary_sources_must_name_a_commit(tmp_path):
    manifest = write_manifest(
        tmp_path / "manifest.json",
        geostd="https://raw.githubusercontent.com/phenix-project/geostd/master/",
    )
    with pytest.raises(ValueError, match="commit"):
        data.Subset(manifest)


def test_the_repository_manifest_is_valid():
    subset = data.Subset()
    assert len(subset.files) == 200
    assert all(
        path.split("/")[0] in data.SOURCES and subset.url(path) for path in subset.files
    )
