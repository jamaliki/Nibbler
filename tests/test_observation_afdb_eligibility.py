"""Match packed AFDB eligibility to native complete-residue chain output."""

import json
import os
from pathlib import Path
import struct

import numpy as np
import pyarrow as pa
import pytest
from token_tables import token_table

observations = pytest.importorskip("nibbler.observations")
if os.environ.get("KAVEH_ATOM32_SCHEMA"):
    ENCODED = Path(os.environ["KAVEH_ATOM32_SCHEMA"]).read_text()
    SPEC = json.loads(ENCODED)
else:
    schema = pytest.importorskip("kaveh_data.atom32_schema")
    ENCODED, SPEC = schema.ENCODED, schema.SPEC


def arrays(buffers):
    return {name: np.frombuffer(buffer, dtype=dtype).reshape(shape)
            for name, (dtype, shape, buffer) in buffers.items()}


def decode(raw):
    _, _, _, _, _, count, residues = struct.unpack_from("<4sBBBBII", raw)
    offsets = np.frombuffer(raw, dtype="<u4", count=count + 1, offset=16)
    start = 16 + 4 * (count + 1)
    sequence = raw[start:start + residues].decode()
    aa = np.array(["ARNDCQEGHILKMFPSTWYV".index(letter) for letter in sequence], dtype=np.int64)
    xyz = np.frombuffer(raw, dtype=np.uint8, offset=start + residues)
    xyz = xyz.reshape(4, -1).T.copy().view("<f4").reshape(residues, 14, 3)
    physical = np.array([[bool(name) for name in SPEC["residues"][identity]["names"][:14]] for identity in aa])
    finite = np.isfinite(xyz)
    if (np.isinf(xyz).any() or np.any(xyz[~physical] != 0)
            or np.any(finite.any(-1) != finite.all(-1))):
        raise ValueError("Invalid packed coordinates")
    mask = physical & finite.all(-1)
    if not mask[:, 1].all():
        raise ValueError("Missing CA")
    xyz[~mask] = 0
    return dict(pos=xyz, mask=mask, aatype=aa, sample_offsets=offsets)


def protein(length, sequence=None):
    sequence = "A" * length if sequence is None else sequence
    xyz = np.zeros((length, 14, 3), dtype="<f4")
    for row, letter in enumerate(sequence):
        aa = "ARNDCQEGHILKMFPSTWYV".index(letter)
        names = SPEC["residues"][aa]["names"][:14]
        for slot, name in enumerate(names):
            if name:
                xyz[row, slot] = [row * 3.8 + 1.4, slot * .1, 0]
        xyz[row, 0] = [row * 3.8, 0, 0]
        xyz[row, 2] = [row * 3.8 + 2.5, 0, 0]
    return sequence, xyz


def encode(proteins):
    lengths = [len(sequence) for sequence, _ in proteins]
    coords = np.concatenate([xyz for _, xyz in proteins])
    return (struct.pack("<4sBBBBII", b"A14P", 1, 32, 1, 0, len(lengths), sum(lengths))
            + np.asarray([0, *np.cumsum(lengths)], dtype="<u4").tobytes()
            + "".join(sequence for sequence, _ in proteins).encode()
            + coords.view(np.uint8).reshape(-1, 4).T.tobytes())


def native():
    return observations.ObservationPacker(ENCODED, token_table(), [])


def scan(raw, **kwargs):
    buffer = memoryview(pa.py_buffer(raw)).toreadonly()
    return arrays(native().afdb_eligibility(buffer, **kwargs))


def reference(raw):
    decoded = decode(raw)
    worker = native()
    result = {key: [] for key in ("eligible", "complete_residues", "max_chain_residues", "reconstructed_chains")}
    for index, (start, end) in enumerate(zip(decoded["sample_offsets"][:-1], decoded["sample_offsets"][1:])):
        document, chains, _ = worker.protein_document(decoded["pos"][start:end].tobytes(),
            decoded["mask"][start:end].tobytes(), decoded["aatype"][start:end].tobytes(), f"AF-{index}")
        buffers, _ = worker.pack(document, protein_only=True, complete_tokens=True, validate=True)
        tokens = arrays(buffers)
        _, counts = np.unique(tokens["chain_idx"], return_counts=True)
        maximum = int(counts.max(initial=0))
        result["eligible"].append(maximum >= 50)
        result["complete_residues"].append(len(tokens["aatype"]))
        result["max_chain_residues"].append(maximum)
        result["reconstructed_chains"].append(chains)
    return result


def test_full_packed_shard_matches_native_document_and_complete_chain_counts():
    sequence = ("ARNDCQEGHILKMFPSTWYV" * 20)[:383]
    raw = encode([protein(33), protein(49), protein(50), protein(383, sequence)])
    actual, expected = scan(raw), reference(raw)
    for key in actual:
        np.testing.assert_array_equal(actual[key], expected[key], err_msg=key)
    np.testing.assert_array_equal(actual["eligible"], [False, False, True, True])
    assert actual["eligible"].dtype == bool
    assert all(value.dtype == np.int64 for key, value in actual.items() if key != "eligible")


@pytest.mark.parametrize("case", ["sidechain_gap", "backbone_gap", "split_chains", "all_incomplete", "exact_crop"])
def test_missing_atoms_and_chain_breaks_match_complete_token_output(case):
    sequence, xyz = protein(384 if case == "exact_crop" else 60)
    if case == "sidechain_gap":
        xyz[::6, 4] = np.nan
    elif case == "backbone_gap":
        xyz[30, 2] = np.nan
    elif case == "split_chains":
        xyz[30:, :5, 0] += 20
    elif case == "all_incomplete":
        xyz[:, 4] = np.nan
    raw = encode([(sequence, xyz)])
    actual, expected = scan(raw), reference(raw)
    for key in actual:
        np.testing.assert_array_equal(actual[key], expected[key], err_msg=key)
    if case == "sidechain_gap":
        assert actual["eligible"].tolist() == [True]
        assert actual["max_chain_residues"].tolist() == [50]
        assert actual["reconstructed_chains"].tolist() == [1]


@pytest.mark.parametrize("value", [.8, 1.8, np.nextafter(np.float32(.8), np.float32(1)),
                                  np.nextafter(np.float32(.8), np.float32(0)),
                                  np.nextafter(np.float32(1.8), np.float32(1)),
                                  np.nextafter(np.float32(1.8), np.float32(2))])
def test_connectivity_uses_native_float32_strict_distance_bounds(value):
    sequence, xyz = protein(60)
    xyz[29, 2] = [0, 0, 0]
    xyz[30, 0] = [value, 0, 0]
    raw = encode([(sequence, xyz)])
    actual, expected = scan(raw), reference(raw)
    for key in actual:
        np.testing.assert_array_equal(actual[key], expected[key], err_msg=key)


@pytest.mark.parametrize("bits", [0, 0x80000000, 1, 0x80000001, 0x007fffff, 0x00800000,
                                 0x3f800001, 0xbf800001, 0x7f7fffff, 0xff7fffff,
                                 0x7fc00001, 0xffc01234, 0x7f800001])
def test_exponent_finite_checks_match_decoder_and_native_tokens(bits):
    sequence, xyz = protein(50)
    xyz[0, 4] = np.uint32(bits).view(np.float32)
    raw = encode([(sequence, xyz)])
    actual, expected = scan(raw), reference(raw)
    for key in actual:
        np.testing.assert_array_equal(actual[key], expected[key], err_msg=key)
    assert bool(actual["eligible"][0]) == bool(np.isfinite(np.uint32(bits).view(np.float32)))


@pytest.mark.parametrize("case", ["infinite", "negative_infinite", "partial_nan", "missing_ca", "padding"])
def test_invalid_coordinate_fields_fail_like_the_full_decoder(case):
    sequence, xyz = protein(50)
    if case in {"infinite", "negative_infinite"}:
        xyz[0, 4] = np.inf if case == "infinite" else -np.inf
    elif case == "partial_nan":
        xyz[0, 4, 0] = np.nan
    elif case == "missing_ca":
        xyz[0, 1] = np.nan
    else:
        xyz[0, 13, 0] = 1
    raw = encode([(sequence, xyz)])
    with pytest.raises(ValueError):
        decode(raw)
    with pytest.raises(ValueError, match="A14P sample 0 residue 0"):
        scan(raw)


def test_unsupported_spatial_crop_lengths_fail_instead_of_filtering():
    with pytest.raises(ValueError, match="length 385 exceeds crop_size 384"):
        scan(encode([protein(385)]))


@pytest.mark.parametrize("unsigned", [False, True])
def test_scan_reads_the_arrow_buffer_without_mutating_or_retaining_it(unsigned):
    raw = encode([protein(50)])
    buffer = pa.Codec("zstd").decompress(pa.Codec("zstd").compress(raw), len(raw))
    view = memoryview(buffer).toreadonly()
    if unsigned:
        view = view.cast("B")
    assert any(byte > 127 for byte in raw)
    assert pa.py_buffer(view).address == buffer.address
    result = arrays(native().afdb_eligibility(view))
    assert view.cast("B") == raw
    del view, buffer
    assert result["eligible"].tolist() == [True]


@pytest.mark.parametrize("buffer", [lambda raw: bytearray(raw), lambda raw: memoryview(raw)[::2]])
def test_scan_rejects_mutable_or_noncontiguous_buffers(buffer):
    with pytest.raises(ValueError, match="read-only contiguous"):
        native().afdb_eligibility(buffer(encode([protein(50)])))
