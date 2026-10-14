"""Native crop plans must match the frozen Python spatial-crop law."""

import numpy as np
import pytest

from nibbler.observations import ObservationPacker

tokens_module = pytest.importorskip("kaveh_data.pdb_tokens")
native = pytest.importorskip("kaveh_data.native_observations")
schema = pytest.importorskip("kaveh_data.atom32_schema")


def packer():
    return ObservationPacker(schema.ENCODED, list(tokens_module.LIGAND_ELEMENTS), [])


def fixture(kind):
    rng = np.random.default_rng(731)
    sizes = {"protein": [700], "ligand": [500, 24], "partial": [500, 70],
             "fallback": [500, *([5] * 12)], "many-chains": [50] * 70,
             "ties": [500]}[kind]
    chains = np.repeat(np.arange(len(sizes)), sizes).astype(np.int64)
    count = len(chains)
    pos = np.zeros((count, 32, 3), dtype=np.float32)
    pos[:, 0] = rng.normal(size=(count, 3))
    if kind == "many-chains":
        pos[:, 0] *= .02
        pos[:, 0, 0] += chains * 5
    if kind == "ties":
        pos[:, 0] = np.arange(count)[:, None] % 11
    ligand = chains > 0 if kind in {"ligand", "partial", "fallback"} else np.zeros(count, dtype=bool)
    mask = np.ones((count, 32), dtype=bool)
    mask[rng.choice(np.flatnonzero(~ligand), 9, replace=False), 1] = False
    return dict(pos=pos, chain_idx=chains, is_ligand=ligand, mask=mask)


@pytest.mark.parametrize("size", [50, 384, 512])
@pytest.mark.parametrize("kind", ["protein", "ligand", "partial", "fallback", "many-chains", "ties"])
def test_native_spatial_plans_match_python(kind, size):
    tokens = fixture(kind)
    expected, expected_mass = tokens_module.crop_plans(tokens, crop_size=size)
    packed, mass = packer().crop_plans(np.ascontiguousarray(tokens["pos"][:, 0]).tobytes(),
        tokens["chain_idx"].tobytes(), tokens["is_ligand"].tobytes(),
        (tokens["is_ligand"] | tokens["mask"][:, 1]).tobytes(), crop_size=size)
    actual = native.arrays(packed)
    for name in ("crop_center", "crop_radius", "crop_tie", "crop_chain_bits"):
        np.testing.assert_array_equal(actual[name], expected[name].astype(actual[name].dtype).reshape(-1), err_msg=name)
    np.testing.assert_allclose(actual["crop_probability"], expected["crop_probability"], rtol=1e-13, atol=1e-15)
    assert mass == pytest.approx(expected_mass, rel=1e-13, abs=1e-15)
    np.testing.assert_array_equal(actual["crop_offsets"], [0, len(expected["crop_center"])])


def test_native_permutation_rejects_heavy_moves_and_cross_chain_swaps():
    before = np.array([[0., 0., 0.], [1., 0., 0.], [2., 0., 0.]], dtype=np.float32)
    after = before[[0, 2, 1]]
    chains = np.array([0, 1, 1], dtype=np.int64)
    ligand = np.array([False, True, True])
    worker = packer()
    assert worker.crop_permutation(before.tobytes(), after.tobytes(), chains.tobytes(), ligand.tobytes()) == [0, 2, 1]
    after[0, 0] += 1
    with pytest.raises(ValueError, match="protein crop center"):
        worker.crop_permutation(before.tobytes(), after.tobytes(), chains.tobytes(), ligand.tobytes())
    after[0] = before[0]
    chains[2] = 2
    with pytest.raises(ValueError, match="coordinate set"):
        worker.crop_permutation(before.tobytes(), after.tobytes(), chains.tobytes(), ligand.tobytes())
