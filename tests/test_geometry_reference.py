"""Compare native primitives against CCTBX, without chemistry interpretation."""

import json
import math
import random
from pathlib import Path

import pytest

from nibbler.geometry import RestraintGraph

restraints = pytest.importorskip("cctbx.geometry_restraints")
flex = pytest.importorskip("scitbx.array_family.flex")


def graph(sites, *, bonds=(), angles=(), chiralities=(), planes=()):
    return RestraintGraph(
        [str(i) for i in range(len(sites))],
        [False] * len(sites),
        list(bonds),
        list(angles),
        list(chiralities),
        list(planes),
    )


def assert_scalar(row, reference, model, delta, weight):
    assert row["model"] == pytest.approx(model, rel=1e-12, abs=1e-12)
    assert row["delta"] == pytest.approx(delta, rel=1e-12, abs=1e-12)
    assert row["score"] == pytest.approx(
        abs(delta) * math.sqrt(weight), rel=1e-12, abs=1e-12
    )
    assert row["residual"] == pytest.approx(reference.residual(), rel=1e-12, abs=1e-12)


def test_randomized_bond_angle_chirality_parity():
    rng = random.Random(20261011)
    for _ in range(150):
        sites = [[rng.uniform(-5.0, 5.0) for _ in range(3)] for _ in range(4)]
        weight = 10.0 ** rng.uniform(-1.0, 4.0)
        target, slack = rng.uniform(0.1, 5.0), rng.uniform(0.0, 0.5)
        proxy = restraints.bond_simple_proxy((0, 1), target, weight, slack=slack)
        reference = restraints.bond(sites_cart=flex.vec3_double(sites), proxy=proxy)
        row = graph(sites, bonds=[(0, 1, target, weight, slack, True)]).evaluate(
            sites, sigma_cutoff=0.0
        )["failures"][0]
        assert_scalar(row, reference, reference.distance_model, reference.delta, weight)
        target = rng.uniform(0.0, 180.0)
        proxy = restraints.angle_proxy((0, 1, 2), target, weight)
        reference = restraints.angle(sites_cart=flex.vec3_double(sites), proxy=proxy)
        row = graph(sites, angles=[(0, 1, 2, target, weight)]).evaluate(
            sites, sigma_cutoff=0.0
        )["failures"][0]
        assert_scalar(row, reference, reference.angle_model, reference.delta, weight)
        for both_signs in (False, True):
            target = rng.uniform(-10.0, 10.0)
            proxy = restraints.chirality_proxy((0, 1, 2, 3), target, both_signs, weight)
            reference = restraints.chirality(
                sites_cart=flex.vec3_double(sites), proxy=proxy
            )
            row = graph(
                sites, chiralities=[(0, 1, 2, 3, target, weight, both_signs, False)]
            ).evaluate(sites, sigma_cutoff=0.0)["failures"][0]
            assert_scalar(
                row, reference, reference.volume_model, reference.delta, weight
            )


def test_randomized_weighted_plane_parity():
    rng = random.Random(11)
    for _ in range(150):
        count = rng.randint(4, 15)
        sites = [[rng.uniform(-10.0, 10.0) for _ in range(3)] for _ in range(count)]
        weights = [10.0 ** rng.uniform(-1.0, 5.0) for _ in sites]
        indices = list(range(count))
        proxy = restraints.planarity_proxy(flex.size_t(indices), flex.double(weights))
        reference = restraints.planarity(
            sites_cart=flex.vec3_double(sites), proxy=proxy
        )
        row = graph(sites, planes=[(indices, weights)]).evaluate(
            sites, sigma_cutoff=0.0
        )["failures"][0]
        deltas = list(reference.deltas())
        score = max(abs(d) * math.sqrt(w) for d, w in zip(deltas, weights, strict=True))
        assert row["score"] == pytest.approx(score, rel=1e-10, abs=1e-10)
        assert row["model"] == pytest.approx(
            reference.rms_deltas(), rel=1e-10, abs=1e-10
        )
        assert row["delta"] == pytest.approx(
            max(map(abs, deltas)), rel=1e-10, abs=1e-10
        )
        assert row["residual"] == pytest.approx(
            reference.residual(), rel=1e-10, abs=1e-10
        )


def test_recorded_plane_uses_absolute_primitive_deviation():
    control = json.loads(
        (
            Path(__file__).with_name("native_validation_reference_cases")
            / "plane-signed-wrapper.json"
        ).read_text()
    )
    sites, weights = control["sites"], control["weights"]
    indices = list(range(len(sites)))
    reference = restraints.planarity(
        sites_cart=flex.vec3_double(sites),
        proxy=restraints.planarity_proxy(flex.size_t(indices), flex.double(weights)),
    )
    score = max(
        abs(d) * math.sqrt(w) for d, w in zip(reference.deltas(), weights, strict=True)
    )
    row = graph(sites, planes=[(indices, weights)]).evaluate(sites)["failures"][0]
    assert score == pytest.approx(control["score"], rel=1e-10)
    assert row["score"] == pytest.approx(score, rel=1e-10)
    assert row["residual"] == pytest.approx(reference.residual(), rel=1e-10)
    assert row["score"] > control["signed_wrapper_score"] + 1


@pytest.mark.parametrize("delta", [0.49999999, 0.5, 0.50000001])
def test_bond_decisions_at_cutoff(delta):
    sites = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
    target, weight = 1.0 + delta, 100.0
    proxy = restraints.bond_simple_proxy((0, 1), target, weight)
    reference = restraints.bond(sites_cart=flex.vec3_double(sites), proxy=proxy)
    candidate = graph(sites, bonds=[(0, 1, target, weight, 0.0, True)]).evaluate(sites)
    expected = abs(reference.delta) * math.sqrt(weight) > 5.0
    assert bool(candidate["failures"]) == expected


@pytest.mark.parametrize("offset", [-1e-7, 0.0, 1e-7])
def test_angle_and_chirality_decisions_at_cutoff(offset):
    sites = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    target = 95.0 + offset
    reference = restraints.angle(
        sites_cart=flex.vec3_double(sites),
        proxy=restraints.angle_proxy((1, 0, 2), target, 1.0),
    )
    result = graph(sites, angles=[(1, 0, 2, target, 1.0)]).evaluate(sites)
    assert bool(result["failures"]) == (abs(reference.delta) > 5.0)
    assert result["near_cutoff"] == (offset == 0.0)
    target = 1.5 + offset
    reference = restraints.chirality(
        sites_cart=flex.vec3_double(sites),
        proxy=restraints.chirality_proxy((0, 1, 2, 3), target, False, 100.0),
    )
    result = graph(
        sites, chiralities=[(0, 1, 2, 3, target, 100.0, False, False)]
    ).evaluate(sites)
    assert bool(result["failures"]) == (abs(reference.delta) * 10.0 > 5.0)
    assert result["near_cutoff"] == (offset == 0.0)


@pytest.mark.parametrize("offset", [-1e-7, 0.0, 1e-7])
def test_plane_decisions_near_cutoff(offset):
    sites = [[-1.0, -1.0, 0.0], [-1.0, 1.0, 0.0], [1.0, -1.0, 0.0], [1.0, 1.0, 0.3]]
    weights = [10000.0, 100.0, 400.0, 900.0]
    reference = restraints.planarity(
        sites_cart=flex.vec3_double(sites),
        proxy=restraints.planarity_proxy(flex.size_t(range(4)), flex.double(weights)),
    )
    score = max(
        abs(d) * math.sqrt(w) for d, w in zip(reference.deltas(), weights, strict=True)
    )
    result = graph(sites, planes=[(list(range(4)), weights)]).evaluate(
        sites, sigma_cutoff=score + offset
    )
    if offset:
        assert bool(result["failures"]) == (offset < 0.0)
    assert result["near_cutoff"] == (offset == 0.0)
