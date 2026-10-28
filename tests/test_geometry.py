"""Tests for explicit restraint evaluation, independent of chemistry assignment."""

import math

import pytest

from nibbler.geometry import RestraintGraph


def test_reusable_graph_preserves_atom_ids_and_strict_cutoff():
    graph = RestraintGraph(
        ["C", "H"], [False, True], [(0, 1, 1.5, 100.0, 0.4, True)], [], [], []
    )
    at = graph.evaluate([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    assert at["failures"] == []
    assert at["near_cutoff"] is True
    below = graph.evaluate([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], sigma_cutoff=4.999)
    assert below["checks"]["bonds"] == {"status": "checked", "count": 1}
    row = below["failures"][0]
    assert row["atom_ids"] == ["C", "H"]
    assert row["proxy_index"] == 0
    assert row["score"] == 5.0
    assert row["residual"] == pytest.approx(1.0)
    assert below["near_cutoff"] is False


@pytest.mark.parametrize("offset", [-1e-9, 1e-9])
def test_near_cutoff_covers_both_sides(offset):
    graph = RestraintGraph(
        ["C", "H"], [False, True], [(0, 1, 1.5 + offset, 100.0, 0.0, True)], [], [], []
    )
    result = graph.evaluate([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    assert result["near_cutoff"] is True
    assert bool(result["failures"]) == (offset > 0.0)


@pytest.mark.parametrize("sites", [[[math.nan, 0.0, 0.0]], [], [[math.inf, 0.0, 0.0]]])
def test_bad_coordinates_fail(sites):
    graph = RestraintGraph(["C"], [False], [], [], [], [])
    with pytest.raises(ValueError, match="Coordinates"):
        graph.evaluate(sites)


@pytest.mark.parametrize("weight", [0.0, -1.0, math.nan, math.inf])
def test_bad_weights_fail(weight):
    with pytest.raises(ValueError, match="weights"):
        RestraintGraph(
            ["C", "H"], [False, True], [(0, 1, 1.0, weight, 0.0, True)], [], [], []
        )


def test_missing_hydrogen_parent_is_a_finding():
    graph = RestraintGraph(["H"], [True], [], [], [], [])
    result = graph.evaluate([[0.0, 0.0, 0.0]])
    assert result["failures"] == [
        {
            "check": "hydrogens",
            "atom": "H",
            "indices": [0],
            "parents": [],
            "reason": "hydrogen_parent_count",
        }
    ]
