"""Jaccard estimation for sketches that support jaccard_struct."""

import random

import pytest
from weighted_cardinality_estimation.stat import (
    jaccard_streams,
    elements_stream,
    weighted_stream,
)


def test_jaccard_same_sketch_is_one(jaccard_spec) -> None:
    s = jaccard_spec.factory(64)
    s.add_many(elements_stream(20))
    assert s.jaccard_struct(s) == 1.0


def test_jaccard_same_sketch_is_one_weighted(weighted_jaccard_spec) -> None:
    s = weighted_jaccard_spec.factory(64)
    elems, weights = weighted_stream(20, 100.0)
    s.add_many(elems, weights)
    assert s.jaccard_struct(s) == 1.0


def test_jaccard_disjoint_near_zero(jaccard_spec) -> None:
    s1 = jaccard_spec.factory(64, seed=0)
    s2 = jaccard_spec.factory(64, seed=0)
    s1.add_many(elements_stream(50, seed=1))
    s2.add_many(elements_stream(50, seed=2))
    j = s1.jaccard_struct(s2)
    assert 0 <= j <= 0.05


def test_jaccard_in_valid_range(weighted_jaccard_spec) -> None:
    (ea, wa), (eb, wb) = jaccard_streams(100, 500.0, 0.5)
    s1 = weighted_jaccard_spec.factory(64, seed=0)
    s2 = weighted_jaccard_spec.factory(64, seed=0)
    s1.add_many(ea, wa)
    s2.add_many(eb, wb)
    j = s1.jaccard_struct(s2)
    assert 0 <= j <= 1.0


def test_jaccard_shuffled_stream_near_one(jaccard_spec) -> None:
    """Same elements in different order should yield Jaccard close to 1."""
    elems = elements_stream(200)
    s1 = jaccard_spec.factory(64, seed=0)
    s1.add_many(elems)

    shuffled = elems.copy()
    random.Random(42).shuffle(shuffled)
    s2 = jaccard_spec.factory(64, seed=0)
    s2.add_many(shuffled)

    assert s1.jaccard_struct(s2) > 0.9


class TestJaccardSweep:
    M = 100
    N = 500
    W = 1000.0
    JACCARDS = [0, 0.05, 0.10, 0.15, 0.20, 0.30, 0.50, 0.70, 0.80, 0.85, 0.90, 0.95, 1.0]

    @pytest.mark.parametrize("j_true", JACCARDS)
    def test_unweighted(self, jaccard_spec, j_true) -> None:
        for seed in (0, 1, 2):
            (ea, _wa), (eb, _wb) = jaccard_streams(self.N, self.W, j_true)
            s_a = jaccard_spec.factory(self.M, seed=seed)
            s_b = jaccard_spec.factory(self.M, seed=seed)
            s_a.add_many(ea)
            s_b.add_many(eb)
            est = s_a.jaccard_struct(s_b)
            assert abs(est - j_true) <= jaccard_spec.jaccard_atol

    @pytest.mark.parametrize("j_true", JACCARDS)
    def test_weighted(self, weighted_jaccard_spec, j_true) -> None:
        for seed in (0, 1, 2):
            (ea, wa), (eb, wb) = jaccard_streams(self.N, self.W, j_true)
            s_a = weighted_jaccard_spec.factory(self.M, seed=seed)
            s_b = weighted_jaccard_spec.factory(self.M, seed=seed)
            s_a.add_many(ea, wa)
            s_b.add_many(eb, wb)
            est = s_a.jaccard_struct(s_b)
            assert abs(est - j_true) <= weighted_jaccard_spec.jaccard_atol
