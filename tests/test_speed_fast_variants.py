"""Speed tests: fast sketch variants must be >= 8x faster than ExpSketch after warmup."""

import time
from statistics import median

import pytest
from conftest import FAST_SPECS
from weighted_cardinality_estimation import ExpSketch
from weighted_cardinality_estimation.stat import elements_stream, weighted_stream

pytestmark = pytest.mark.benchmark

M = 128
REPS = 20
ELEMS = elements_stream(2000)
WARMUP_ELEMS = elements_stream(10_000)
ELEMS_W, WEIGHTS = weighted_stream(2000, 1e4, seed=0)
MIN_SPEEDUP = 8


def _median_add_time(make_sketch, elems, weights=None):
    """Median time to add_many on a fresh sketch from make_sketch()."""
    def run():
        s = make_sketch()
        t0 = time.perf_counter()
        if weights is None:
            s.add_many(elems)
        else:
            s.add_many(elems, weights)
        return time.perf_counter() - t0

    return median(run() for _ in range(REPS))


def _expsketch_warmed():
    s = ExpSketch(M, seed=42)
    s.add("warmup", 1e4)
    return s


class TestSpeedupWarmed:
    """All fast variants must beat ExpSketch by at least MIN_SPEEDUP after warmup."""

    @pytest.mark.parametrize("spec", FAST_SPECS, ids=lambda s: s.name)
    def test_unweighted(self, spec) -> None:
        base_time = _median_add_time(_expsketch_warmed, ELEMS)

        def make():
            s = spec.factory(M)
            s.add_many(WARMUP_ELEMS)
            return s

        speedup = base_time / _median_add_time(make, ELEMS)
        assert speedup >= MIN_SPEEDUP, f"{spec.name}: {speedup:.1f}x (need >={MIN_SPEEDUP}x)"

    def test_weighted(self, fast_weighted_spec) -> None:
        base_time = _median_add_time(_expsketch_warmed, ELEMS_W, WEIGHTS)

        def make():
            s = fast_weighted_spec.factory(M)
            s.add("warmup", 1e4)
            return s

        speedup = base_time / _median_add_time(make, ELEMS_W, WEIGHTS)
        assert speedup >= MIN_SPEEDUP, f"{fast_weighted_spec.name}: {speedup:.1f}x (need >={MIN_SPEEDUP}x)"
