"""Benchmarks: update throughput and estimate time for all sketches."""

import pytest
from conftest import SKETCH_SPECS, WEIGHTED_SPECS
from weighted_cardinality_estimation.stat import elements_stream, weighted_stream

pytestmark = pytest.mark.benchmark

M = 128
N = 500
ELEMS = elements_stream(N)
ELEMS_W, WEIGHTS = weighted_stream(N, 1e4, seed=0)

WEIGHTED_NAMES = {s.name for s in WEIGHTED_SPECS}


@pytest.mark.parametrize("spec", SKETCH_SPECS, ids=lambda s: s.name)
def test_bench_update(benchmark, spec) -> None:
    """Benchmark: add_many (weighted if supported, else unweighted)."""
    if spec.name in WEIGHTED_NAMES:
        benchmark.pedantic(
            lambda s: s.add_many(ELEMS_W, WEIGHTS),
            setup=lambda: ((spec.factory(M),), {}),
            rounds=100,
        )
    else:
        benchmark.pedantic(
            lambda s: s.add_many(ELEMS),
            setup=lambda: ((spec.factory(M),), {}),
            rounds=100,
        )


@pytest.mark.parametrize("spec", SKETCH_SPECS, ids=lambda s: s.name)
def test_bench_estimate(benchmark, spec) -> None:
    """Benchmark: estimate() on a populated sketch."""
    s = spec.factory(M)
    if spec.name in WEIGHTED_NAMES:
        s.add_many(ELEMS_W, WEIGHTS)
    else:
        s.add_many(ELEMS)
    benchmark.pedantic(s.estimate, rounds=100)
