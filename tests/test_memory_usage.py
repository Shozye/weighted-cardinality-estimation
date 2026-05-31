import math
import pytest
from weighted_cardinality_estimation import (
    ExpSketch, ExpSketchFloat, FastExpSketch, FastGMExpSketch,
    BaseQSketch, FastQSketch, QSketch, QSketchDyn,
    BaseLogExpSketch, FastLogExpSketch,
    BaseShiftedLogExpSketch, FastShiftedLogExpSketch,
    BaseLogExpSketchJacc,
)

M = 1000
SEEDS = list(range(1, M + 1))


def compact_vector_bytes(bits: int, size: int) -> int:
    return 8 * math.ceil(size * bits / 64)


ALL_SKETCHES = [
    pytest.param(lambda: ExpSketch(M, SEEDS), id="ExpSketch"),
    pytest.param(lambda: ExpSketchFloat(M, SEEDS), id="ExpSketchFloat"),
    pytest.param(lambda: FastExpSketch(M, SEEDS), id="FastExpSketch"),
    pytest.param(lambda: FastGMExpSketch(M, SEEDS), id="FastGMExpSketch"),
    pytest.param(lambda: BaseQSketch(M, SEEDS, 8), id="BaseQSketch"),
    pytest.param(lambda: FastQSketch(M, SEEDS, amount_bits=8), id="FastQSketch"),
    pytest.param(lambda: QSketch(M, SEEDS, amount_bits=8), id="QSketch"),
    pytest.param(lambda: QSketchDyn(M, SEEDS, amount_bits=8, g_seed=42), id="QSketchDyn"),
    pytest.param(lambda: BaseLogExpSketch(M, SEEDS, amount_bits=8, logarithm_base=2), id="BaseLogExpSketch"),
    pytest.param(lambda: FastLogExpSketch(M, SEEDS, amount_bits=8, logarithm_base=2), id="FastLogExpSketch"),
    pytest.param(lambda: BaseShiftedLogExpSketch(M, SEEDS, amount_bits=8, logarithm_base=2), id="BaseShiftedLogExpSketch"),
    pytest.param(lambda: FastShiftedLogExpSketch(M, SEEDS, amount_bits=8, logarithm_base=2), id="FastShiftedLogExpSketch"),
    pytest.param(lambda: BaseLogExpSketchJacc(M, SEEDS, amount_bits=8, logarithm_base=2, amount_bits_jaccard=8), id="BaseLogExpSketchJacc"),
]


@pytest.mark.parametrize("make_sketch", ALL_SKETCHES)
def test_total_geq_write(make_sketch):
    sketch = make_sketch()
    sketch.add("elem", weight=1.0)
    assert sketch.memory_usage_total() >= sketch.memory_usage_write()


@pytest.mark.parametrize("make_sketch", ALL_SKETCHES)
def test_total_geq_estimate(make_sketch):
    sketch = make_sketch()
    sketch.add("elem", weight=1.0)
    assert sketch.memory_usage_total() >= sketch.memory_usage_estimate()


def test_exp_sketch_total():
    sketch = ExpSketch(M, SEEDS)
    # seeds [1..m] are sequential → stored as empty vector (0 bytes)
    # total = sizeof(size_t)=8 + seeds=0 + registers=8*m = 8m + 8
    assert sketch.memory_usage_total() == 8 * M + 8


def test_base_q_sketch_write():
    sketch = BaseQSketch(M, SEEDS, 8)
    # write = compact_vector(8, 1000).bytes()
    assert sketch.memory_usage_write() == compact_vector_bytes(8, M)
