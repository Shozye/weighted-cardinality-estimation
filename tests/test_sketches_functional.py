import random
import pytest
from tests.utils import assert_error
from weighted_cardinality_estimation import FastExpSketch, ExpSketch, FastQSketch

M_SIZE = 400
AMOUNT_ELEMENTS = 1000
ELEMENTS_WEIGHT = 10.0
AMOUNT_TEST_RUNS = 100


SKETCH_PARAMS = [
    pytest.param(ExpSketch, 0.05, id="ExpSketch"),
    pytest.param(FastExpSketch, 0.05, id="FastExpSketch"),
    # Dla FastQSketch używamy lambdy, aby przekazać dodatkowy argument `amount_bits`.
    pytest.param(lambda m, seeds: FastQSketch(m, seeds, 8), 0.1, id="FastQSketch"),
]

@pytest.mark.parametrize("sketch_cls, allowed_error", SKETCH_PARAMS)
def test_sketch_functional_accuracy(
        sketch_cls,
        allowed_error: float,
    ):
    # here I want to assert some level of error on the structures to make sure they are any good

    total_weight = AMOUNT_ELEMENTS * ELEMENTS_WEIGHT
    estimates = []
    for _ in range(AMOUNT_TEST_RUNS):
        seeds = [random.randint(1,10000000) for _ in range(M_SIZE)]
        s = sketch_cls(M_SIZE, seeds)
        elements = [f"e{i}" for i in range(AMOUNT_ELEMENTS)]
        weights = [ELEMENTS_WEIGHT] * AMOUNT_ELEMENTS
        s.add_many(elements, weights)
        estimates.append(s.estimate())
    
    average_estimate = sum(estimates)/len(estimates)
    print(average_estimate)
    assert_error(total_weight, average_estimate, allowed_error)


