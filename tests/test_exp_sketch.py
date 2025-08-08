from weighted_cardinality_estimation import ExpSketch
import random
from tests.utils import assert_error

def test_basic_estimate_bias_20_percent():
    ERR = 0.2
    M = 400
    AMOUNT_ELEMENTS = 1000
    ELEMENT_WEIGHT = 10
    AM_TESTS = 10
    
    total_weight = 0
    estimates = []
    for _ in range(AM_TESTS):
        seeds = [random.randint(1,10000000) for _ in range(M)]
        s = ExpSketch(M, seeds)
        for i in range(AMOUNT_ELEMENTS):
            elem = f"e{i}"
            total_weight += ELEMENT_WEIGHT
            s.add(elem, weight=ELEMENT_WEIGHT)
        estimates.append(s.estimate())
    
    average_estimate = sum(estimates)/len(estimates)
    assert_error(total_weight/AM_TESTS, average_estimate, ERR)

