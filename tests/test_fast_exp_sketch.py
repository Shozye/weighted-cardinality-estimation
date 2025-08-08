from weighted_cardinality_estimation import FastExpSketch
import random
from tests.utils import assert_error

def test_basic_estimate_bias_20_percent():
    ERR = 0.5
    M = 5
    AMOUNT_ELEMENTS = 1
    ELEMENT_WEIGHT = 1
    AM_TESTS=1
    
    total_weight = 0
    estimates = []
    for _ in range(AM_TESTS):
        seeds = [random.randint(1,10000000) for _ in range(M)]
        s = FastExpSketch(M, seeds)
        for i in range(AMOUNT_ELEMENTS):
            elem = f"e{i}"
            total_weight += ELEMENT_WEIGHT
            s.add(elem, weight=ELEMENT_WEIGHT)
        estimates.append(s.estimate())
    
    average_estimate = sum(estimates)/len(estimates)
    assert_error(total_weight/1, average_estimate, ERR)

