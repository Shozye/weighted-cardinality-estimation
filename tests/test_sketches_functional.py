
import random
from tests.utils import assert_error
from weighted_cardinality_estimation import FastExpSketch, ExpSketch, FastQSketch

def functional_test(
        sketch_cls,
        err: float,
        m: int,
        amount_elements: int,
        elements_weight: float,
        am_tests: int
    ):

    total_weight = 0.0
    estimates = []
    for _ in range(am_tests):
        seeds = [random.randint(1,10000000) for _ in range(m)]
        s = sketch_cls(m, seeds)
        for i in range(amount_elements):
            elem = f"e{i}"
            total_weight += elements_weight
            s.add(elem, weight=elements_weight)
        estimates.append(s.estimate())
    
    average_estimate = sum(estimates)/len(estimates)
    print(average_estimate)
    assert_error(total_weight/am_tests, average_estimate, err)


def test_exp_sketch_functional():
    functional_test(
        sketch_cls=ExpSketch,
        err=0.05,
        m=400,
        amount_elements=1000,
        elements_weight=10,
        am_tests=100
    )

def test_fast_exp_sketch_functional():
    functional_test(
        sketch_cls=FastExpSketch,
        err=0.05,
        m=400,
        amount_elements=1000,
        elements_weight=10,
        am_tests=100
    )

def test_q_exp_sketch_functional():
    functional_test(
        sketch_cls=lambda m, seeds: FastQSketch(m, seeds, 8),
        err=0.99,
        m=10,
        amount_elements=10,
        elements_weight=10,
        am_tests=10
    )

