import random
from weighted_cardinality_estimation import FastExpSketch, ExpSketch



def unitary_test(sketch_cls):
    M=5
    
    seeds = [random.randint(1,10000000) for _ in range(M)]
    sketch = sketch_cls(M, seeds)
    sketch.add("I am just a simple element.", weight=1)

    estimate = sketch.estimate()
    assert estimate > 0


def test_exp_sketch_functional():
    unitary_test(
        sketch_cls=ExpSketch,
    )

def test_fast_exp_sketch_functional():
    unitary_test(
        sketch_cls=FastExpSketch,
    )


