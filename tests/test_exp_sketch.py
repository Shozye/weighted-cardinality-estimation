from weighted_cardinality_estimation import ExpSketch
import random

def test_basic_estimate_bias_20_percent():
    m = 400
    
    estimates = []
    for trial in range(10):
        seeds = [random.randint(1,10000000) for _ in range(m)]
        s = ExpSketch(m, seeds)
        for i in range(1000):
            elem = f"e{i}"
            s.add(elem, weight=10)
        estimates.append(s.estimate())
    
    average = sum(estimates)/len(estimates)
    assert 9500 < average < 10500

