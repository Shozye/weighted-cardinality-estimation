from weighted_cardinality_estimation import FastExpSketch
import random


def test_basic_estimate_bias_percent():
    m = 10
    err = 0.5
    
    seeds = [random.randint(1,10000000) for _ in range(m)]
    s = FastExpSketch(m, seeds)

    total_sum = 0
    for i in range(1):
        elem = f"e{i}"
        w = 10
        s.add(elem, weight=w)
        total_sum += w
    

    
    esimate = s.estimate()
    assert (total_sum*(1-err)) < esimate < (total_sum*(1+err))



# def test_basic_estimate_bias_20_percent():
#     m = 400
    
#     estimates = []
#     for trial in range(10):
#         seeds = [random.randint(1,10000000) for _ in range(m)]
#         s = FastExpSketch(m, seeds)
#         for i in range(1000):
#             elem = f"e{i}"
#             s.add(elem, weight=10)
#         estimates.append(s.estimate())
    
#     average = sum(estimates)/len(estimates)
#     assert 9500 < average < 10500

