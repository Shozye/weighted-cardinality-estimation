from weighted_cardinality_estimation import ExpSketch
import random

def test_add_speed(benchmark):
    seeds = list(range(400))
    sketch = ExpSketch(400, seeds)

    @benchmark
    def _():
        x = random.random()
        sketch.add(str(x), 1.0)
