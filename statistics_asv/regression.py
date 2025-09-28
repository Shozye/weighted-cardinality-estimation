# regression.py
# This file is used for checking if something changed between commits
# It tracks accuracy, time and memory used by sketches
import random
from typing import Literal
import numpy as np

from weighted_cardinality_estimation import ExpSketch, FastExpSketch, FastQSketch

M_SIZE = 512
NUM_ELEMENTS = 100_000
STATISTICAL_RUNS = 50

def generate_data(n: int) -> tuple[list[str], list[float], float]:
    random.seed(n)
    elems = [f"element_{i}" for i in range(n)]
    weights = [random.uniform(1, 10.0) for _ in range(n)]
    true_cardinality = sum(weights)
    return elems, weights, true_cardinality

def get_seeds(m: int):
    return [random.randint(1, 1_000_000) for _ in range(m)]

class RegressionSuite:
    param_names = ['implementation']
    params = [['ExpSketch', 'FastExpSketch', 'FastQSketch']]
    
    IMPLS = {
        "ExpSketch": lambda m, seeds: ExpSketch(m, seeds),
        "FastExpSketch": lambda m, seeds: FastExpSketch(m, seeds),
        "FastQSketch": lambda m, seeds: FastQSketch(m, seeds, amount_bits=8)
    }

    def setup_cache(self):
        self.elems, self.weights, self.true_cardinality = generate_data(NUM_ELEMENTS)

    def setup(self, impl_name):
        self.cls_factory = self.IMPLS[impl_name]
        self.fresh_instance = self.cls_factory(M_SIZE, get_seeds(M_SIZE))
        
        self.filled_instance = self.cls_factory(M_SIZE, get_seeds(M_SIZE))
        self.filled_instance.add_many(self.elems, self.weights)


    def mem_instance(self, _):
        return self.fresh_instance

    def time_add_many(self, _):
        self.fresh_instance.add_many(self.elems, self.weights)
        
    def time_estimate(self, _):
        self.filled_instance.estimate()

    def _get_estimates(self) -> list[float]:
        estimates = []
        for _ in range(STATISTICAL_RUNS):
            s = self.cls_factory(M_SIZE, get_seeds(M_SIZE))
            s.add_many(self.elems, self.weights)
            estimates.append(s.estimate())
        return estimates

    def track_relative_error(self, _):
        estimates = self._get_estimates()
        mean_estimate = np.mean(estimates)
        error = abs(mean_estimate - self.true_cardinality) / self.true_cardinality
        return error * 100

    def track_variance(self, _):
        estimates = self._get_estimates()
        return np.var(estimates)
    
    track_relative_error.unit = '%' # type: ignore
    track_variance.unit = 'value^2' # type: ignore 