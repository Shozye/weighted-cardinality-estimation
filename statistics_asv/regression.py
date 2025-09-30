# regression.py
# This file is used for checking if something changed between commits
# It tracks accuracy, time and memory used by sketches
from functools import cache
import random
from typing import Callable, Union
import numpy as np

from weighted_cardinality_estimation import ExpSketch, FastExpSketch, FastQSketch

SketchType = Union[ExpSketch, FastExpSketch, FastQSketch]
IMPLS: dict[str, Callable[..., SketchType]] = {
    "ExpSketch": lambda m, seeds: ExpSketch(m, seeds),
    # "FastExpSketch": lambda m, seeds: FastExpSketch(m, seeds),
    # "FastQSketch": lambda m, seeds: FastQSketch(m, seeds, amount_bits=8)
}

M_SIZE = 100
NUM_ELEMENTS = 1000
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
    cls_factory: Callable[..., SketchType]
    fresh_instance: SketchType
    filled_instance: SketchType
    elems: list[str]
    weights: list[float]
    true_cardinality: float

    param_names = ['implementation']
    params = [list(IMPLS.keys())]

    def setup_cache(self):
        return generate_data(NUM_ELEMENTS)

    def setup(self, cached_data, impl_name: str):
        self.elems, self.weights, self.true_cardinality = cached_data
        self.cls_factory = IMPLS[impl_name]
        self.fresh_instance = self.cls_factory(M_SIZE, get_seeds(M_SIZE))
        
        self.filled_instance = self.cls_factory(M_SIZE, get_seeds(M_SIZE))
        self.filled_instance.add_many(self.elems, self.weights)

    def mem_instance(self, cached_data, impl_name):
        return self.fresh_instance
    
    def time_add_many(self, cached_data, impl_name):
        self.fresh_instance.add_many(self.elems, self.weights)
        
    def time_estimate(self, cached_data, impl_name):
        self.filled_instance.estimate()

    def _get_estimates(self) -> list[float]:
        estimates = []
        for _ in range(STATISTICAL_RUNS):
            s = self.cls_factory(M_SIZE, get_seeds(M_SIZE))
            s.add_many(self.elems, self.weights)
            estimates.append(s.estimate())
        return estimates

    def track_relative_error(self, cached_data, impl_name):
        estimates = self._get_estimates()
        mean_estimate = np.mean(estimates)
        error = abs(mean_estimate - self.true_cardinality) / self.true_cardinality
        return error * 100

    def track_variance(self, cached_data, impl_name):
        estimates = self._get_estimates()
        return np.var(estimates)
    
    track_relative_error.unit = '%' # type: ignore
    track_variance.unit = 'value^2' # type: ignore 