import random
from typing import Callable, Union
import numpy as np

from weighted_cardinality_estimation import ExpSketch, FastExpSketch, FastQSketch

SketchType = Union[ExpSketch, FastExpSketch, FastQSketch]
IMPLS: dict[str, Callable[..., SketchType]] = {
    "ExpSketch": lambda m, seeds: ExpSketch(m, seeds),
    "FastExpSketch": lambda m, seeds: FastExpSketch(m, seeds),
    "FastQSketch(b=8)": lambda m, seeds: FastQSketch(m, seeds, amount_bits=8)
}

def get_seeds(m: int):
    return [random.randint(1, 1_000_000) for _ in range(m)]