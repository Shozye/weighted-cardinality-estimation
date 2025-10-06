import random
from typing import Callable, Union


from weighted_cardinality_estimation import ExpSketch, FastExpSketch, BaseQSketch, FastQSketch, QSketchDyn

SketchType = Union[ExpSketch, FastExpSketch, BaseQSketch, FastQSketch, QSketchDyn]
IMPLS: dict[str, Callable[..., SketchType]] = {
    "ExpSketch": lambda m, seeds: ExpSketch(m, seeds),
    "FastExpSketch": lambda m, seeds: FastExpSketch(m, seeds),
    "BaseQSketch(b=8)": lambda m, seeds: BaseQSketch(m, seeds, amount_bits=8),
    "FastQSketch(b=8)": lambda m, seeds: FastQSketch(m, seeds, amount_bits=8),
    "QSketchDyn(b=8)": lambda m, seeds: QSketchDyn(m, seeds, amount_bits=8, g_seed=42)
}

def get_seeds(m: int):
    return [random.randint(1, 1_000_000) for _ in range(m)]