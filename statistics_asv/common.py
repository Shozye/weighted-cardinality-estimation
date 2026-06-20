import random
from typing import Callable, Union

from weighted_cardinality_estimation import (
    ExpSketch,
    FastExpSketchFloat32,
    FastExpSketch,
    FastGMExpSketch,
    kQSketch,
    QSketch,
    QSketchDyn,
    WeightedMinHash,
)

SketchType = Union[
    ExpSketch,
    FastExpSketchFloat32,
    FastExpSketch,
    FastGMExpSketch,
    QSketchDyn,
    QSketch,
    kQSketch,
    WeightedMinHash,
]

IMPLS: dict[str, Callable[..., SketchType]] = {
    "ExpSketch": lambda m, seeds: ExpSketch(m, seeds),
    "FastExpSketchFloat32": lambda m, seeds: FastExpSketchFloat32(m, seeds),
    "FastExpSketch": lambda m, seeds: FastExpSketch(m, seeds),
    "FastGMExpSketch": lambda m, seeds: FastGMExpSketch(m, seeds),
    "WeightedMinHash": lambda m, seeds: WeightedMinHash(m, seeds),
    "QSketchDyn(b=8)": lambda m, seeds: QSketchDyn(m, seeds, amount_bits=8, g_seed=42),
    "QSketch(b=8)": lambda m, seeds: QSketch(m, seeds, amount_bits=8),
    "kQSketch(b=8,k=2)": lambda m, seeds: kQSketch(m, seeds, amount_bits=8, logarithm_base=2),
}


def get_seeds(m: int) -> list[int]:
    return [random.randint(1, 1_000_000) for _ in range(m)]
