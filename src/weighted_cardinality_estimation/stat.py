"""Experiment utilities for statistics and thesis experiments."""
from typing import Literal
import numpy as np

WeightDist = Literal["uniform", "exponential", "geometric", "constant"]


def make_seeds(m: int, seed: int = 42) -> list[int]:
    """Generate m reproducible seeds from a single integer seed."""
    return np.random.default_rng(seed).integers(1, 2**31, size=m).tolist()


def stream(
    n: int,
    total_weight: float,
    dist: WeightDist = "uniform",
    seed: int | None = None,
    low: float = 0.5,
    high: float = 1.5,
) -> tuple[list[str], list[float]]:
    """Generate n distinct elements with weights summing to total_weight.

    Args:
        n: Number of distinct elements.
        total_weight: Desired sum of all weights (Lambda).
        dist: Weight distribution before normalisation.
            - "uniform":     weights ~ Uniform(low, high)  [default]
            - "exponential": weights ~ Exponential(1)
            - "geometric":   weights geometrically decreasing (deterministic)
            - "constant":    all weights equal total_weight / n
        seed: RNG seed (only used for stochastic distributions).
        low, high: Bounds for the uniform distribution.

    Returns:
        (elements, weights) where elements are strings "x0", "x1", ...
        and weights sum to total_weight.
    """
    rng = np.random.default_rng(seed)
    if dist == "uniform":
        raw = rng.uniform(low, high, size=n)
    elif dist == "exponential":
        raw = rng.exponential(1.0, size=n)
    elif dist == "geometric":
        raw = np.exp(-np.linspace(0, 4, n))
    elif dist == "constant":
        raw = np.ones(n)
    else:
        raise ValueError(f"Unknown dist {dist!r}. Choose from: uniform, exponential, geometric, constant")

    weights = (raw / raw.sum() * total_weight).tolist()
    elements = [f"x{i}" for i in range(n)]
    return elements, weights


def relative_error(estimate: float, true_value: float) -> float:
    """Compute |estimate - true_value| / true_value."""
    return abs(estimate - true_value) / true_value


def memory_breakdown(sketch) -> dict[str, int]:
    """Return total/write/estimate memory usage in bytes as a dict."""
    return {
        "total": sketch.memory_usage_total(),
        "write": sketch.memory_usage_write(),
        "estimate": sketch.memory_usage_estimate(),
    }
