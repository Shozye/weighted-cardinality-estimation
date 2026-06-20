"""Experiment utilities for statistics and thesis experiments."""

from dataclasses import dataclass
from enum import Enum

import numpy as np

# ─── Weight distributions ─────────────────────────────────────────────


@dataclass(frozen=True)
class Uniform:
    low: float = 0.5
    high: float = 1.5


@dataclass(frozen=True)
class Exponential:
    scale: float = 1.0


@dataclass(frozen=True)
class Geometric:
    pass


@dataclass(frozen=True)
class Constant:
    pass


@dataclass(frozen=True)
class Pareto:
    alpha: float = 0.2


@dataclass(frozen=True)
class Monotonic:
    pass


@dataclass(frozen=True)
class Zipf:
    pass


WeightDist = Uniform | Exponential | Geometric | Constant | Pareto | Monotonic | Zipf

_DIST_ALIASES: dict[str, WeightDist] = {
    "uniform": Uniform(),
    "exponential": Exponential(),
    "geometric": Geometric(),
    "constant": Constant(),
    "pareto": Pareto(),
    "monotonic": Monotonic(),
    "zipf": Zipf(),
}


def _resolve_dist(dist: WeightDist | str) -> WeightDist:
    if isinstance(dist, str):
        if dist not in _DIST_ALIASES:
            msg = f"Unknown dist alias '{dist}'. Valid: {list(_DIST_ALIASES)}"
            raise ValueError(msg)
        return _DIST_ALIASES[dist]
    return dist


# ─── Common placement ─────────────────────────────────────────────────


class CommonPlacement(Enum):
    FIRST = "common_first"
    LAST = "common_last"
    RANDOM = "common_random"


def _resolve_common(common: "CommonPlacement | str") -> "CommonPlacement":
    if isinstance(common, str):
        return CommonPlacement(common)
    return common


# ─── Functions ─────────────────────────────────────────────────────────


def make_seeds(m: int, seed: int | None = None) -> list[int]:
    """Generate m reproducible seeds from a single integer seed."""
    return np.random.default_rng(seed).integers(1, 2**31, size=m).tolist()


def _raw_weights(dist: WeightDist, n: int, rng: np.random.Generator) -> np.ndarray:
    match dist:
        case Uniform(low, high):
            return rng.uniform(low, high, size=n)
        case Exponential(scale):
            return rng.exponential(scale, size=n)
        case Pareto(alpha):
            return rng.pareto(alpha, size=n) + 1
        case Geometric():
            return np.exp(-np.linspace(0, 12, n))
        case Monotonic():
            return np.arange(1, n + 1, dtype=float)
        case Zipf():
            return 1.0 / np.arange(1, n + 1)
        case Constant():
            return np.ones(n)


def elements_stream(n: int, seed: int | None = None) -> list[str]:
    """Generate n distinct randomized element strings."""
    rng = np.random.default_rng(seed)
    ids = rng.integers(0, 2**63, size=n)
    return [f"x{i}" for i in ids]


def weighted_stream(
    n: int,
    total_weight: float,
    dist: WeightDist = Uniform(),
    seed: int | None = None,
) -> tuple[list[str], list[float]]:
    """Generate n distinct elements with weights summing to total_weight.

    Args:
        n: Number of distinct elements.
        total_weight: Desired sum of all weights (Lambda).
        dist: Weight distribution before normalisation.
        seed: RNG seed.

    Returns:
        (elements, weights) where elements are randomized strings
        and weights sum to total_weight.

    """
    dist = _resolve_dist(dist)
    rng = np.random.default_rng(seed)
    raw = _raw_weights(dist, n, rng)
    weights = (raw / raw.sum() * total_weight).tolist()
    elements = elements_stream(n, seed)
    return elements, weights


def jaccard_streams(
    n: int,
    total_weight: float,
    jaccard: float,
    dist: WeightDist = Uniform(),
    common: CommonPlacement = CommonPlacement.RANDOM,
    seed: int | None = None,
) -> tuple[tuple[list[str], list[float]], tuple[list[str], list[float]]]:
    """Generate two weighted streams with a prescribed weighted Jaccard similarity.

    Weighted Jaccard: J_w = Σ min(w_A(x), w_B(x)) / Σ max(w_A(x), w_B(x))

    Construction:
    - Both streams have n elements, each summing to total_weight.
    - k elements are shared (identical id + weight in both streams).
    - (n-k) elements are unique to each stream.
    - J_w = sum(shared) / (2*total_weight - sum(shared))
    - Weights are drawn from `dist` and normalised per group.

    Args:
        n: Number of distinct elements per stream.
        total_weight: Desired sum of weights per stream.
        jaccard: Target weighted Jaccard similarity in [0, 1].
        dist: Weight distribution.
        common: Where shared elements appear in insertion order.
        seed: RNG seed for reproducibility.

    Returns:
        ((elems_a, weights_a), (elems_b, weights_b))

    """
    if not 0.0 <= jaccard <= 1.0:
        msg = f"jaccard must be in [0, 1], got {jaccard}"
        raise ValueError(msg)

    dist = _resolve_dist(dist)
    common = _resolve_common(common)
    rng = np.random.default_rng(seed)

    # J = S / (2W - S) => S = J * 2W / (1 + J)
    target_shared_sum = jaccard * 2 * total_weight / (1 + jaccard)
    unique_sum = total_weight - target_shared_sum

    k = round(n * target_shared_sum / total_weight) if total_weight > 0 else 0
    k = max(0, min(k, n))
    if jaccard == 1.0:
        k = n
    elif jaccard == 0.0:
        k = 0
    n_unique = n - k

    if k > 0:
        shared_w = _raw_weights(dist, k, rng)
        shared_w = shared_w / shared_w.sum() * target_shared_sum
    else:
        shared_w = np.empty(0)

    if n_unique > 0:
        a_only_w = _raw_weights(dist, n_unique, rng)
        a_only_w = a_only_w / a_only_w.sum() * unique_sum
        b_only_w = _raw_weights(dist, n_unique, rng)
        b_only_w = b_only_w / b_only_w.sum() * unique_sum
    else:
        a_only_w = np.empty(0)
        b_only_w = np.empty(0)

    shared_elems = elements_stream(k, seed)
    a_only_elems = elements_stream(n_unique, seed + 1 if seed is not None else None)
    b_only_elems = elements_stream(n_unique, seed + 2 if seed is not None else None)

    match common:
        case CommonPlacement.FIRST:
            elems_a = shared_elems + a_only_elems
            wa_out = shared_w.tolist() + a_only_w.tolist()
            elems_b = shared_elems + b_only_elems
            wb_out = shared_w.tolist() + b_only_w.tolist()
        case CommonPlacement.LAST:
            elems_a = a_only_elems + shared_elems
            wa_out = a_only_w.tolist() + shared_w.tolist()
            elems_b = b_only_elems + shared_elems
            wb_out = b_only_w.tolist() + shared_w.tolist()
        case CommonPlacement.RANDOM:
            idx_a = rng.permutation(n).tolist()
            idx_b = rng.permutation(n).tolist()
            all_elems_a = shared_elems + a_only_elems
            all_wa = shared_w.tolist() + a_only_w.tolist()
            all_elems_b = shared_elems + b_only_elems
            all_wb = shared_w.tolist() + b_only_w.tolist()
            elems_a = [all_elems_a[i] for i in idx_a]
            wa_out = [all_wa[i] for i in idx_a]
            elems_b = [all_elems_b[i] for i in idx_b]
            wb_out = [all_wb[i] for i in idx_b]

    return (elems_a, wa_out), (elems_b, wb_out)


def relative_error(estimate: float, true_value: float) -> float:
    """Compute |estimate - true_value| / true_value."""
    return abs(estimate - true_value) / true_value


def safe_relative_error(estimate: float, true_value: float, cap: float = 10.0) -> float:
    """Relative error capped at `cap` for non-finite or non-positive estimates."""
    if not np.isfinite(estimate) or estimate <= 0:
        return cap
    err = abs(estimate - true_value) / true_value
    return min(err, cap) if np.isfinite(err) else cap


def compute_rse(errors: list[float] | np.ndarray) -> float:
    """Root mean square relative error, ignoring NaN values."""
    arr = np.asarray(errors, dtype=float)
    return float(np.sqrt(np.nanmean(arr**2)))


def compact_vector_bytes(size: int, bits: int) -> int:
    """Bytes needed to store `size` values of `bits` bits each (bit-packed)."""
    return (size * bits + 7) // 8


def jaccard_variance(j: float, m: int) -> float:
    """Theoretical variance of Jaccard estimator for bottom-k / ExpSketch sketches.

    Var(Ĵ) = J(1-J) / m for all sketches based on exponential-minimum hashing.
    """
    return j * (1 - j) / m

