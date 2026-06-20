"""Viable weight range validation for each sketch family.

Each sketch has a representable range of weights [min_weight, max_weight].
Tests use m=400 (expected RSE ~5%) and check:
1. Single element at min_weight → estimate within 20% of min_weight
2. Single element at max_weight → estimate within 20% of max_weight
3. Sweep from min to max in powers of 10 — each step within 20% of weight
"""

import math
from collections.abc import Callable
from dataclasses import dataclass

import pytest
from weighted_cardinality_estimation import (

    ExpSketch,
    ExpSketchFloat32,
    LogExpSketchSlowNoShifted,
    FastExpSketch,
    FastExpSketchCustomFloat,
    FastExpSketchFloat32,
    LogExpSketchFastNoShifted,
    FastGMExpSketch,

    QSketch,
    QSketchDyn,
    WeightedHyperLogLog,
    WeightedHyperLogLogCustomFloat,
    WeightedHyperLogLogFloat32,
    WeightedMinHash,
    kQSketch,
    kQSketchRounding,
    kQSketchShifted,
)

from custom_float_configs import ALL_CONFIGS
M = 400
REL_TOL = 0.20


@dataclass
class RangeSpec:
    name: str
    factory: Callable[[], object]
    min_weight: float
    max_weight: float
    estimator: str = "estimate"
    n_elems: int = 1


RANGE_SPECS: list[RangeSpec] = [
    # --- ExpSketch float64 family ---
    RangeSpec("ExpSketch", lambda: ExpSketch(M, seed=0), 1e-305, 1e308),
    RangeSpec("FastExpSketch", lambda: FastExpSketch(M, seed=0), 1e-305, 1e308),
    RangeSpec("FastGMExpSketch", lambda: FastGMExpSketch(M, seed=0), 1e-305, 1e308),
    RangeSpec("WeightedMinHash", lambda: WeightedMinHash(M, seed=0), 1e1, 1e16),
    RangeSpec("WeightedHyperLogLog", lambda: WeightedHyperLogLog(M, seed=0), 1e-305, 1e304,
              n_elems=5000),
    RangeSpec("WeightedHyperLogLogFloat32",
              lambda: WeightedHyperLogLogFloat32(M, seed=0), 1e-37, 1e38, n_elems=5000),
    # --- ExpSketch float32 family ---
    RangeSpec("ExpSketchFloat32", lambda: ExpSketchFloat32(M, seed=0), 1e-37, 1e45),
    RangeSpec("FastExpSketchFloat32", lambda: FastExpSketchFloat32(M, seed=0), 1e-37, 1e45),
    # --- QSketch 8-bit family ---

    RangeSpec("QSketch_b8", lambda: QSketch(M, seed=0, amount_bits=8), 1e-37, 1e38),
    RangeSpec("QSketchDyn_b8", lambda: QSketchDyn(M, seed=0, amount_bits=8, g_seed=42), 1e-37, 1e308),
    # --- kQSketch 8-bit base-2 ---
    RangeSpec("kQSketch_b8_k2", lambda: kQSketch(M, seed=0, amount_bits=8, logarithm_base=2.0), 1e-37, 1e38),
    RangeSpec("kQSketchRounding_b8_k2", lambda: kQSketchRounding(M, seed=0, amount_bits=8, logarithm_base=2.0), 1e-37, 1e38),
    # --- kQSketchShifted variants ---
    RangeSpec("kQSketchShifted_b4_k2", lambda: kQSketchShifted(M, seed=0, amount_bits=4, logarithm_base=2.0), 1e-305, 1e307,
              estimator="estimate_direct"),
    RangeSpec("kQSketchShifted_b5_k2", lambda: kQSketchShifted(M, seed=0, amount_bits=5, logarithm_base=2.0), 1e-305, 1e308,
              estimator="estimate_direct"),
    RangeSpec("kQSketchShifted_b6_k2", lambda: kQSketchShifted(M, seed=0, amount_bits=6, logarithm_base=2.0), 1e-305, 1e308,
              estimator="estimate_direct"),
    RangeSpec("kQSketchShifted_b8_k2", lambda: kQSketchShifted(M, seed=0, amount_bits=8, logarithm_base=2.0), 1e-305, 1e308,
              estimator="estimate_direct"),
    RangeSpec(
        "kQSketchShifted_b10_k2", lambda: kQSketchShifted(M, seed=0, amount_bits=10, logarithm_base=2.0), 1e-305, 1e308,
        estimator="estimate_direct",
    ),
    RangeSpec(
        "kQSketchShifted_b12_k2", lambda: kQSketchShifted(M, seed=0, amount_bits=12, logarithm_base=2.0), 1e-305, 1e308,
        estimator="estimate_direct",
    ),
    # --- FastExpSketchCustomFloat (5+10, ALL_NORMAL) ---
    # --- LogExpSketchSlowNoShifted (10 bits, v_max=1e5) ---
    RangeSpec(
        "LogExpSketchSlowNoShifted_b10_vmax1e5", lambda: LogExpSketchSlowNoShifted(M, seed=0, amount_bits=10, v_max=1e5),
        1e-4, 1e4,
    ),
    RangeSpec(
        "LogExpSketchFastNoShifted_b10_vmax1e5", lambda: LogExpSketchFastNoShifted(M, seed=0, amount_bits=10, v_max=1e5),
        1e-4, 1e4,
    ),
    # --- Custom float standard configs (all 8) ---
    *[
        RangeSpec(
            f"FastExpSketchCustomFloat_{cfg.label}",
            lambda c=cfg: FastExpSketchCustomFloat(M, seed=0, exp_bits=c.exp_bits, mant_bits=c.mant_bits),
            cfg.min_weight, cfg.max_weight,
        )
        for cfg in ALL_CONFIGS
    ],
    *[
        RangeSpec(
            f"WeightedHyperLogLogCustomFloat_{cfg.label}",
            lambda c=cfg: WeightedHyperLogLogCustomFloat(M, seed=0, exp_bits=c.exp_bits, mant_bits=c.mant_bits),
            cfg.min_weight,
            # WHLL uses n_elems=5000, so true cardinality = 5000*weight.
            # Must keep 5000*max_weight < float64_max (~1.7e308).
            1e36 if cfg.exp_bits == 8 else 1e304,
            n_elems=5000,
        )
        for cfg in ALL_CONFIGS
    ],
]


@pytest.fixture(params=RANGE_SPECS, ids=lambda s: s.name)
def range_spec(request):
    return request.param


def test_min_weight(range_spec) -> None:
    sketch = range_spec.factory()
    for i in range(range_spec.n_elems):
        sketch.add(f"e{i}", range_spec.min_weight)
    expected = range_spec.min_weight * range_spec.n_elems
    est = getattr(sketch, range_spec.estimator)()
    assert math.isfinite(est) and est > 0, f"estimate={est}"
    assert abs(est - expected) / expected <= REL_TOL, (
        f"estimate={est}, expected≈{expected}"
    )


def test_max_weight(range_spec) -> None:
    sketch = range_spec.factory()
    for i in range(range_spec.n_elems):
        sketch.add(f"e{i}", range_spec.max_weight)
    expected = range_spec.max_weight * range_spec.n_elems
    est = getattr(sketch, range_spec.estimator)()
    assert math.isfinite(est) and est > 0, f"estimate={est}"
    assert abs(est - expected) / expected <= REL_TOL, (
        f"estimate={est}, expected≈{expected}"
    )


def test_log_sweep(range_spec) -> None:
    log_min = math.ceil(math.log10(range_spec.min_weight))
    log_max = math.floor(math.log10(range_spec.max_weight))
    failures = []
    for exp in range(log_min, log_max + 1):
        weight = 10.0**exp
        expected = weight * range_spec.n_elems
        sketch = range_spec.factory()
        for i in range(range_spec.n_elems):
            sketch.add(f"e{i}", weight)
        est = getattr(sketch, range_spec.estimator)()
        if not (math.isfinite(est) and est > 0):
            failures.append(f"w=1e{exp}: estimate={est} not finite/positive")
        elif abs(est - expected) / expected > REL_TOL:
            failures.append(
                f"w=1e{exp}: estimate={est}, expected≈{expected}, "
                f"rel_error={abs(est - expected) / expected:.4f}"
            )
    assert not failures, "Failing exponents:\n" + "\n".join(failures)
