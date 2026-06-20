"""Declarative sketch registry and shared fixtures."""

from collections.abc import Callable
from dataclasses import dataclass, field

import pytest
from weighted_cardinality_estimation import (

    CardinalitySketch,
    ExpSketch,
    ExpSketchFloat32,
    LogExpSketchSlowNoShifted,
    LogExpSketchSlowShifted,
    FastExpSketch,
    FastExpSketchCustomFloat,
    FastExpSketchFloat32,
    LogExpSketchFastNoShifted,
    LogExpSketchFastShifted,
    FastGMExpSketch,

    HyperLogLog,
    MartingaleMinHash,
    MemoryFlag,
    MinHash,
    QSketch,
    QSketchDyn,
    QuantizationMode,
    WeightedHyperLogLog,
    WeightedHyperLogLogCustomFloat,
    WeightedHyperLogLogFloat32,
    WeightedMinHash,
    kQSketch,
    kQSketchRounding,
    kQSketchRoundedDyn,
    kQSketchShifted,
)
from weighted_cardinality_estimation.stat import elements_stream

from custom_float_configs import SKETCH_SPEC_CONFIGS

M = 64


def _add_elements_stream(sketch, n: int) -> None:
    """Add n distinct elements to a sketch."""
    sketch.add_many(elements_stream(n))


def _probe_weighted(instance) -> bool:
    try:
        instance.add("__probe__", 1.0)
        return True
    except TypeError:
        return False


@dataclass
class SketchSpec:
    name: str
    factory: Callable[[int, list[int] | None], CardinalitySketch]
    min_weight: float
    max_weight: float
    is_fast: bool = False
    jaccard_atol: float = 0.25
    # Maximum acceptable relative error |estimate - true| / true for cardinality tests.
    # Tuned per sketch for m=400, n=500, over 20 random seeds. Set to observed max + ~2% margin.
    estimate_rel_error: float = 0.5
    # Auto-detected from factory instance — not settable.
    has_jaccard: bool = field(init=False)
    has_newton: bool = field(init=False)
    has_merge: bool = field(init=False)
    is_weighted: bool = field(init=False)

    def __post_init__(self):
        instance = self.factory(100)
        self.has_jaccard = hasattr(instance, "jaccard_struct")
        self.has_merge = hasattr(instance, "merge")
        self.has_newton = hasattr(instance, "estimate_newton_cold")
        self.is_weighted = _probe_weighted(instance)
        # FisherYates stores two permutation arrays → extra write memory
        cls_name = type(instance).__name__.lower()
        self.is_fast = self.is_fast or "fast" in cls_name


SKETCH_SPECS: list[SketchSpec] = [
    SketchSpec("ExpSketch", lambda m, seed=42: ExpSketch(m, seed=seed),
               estimate_rel_error=0.12, min_weight=1e-305, max_weight=1e307),
    SketchSpec("ExpSketchFloat32", lambda m, seed=42: ExpSketchFloat32(m, seed=seed),
               estimate_rel_error=0.10, min_weight=1e-37, max_weight=1e45),
    SketchSpec(
        "FastExpSketchFloat32", lambda m, seed=42: FastExpSketchFloat32(m, seed=seed),
        estimate_rel_error=0.10, min_weight=1e-37, max_weight=1e45,
    ),
    SketchSpec(
        "FastExpSketch", lambda m, seed=42: FastExpSketch(m, seed=seed),
        estimate_rel_error=0.12, min_weight=1e-305, max_weight=1e307,
    ),
    SketchSpec(
        "FastGMExpSketch", lambda m, seed=42: FastGMExpSketch(m, seed=seed),
        estimate_rel_error=0.12, min_weight=1e-305, max_weight=1e307,
    ),

    SketchSpec("QSketch", lambda m, seed=42: QSketch(m, seed=seed, amount_bits=8),
               estimate_rel_error=0.14, min_weight=1e-37, max_weight=1e38),
    SketchSpec(
        "QSketchDyn", lambda m, seed=42: QSketchDyn(m, seed=seed, amount_bits=8, g_seed=42),
        is_fast=True, estimate_rel_error=0.07, min_weight=1e-37, max_weight=1e307,
    ),
    SketchSpec(
        "kQSketch", lambda m, seed=42: kQSketch(m, seed=seed, amount_bits=8, logarithm_base=2),
        estimate_rel_error=0.12, min_weight=1e-37, max_weight=1e38,
    ),
    SketchSpec(
        "kQSketchRounding",
        lambda m, seed=42: kQSketchRounding(m, seed=seed, amount_bits=8, logarithm_base=2),
        estimate_rel_error=0.11, min_weight=1e-37, max_weight=1e38,
    ),
    SketchSpec(
        "kQSketchRoundedDyn",
        lambda m, seed=42: kQSketchRoundedDyn(m, seed=seed, amount_bits=8, logarithm_base=2, g_seed=42),
        is_fast=True, estimate_rel_error=0.07, min_weight=1e-37, max_weight=1e307,
    ),
    SketchSpec(
        "kQSketchShifted",
        lambda m, seed=42: kQSketchShifted(m, seed=seed, amount_bits=8, logarithm_base=2),
        estimate_rel_error=0.11, min_weight=1e-305, max_weight=1e307,
    ),
    SketchSpec("WeightedMinHash", lambda m, seed=42: WeightedMinHash(m, seed=seed),
               estimate_rel_error=0.11, min_weight=1e1, max_weight=1e16),
    SketchSpec("WeightedHyperLogLog", lambda m, seed=42: WeightedHyperLogLog(m, seed=seed),
               is_fast=True, estimate_rel_error=0.59, min_weight=1e-305, max_weight=1e307),
    SketchSpec("WeightedHyperLogLogFloat32",
               lambda m, seed=42: WeightedHyperLogLogFloat32(m, seed=seed),
               is_fast=True, estimate_rel_error=0.59, min_weight=1e-37, max_weight=1e38),
    SketchSpec("MinHash", lambda m, seed=42: MinHash(m, seed=seed),
               min_weight=0.0, max_weight=0.0, estimate_rel_error=0.12),
    SketchSpec("HyperLogLog", lambda m, seed=42: HyperLogLog(m, seed=seed),
               min_weight=0.0, max_weight=0.0, is_fast=True, estimate_rel_error=0.09),
    SketchSpec("MartingaleMinHash", lambda m, seed=42: MartingaleMinHash(m, seed=seed),
               min_weight=0.0, max_weight=0.0, estimate_rel_error=0.06),
    # Custom float sketches — parametrized over standard configs
    *[
        SketchSpec(
            f"FastExpSketchCustomFloat_{cfg.label}",
            lambda m, seed=42, c=cfg: FastExpSketchCustomFloat(
                m, seed=seed, exp_bits=c.exp_bits, mant_bits=c.mant_bits,
            ),
            estimate_rel_error=0.12, min_weight=cfg.min_weight, max_weight=cfg.max_weight,
        )
        for cfg in SKETCH_SPEC_CONFIGS
    ],
    *[
        SketchSpec(
            f"WeightedHyperLogLogCustomFloat_{cfg.label}",
            lambda m, seed=42, c=cfg: WeightedHyperLogLogCustomFloat(
                m, seed=seed, exp_bits=c.exp_bits, mant_bits=c.mant_bits,
            ),
            is_fast=True,
            estimate_rel_error=0.60, min_weight=cfg.min_weight, max_weight=cfg.max_weight,
        )
        for cfg in SKETCH_SPEC_CONFIGS
    ],
    SketchSpec(
        "LogExpSketchSlowNoShifted",
        lambda m, seed=42: LogExpSketchSlowNoShifted(m, seed=seed, amount_bits=10, v_max=1e5),
        estimate_rel_error=0.15, min_weight=1e-4, max_weight=1e4,
    ),
    SketchSpec(
        "LogExpSketchSlowShifted",
        lambda m, seed=42: LogExpSketchSlowShifted(m, seed=seed, amount_bits=10, v_max=1e5),
        estimate_rel_error=0.15, min_weight=1e-4, max_weight=1e307,
    ),
    SketchSpec(
        "LogExpSketchFastNoShifted",
        lambda m, seed=42: LogExpSketchFastNoShifted(m, seed=seed, amount_bits=10, v_max=1e5),
        estimate_rel_error=0.11, min_weight=1e-4, max_weight=1e4,
    ),
    SketchSpec(
        "LogExpSketchFastShifted",
        lambda m, seed=42: LogExpSketchFastShifted(m, seed=seed, amount_bits=10, v_max=1e5),
        estimate_rel_error=0.15, min_weight=1e-4, max_weight=1e307,
    ),
]

JACCARD_SPECS = [s for s in SKETCH_SPECS if s.has_jaccard]
WEIGHTED_JACCARD_SPECS = [s for s in SKETCH_SPECS if s.has_jaccard and s.is_weighted]
NEWTON_SPECS = [s for s in SKETCH_SPECS if s.has_newton]
MERGE_SPECS = [s for s in SKETCH_SPECS if s.has_merge]
WEIGHTED_SPECS = [s for s in SKETCH_SPECS if s.is_weighted]
FAST_SPECS = [s for s in SKETCH_SPECS if s.is_fast]
FAST_WEIGHTED_SPECS = [s for s in SKETCH_SPECS if s.is_fast and s.is_weighted]


@pytest.fixture(params=SKETCH_SPECS, ids=lambda s: s.name)
def spec(request):
    """Parametrized over all sketch types."""
    return request.param


@pytest.fixture
def sketch(spec):
    """A fresh sketch with standard M and sequential seeds."""
    return spec.factory(M)


@pytest.fixture(params=JACCARD_SPECS, ids=lambda s: s.name)
def jaccard_spec(request):
    return request.param


@pytest.fixture(params=WEIGHTED_JACCARD_SPECS, ids=lambda s: s.name)
def weighted_jaccard_spec(request):
    return request.param


@pytest.fixture(params=NEWTON_SPECS, ids=lambda s: s.name)
def newton_spec(request):
    return request.param


@pytest.fixture(params=MERGE_SPECS, ids=lambda s: s.name)
def merge_spec(request):
    return request.param


@pytest.fixture(params=WEIGHTED_SPECS, ids=lambda s: s.name)
def weighted_spec(request):
    return request.param


@pytest.fixture(params=FAST_SPECS, ids=lambda s: s.name)
def fast_spec(request):
    return request.param


@pytest.fixture(params=[
    MemoryFlag.REGISTERS,
    MemoryFlag.SEEDS,
    MemoryFlag.FISHER_YATES_PERM_INIT,
    MemoryFlag.FISHER_YATES_NON_PERM_INIT,
], ids=["REGISTERS", "SEEDS", "FY_PERM", "FY_NON_PERM"])
def memory_flag(request):
    return request.param


WEIGHTED_MERGE_SPECS = [s for s in MERGE_SPECS if s.is_weighted]
ALL_QUANTIZATION_MODES = [
    QuantizationMode.ALL_NORMAL,
    QuantizationMode.WITH_SUBNORMALS,
    QuantizationMode.LINEAR,
    QuantizationMode.LOGARITHMIC,
]
IEEE_QUANTIZATION_MODES = [QuantizationMode.ALL_NORMAL, QuantizationMode.WITH_SUBNORMALS]


def make_sketches(spec, m: int, seed: int, n: int = 2):
    """Create n fresh sketches from spec with the same seeds."""
    return [spec.factory(m, seed) for _ in range(n)]


@pytest.fixture(params=WEIGHTED_MERGE_SPECS, ids=lambda s: s.name)
def weighted_merge_spec(request):
    return request.param


@pytest.fixture(params=FAST_WEIGHTED_SPECS, ids=lambda s: s.name)
def fast_weighted_spec(request):
    return request.param


@pytest.fixture(params=ALL_QUANTIZATION_MODES, ids=lambda m: str(m).split(".")[-1])
def quantization_mode(request):
    return request.param


@pytest.fixture(params=IEEE_QUANTIZATION_MODES, ids=lambda m: str(m).split(".")[-1])
def ieee_quantization_mode(request):
    return request.param