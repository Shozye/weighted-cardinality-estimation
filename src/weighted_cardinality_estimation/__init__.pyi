from . import MemoryFlag as MemoryFlag
from . import stat as stat

# ─── MemoryFlag constants ─────────────────────────────────────────────────────

class MemoryFlag:
    NOTHING: int
    REGISTERS: int
    ALL_WRITE_NO_REGISTERS: int
    FISHER_YATES_PERM_INIT: int
    FISHER_YATES_NON_PERM_INIT: int
    SEEDS: int
    TOTAL: int
    ALL_WRITE: int
    FISHER_YATES: int


# ─── RngEngine enum ──────────────────────────────────────────────────────────

class RngEngine:
    PCG64: RngEngine
    MT19937: RngEngine
    XOSHIRO128PP: RngEngine
    XOSHIRO256PP: RngEngine


# ─── Base ─────────────────────────────────────────────────────────────────────

class CardinalitySketch:


    def add(self, x: str) -> None:
        ...
    def add_many(self, elems: list[str]) -> None:
        ...
    def estimate(self) -> float:
        ...
    def memory_usage(self, flags: int) -> int: ...


# ─── Mixins ───────────────────────────────────────────────────────────────────

class WeightedMixin:


    def add(self, x: str, weight: float = ...) -> None:
        ...
    def add_many(self, elems: list[str], weights: list[float] = ...) -> None:
        ...

class MergeableMixin:


    def merge(self, other: MergeableMixin) -> None:
        ...

class JaccardMixin:


    def jaccard_struct(self, other: JaccardMixin) -> float:
        ...

class NewtonMixin:


    def estimate_direct(self) -> float: ...
    def estimate_newton_cold(self) -> float: ...
    def estimate_newton_warm(self) -> float: ...
    def estimate_newton_cold_iterations(self) -> int: ...
    def estimate_newton_warm_iterations(self) -> int: ...


# ─── Unweighted sketches ──────────────────────────────────────────────────────

class MartingaleMinHash(CardinalitySketch):


    def __init__(self, m: int, seed: int) -> None: ...

class HyperLogLog(MergeableMixin, CardinalitySketch):


    def __init__(self, m: int, seed: int) -> None: ...

class MinHash(JaccardMixin, MergeableMixin, CardinalitySketch):


    def __init__(self, m: int, seed: int) -> None: ...


# ─── Weighted sketches ────────────────────────────────────────────────────────

class MartingaleWeightedMinHash(WeightedMixin, CardinalitySketch):


    def __init__(self, m: int, seed: int) -> None: ...

class WeightedMinHash(MergeableMixin, WeightedMixin, CardinalitySketch):


    def __init__(self, m: int, seed: int) -> None: ...

class WeightedHyperLogLog(MergeableMixin, WeightedMixin, CardinalitySketch):


    def __init__(self, m: int, seed: int) -> None: ...

class WeightedHyperLogLogFloat32(MergeableMixin, WeightedMixin, CardinalitySketch):


    def __init__(self, m: int, seed: int) -> None: ...

class QSketch(MergeableMixin, WeightedMixin, CardinalitySketch):


    def __init__(self, m: int, seed: int, amount_bits: int, rng_engine: RngEngine = ...) -> None: ...

class QSketchDyn(MergeableMixin, WeightedMixin, CardinalitySketch):


    def __init__(self, m: int, seed: int, amount_bits: int, g_seed: int = ...) -> None: ...

class kQSketch(NewtonMixin, MergeableMixin, WeightedMixin, CardinalitySketch):


    def __init__(self, m: int, seed: int, amount_bits: int, logarithm_base: float, rng_engine: RngEngine = ...) -> None: ...

class kQSketchRounding(NewtonMixin, MergeableMixin, WeightedMixin, CardinalitySketch):


    def __init__(self, m: int, seed: int, amount_bits: int, logarithm_base: float, rng_engine: RngEngine = ...) -> None: ...
    def estimate_corrected(self) -> float: ...

class kQSketchRoundedDyn(NewtonMixin, WeightedMixin, CardinalitySketch):


    def __init__(self, m: int, seed: int, amount_bits: int, logarithm_base: float, g_seed: int = ...) -> None: ...

class kQSketchShifted(NewtonMixin, MergeableMixin, WeightedMixin, CardinalitySketch):


    def __init__(self, m: int, seed: int, amount_bits: int, logarithm_base: float, rng_engine: RngEngine = ...) -> None: ...
    def get_offset(self) -> int: ...

class ExpSketch(JaccardMixin, MergeableMixin, WeightedMixin, CardinalitySketch):


    def __init__(self, m: int, seed: int) -> None: ...

class ExpSketchFloat32(JaccardMixin, MergeableMixin, WeightedMixin, CardinalitySketch):


    def __init__(self, m: int, seed: int) -> None: ...

class FastExpSketchFloat32(JaccardMixin, MergeableMixin, WeightedMixin, CardinalitySketch):


    def __init__(self, m: int, seed: int, rng_engine: RngEngine = ...) -> None: ...

class FastExpSketch(JaccardMixin, MergeableMixin, WeightedMixin, CardinalitySketch):


    def __init__(self, m: int, seed: int, rng_engine: RngEngine = ...) -> None: ...

class FastGMExpSketch(JaccardMixin, MergeableMixin, WeightedMixin, CardinalitySketch):


    def __init__(self, m: int, seed: int, rng_engine: RngEngine = ...) -> None: ...


class LogExpSketchSlowNoShifted(JaccardMixin, MergeableMixin, WeightedMixin, CardinalitySketch):


    def __init__(self, m: int, seed: int, amount_bits: int, v_max: float) -> None: ...


class LogExpSketchSlowShifted(JaccardMixin, MergeableMixin, WeightedMixin, CardinalitySketch):


    def __init__(self, m: int, seed: int, amount_bits: int, v_max: float) -> None: ...
    def get_offset(self) -> int: ...


class LogExpSketchFastNoShifted(JaccardMixin, MergeableMixin, WeightedMixin, CardinalitySketch):


    def __init__(self, m: int, seed: int, amount_bits: int, v_max: float, rng_engine: RngEngine = ...) -> None: ...


class LogExpSketchFastShifted(JaccardMixin, MergeableMixin, WeightedMixin, CardinalitySketch):


    def __init__(self, m: int, seed: int, amount_bits: int, v_max: float, rng_engine: RngEngine = ...) -> None: ...
    def get_offset(self) -> int: ...


# ─── Custom float ExpSketch ───────────────────────────────────────────────────

class QuantizationMode:


    ALL_NORMAL: QuantizationMode
    WITH_SUBNORMALS: QuantizationMode
    LINEAR: QuantizationMode
    LOGARITHMIC: QuantizationMode

class FastExpSketchCustomFloat(JaccardMixin, MergeableMixin, WeightedMixin, CardinalitySketch):


    exp_bits: int
    mant_bits: int

    def __init__(self, m: int, seed: int, exp_bits: int, mant_bits: int, rng_engine: RngEngine = ...) -> None: ...
    def clone_with(self, exp_bits: int, mant_bits: int) -> FastExpSketchCustomFloat: ...


class WeightedHyperLogLogCustomFloat(MergeableMixin, WeightedMixin, CardinalitySketch):


    exp_bits: int
    mant_bits: int

    def __init__(self, m: int, seed: int, exp_bits: int, mant_bits: int) -> None: ...
