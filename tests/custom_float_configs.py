"""Standard custom float configurations for testing.

8 configs organized as 2 range levels × 4 accuracy levels.
No sign bit (sketch registers are always positive), so total_bits = exp_bits + mant_bits.

Range levels:
    f32Range: exp_bits=8  (bias=127,  range ≈ 1e-38 to 1e38)
    f64Range: exp_bits=11 (bias=1023, range ≈ 1e-308 to 1e308)

Accuracy levels (mantissa bits):
    XS: mant_bits=7   (bfloat16-level precision)
    S:  mant_bits=10  (float16-level precision)
    M:  mant_bits=15  (medium precision)
    L:  mant_bits=23  (float32-level precision)
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class CustomFloatConfig:
    exp_bits: int
    mant_bits: int
    min_weight: float
    max_weight: float

    @property
    def total_bits(self) -> int:
        return self.exp_bits + self.mant_bits

    @property
    def label(self) -> str:
        return f"q{self.exp_bits}_p{self.mant_bits}"


# ─── f32 Range (exp_bits=8, range ≈ 1e-38 to 1e38) ──────────────────────────

F32_XS = CustomFloatConfig(exp_bits=8, mant_bits=7, min_weight=1e-37, max_weight=1e38)
F32_S = CustomFloatConfig(exp_bits=8, mant_bits=10, min_weight=1e-37, max_weight=1e38)
F32_M = CustomFloatConfig(exp_bits=8, mant_bits=15, min_weight=1e-37, max_weight=1e38)
F32_L = CustomFloatConfig(exp_bits=8, mant_bits=23, min_weight=1e-37, max_weight=1e38)

# ─── f64 Range (exp_bits=11, range ≈ 1e-308 to 1e308) ────────────────────────

F64_XS = CustomFloatConfig(exp_bits=11, mant_bits=7, min_weight=1e-305, max_weight=1e307)
F64_S = CustomFloatConfig(exp_bits=11, mant_bits=10, min_weight=1e-305, max_weight=1e307)
F64_M = CustomFloatConfig(exp_bits=11, mant_bits=15, min_weight=1e-305, max_weight=1e307)
F64_L = CustomFloatConfig(exp_bits=11, mant_bits=23, min_weight=1e-305, max_weight=1e307)

# ─── All configs for iteration ────────────────────────────────────────────────

ALL_CONFIGS = [F32_XS, F32_S, F32_M, F32_L, F64_XS, F64_S, F64_M, F64_L]

# Subset for SKETCH_SPECS (one from each range level, moderate precision)
SKETCH_SPEC_CONFIGS = [F32_XS, F32_M, F64_XS, F64_M]
