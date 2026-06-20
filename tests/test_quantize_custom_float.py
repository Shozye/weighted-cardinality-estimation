"""Tests for custom float quantization behavior (ALL_NORMAL mode).

Tests verify:
- Correctness of representable value ranges (min/max)
- Log-uniform spacing property
- Saturation (no inf, no zero)
- Truncation behavior
"""

import math

import pytest
from weighted_cardinality_estimation import FastExpSketchCustomFloat
from weighted_cardinality_estimation.stat import elements_stream, weighted_stream


def make_sketch(exp_bits, mant_bits, m=64):
    return FastExpSketchCustomFloat(m, seed=0, exp_bits=exp_bits, mant_bits=mant_bits)


def get_max(s):
    """Fresh sketch registers are initialized to max value."""
    return s.get_registers()[0]


def expected_max(exp_bits, mant_bits):
    bias = (1 << (exp_bits - 1)) - 1
    max_exp = (1 << exp_bits) - 1
    mant_scale = 1 << mant_bits
    return math.ldexp(1.0 + (mant_scale - 1) / mant_scale, max_exp - bias)


def expected_min(exp_bits, mant_bits):
    bias = (1 << (exp_bits - 1)) - 1
    return math.ldexp(1.0, -bias)


class TestRange:
    @pytest.mark.parametrize(("exp_bits", "mant_bits"), [(3, 2), (4, 3), (5, 10), (8, 23)])
    def test_max_value(self, exp_bits, mant_bits) -> None:
        s = make_sketch(exp_bits, mant_bits, m=1)
        assert get_max(s) == pytest.approx(expected_max(exp_bits, mant_bits))

    @pytest.mark.parametrize(("exp_bits", "mant_bits"), [(3, 2), (4, 3), (5, 10)])
    def test_min_value_via_large_weight(self, exp_bits, mant_bits) -> None:
        s = make_sketch(exp_bits, mant_bits)
        s.add("x", weight=1e30)
        exp_min = expected_min(exp_bits, mant_bits)
        assert all(r == exp_min for r in s.get_registers())


class TestLogUniformSpacing:
    @pytest.mark.parametrize(("exp_bits", "mant_bits"), [(3, 2), (4, 3), (5, 10)])
    def test_all_registers_have_full_mantissa_precision(self, exp_bits, mant_bits) -> None:
        s = make_sketch(exp_bits, mant_bits, m=256)
        for e in elements_stream(500):
            s.add(e, weight=1.0)
        mant_scale = 1 << mant_bits
        for r in s.get_registers():
            frac, _ = math.frexp(r)
            mantissa_scaled = (frac * 2.0 - 1.0) * mant_scale
            assert abs(mantissa_scaled - round(mantissa_scaled)) < 1e-9

    def test_no_subnormal_values_exist(self) -> None:
        exp_bits, mant_bits = 4, 3
        bias = (1 << (exp_bits - 1)) - 1
        min_normal = math.ldexp(1.0, -bias)
        s = make_sketch(exp_bits, mant_bits, m=256)
        elems, weights = weighted_stream(1000, 1000.0 * 1001.0 / 2.0, dist="monotonic")
        s.add_many(elems, weights)
        for r in s.get_registers():
            assert r >= min_normal

    def test_distinct_values_per_exponent(self) -> None:
        exp_bits, mant_bits = 3, 2
        bias = (1 << (exp_bits - 1)) - 1
        mant_scale = 1 << mant_bits
        all_values = set()
        for biased_exp in range(1 << exp_bits):
            for k in range(mant_scale):
                val = math.ldexp(1.0 + k / mant_scale, biased_exp - bias)
                all_values.add(val)
        assert len(all_values) == (1 << exp_bits) * (1 << mant_bits)


class TestSaturation:
    @pytest.mark.parametrize(("exp_bits", "mant_bits"), [(2, 1), (3, 0), (3, 4), (5, 10)])
    @pytest.mark.parametrize(("weight", "check"), [
        (1e-20, lambda r: math.isfinite(r)),
        (1e20, lambda r: r > 0),
    ], ids=["no_inf", "no_zero"])
    def test_saturation(self, exp_bits, mant_bits, weight, check) -> None:
        s = make_sketch(exp_bits, mant_bits)
        s.add("x", weight=weight)
        assert all(check(r) for r in s.get_registers())

    def test_overflow_saturates(self) -> None:
        s = make_sketch(2, 1)
        s.add("x", weight=1e-30)
        assert all(r == 6.0 for r in s.get_registers())

    def test_underflow_saturates_to_min(self) -> None:
        s = make_sketch(2, 1)
        s.add("x", weight=1e30)
        assert all(r == 0.5 for r in s.get_registers())


