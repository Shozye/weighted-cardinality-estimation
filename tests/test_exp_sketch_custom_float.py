"""Tests for FastExpSketchCustomFloat clone_with() feature.

Basic contract, merge, jaccard, pickle, and memory tests are covered by the
parametric suite via conftest SKETCH_SPECS.
"""

import pytest
from weighted_cardinality_estimation import FastExpSketchCustomFloat
from weighted_cardinality_estimation.stat import elements_stream

M = 64


class TestCloneWith:
    def test_clone_matches_fresh_sketch(self) -> None:
        m = 256
        s_high = FastExpSketchCustomFloat(m, seed=42, exp_bits=7, mant_bits=7)
        s_low = FastExpSketchCustomFloat(m, seed=42, exp_bits=6, mant_bits=5)
        elems = elements_stream(1000, seed=0)
        s_high.add_many(elems)
        s_low.add_many(elems)
        s_cloned = s_high.clone_with(exp_bits=6, mant_bits=5)
        assert s_cloned.estimate() == s_low.estimate()

    def test_clone_same_format_is_identity(self) -> None:
        s = FastExpSketchCustomFloat(M, seed=42, exp_bits=5, mant_bits=10)
        s.add_many(elements_stream(100, seed=0))
        clone = s.clone_with(exp_bits=5, mant_bits=10)
        assert clone.get_registers() == s.get_registers()

    def test_clone_raises_on_expansion(self) -> None:
        s = FastExpSketchCustomFloat(M, seed=42, exp_bits=5, mant_bits=10)
        with pytest.raises(ValueError):
            s.clone_with(exp_bits=8, mant_bits=10)
        with pytest.raises(ValueError):
            s.clone_with(exp_bits=5, mant_bits=11)

    def test_clone_preserves_properties(self) -> None:
        s = FastExpSketchCustomFloat(M, seed=42, exp_bits=7, mant_bits=7)
        s.add("x", weight=2.0)
        clone = s.clone_with(exp_bits=5, mant_bits=3)
        assert clone.exp_bits == 5
        assert clone.mant_bits == 3
        assert len(clone.get_registers()) == M

    def test_clone_jaccard_matches_fresh(self) -> None:
        m = 256
        s_a_high = FastExpSketchCustomFloat(m, seed=42, exp_bits=7, mant_bits=7)
        s_b_high = FastExpSketchCustomFloat(m, seed=42, exp_bits=7, mant_bits=7)
        s_a_low = FastExpSketchCustomFloat(m, seed=42, exp_bits=6, mant_bits=5)
        s_b_low = FastExpSketchCustomFloat(m, seed=42, exp_bits=6, mant_bits=5)

        shared = elements_stream(500, seed=0)
        only_a = elements_stream(300, seed=1)
        only_b = elements_stream(200, seed=2)

        s_a_high.add_many(shared)
        s_b_high.add_many(shared)
        s_a_low.add_many(shared)
        s_b_low.add_many(shared)
        s_a_high.add_many(only_a)
        s_a_low.add_many(only_a)
        s_b_high.add_many(only_b)
        s_b_low.add_many(only_b)

        s_a_cloned = s_a_high.clone_with(exp_bits=6, mant_bits=5)
        s_b_cloned = s_b_high.clone_with(exp_bits=6, mant_bits=5)
        assert s_a_cloned.jaccard_struct(s_b_cloned) == s_a_low.jaccard_struct(s_b_low)
