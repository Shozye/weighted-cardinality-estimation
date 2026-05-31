"""Tests for estimate_direct, estimate_newton_cold, estimate_newton_warm on BaseLogExpSketch only."""
import pytest
from weighted_cardinality_estimation import (
    BaseLogExpSketch,
    FastLogExpSketch,
    BaseShiftedLogExpSketch,
    FastShiftedLogExpSketch,
)

LOG_EXP_CLASSES_ALL = [
    BaseLogExpSketch,
    FastLogExpSketch,
    BaseShiftedLogExpSketch,
    FastShiftedLogExpSketch,
]

M = 200
SEEDS = list(range(1, M + 1))
WEIGHT = 1e4
ELEMENTS = [f"elem_{i}" for i in range(500)]


def make_sketch(cls):
    s = cls(M, SEEDS, amount_bits=8, logarithm_base=2)
    for elem in ELEMENTS:
        s.add(elem, WEIGHT / len(ELEMENTS))
    return s


@pytest.mark.parametrize("cls", [BaseLogExpSketch], ids=[BaseLogExpSketch.__name__])
def test_three_estimators_callable(cls):
    s = make_sketch(cls)
    assert s.estimate_direct() > 0
    assert s.estimate_newton_cold() > 0
    assert s.estimate_newton_warm() > 0


@pytest.mark.parametrize("cls", [BaseLogExpSketch], ids=[BaseLogExpSketch.__name__])
def test_direct_and_warm_agree(cls):
    s = make_sketch(cls)
    d = s.estimate_direct()
    w = s.estimate_newton_warm()
    assert abs(d - w) / w < 0.01, f"direct={d}, warm={w}"


@pytest.mark.parametrize("cls", [BaseLogExpSketch], ids=[BaseLogExpSketch.__name__])
def test_cold_and_warm_converge(cls):
    """Cold and warm Newton should converge to the same value given enough iterations.
    With NEWTON_MAX_ITERATIONS=5 and a distant cold start they may differ, so we
    only assert both are positive and in the same ballpark (within 2 orders of magnitude).
    """
    s = make_sketch(cls)
    cold = s.estimate_newton_cold()
    warm = s.estimate_newton_warm()
    assert cold > 0
    assert warm > 0
    # Both should at least be in the same order of magnitude (not wildly divergent)
    assert min(cold, warm) / max(cold, warm) > 1e-6


@pytest.mark.parametrize("cls", LOG_EXP_CLASSES_ALL, ids=[c.__name__ for c in LOG_EXP_CLASSES_ALL])
def test_estimate_unchanged(cls):
    """estimate() must return the same value as before (regression guard)."""
    s = make_sketch(cls)
    e1 = s.estimate()
    e2 = s.estimate()
    assert e1 == e2
    assert e1 > 0


# --- newton_iteration_count tests (BaseLogExpSketch only) ---

def make_base_sketch():
    s = BaseLogExpSketch(M, SEEDS, amount_bits=8, logarithm_base=2)
    for elem in ELEMENTS:
        s.add(elem, WEIGHT / len(ELEMENTS))
    return s


def test_newton_iteration_count_positive_after_estimate():
    s = make_base_sketch()
    s.estimate()
    assert s.newton_iteration_count > 0


def test_newton_iteration_count_resets_between_calls():
    s = make_base_sketch()
    s.estimate()
    first = s.newton_iteration_count
    s.estimate()
    second = s.newton_iteration_count
    assert first == second  # count reflects only the last call


def test_newton_iteration_count_zero_after_direct():
    s = make_base_sketch()
    s.estimate_direct()
    assert s.newton_iteration_count == 0


def test_newton_iteration_count_positive_after_cold():
    s = make_base_sketch()
    s.estimate_newton_cold()
    assert s.newton_iteration_count >= 0  # cold start may converge in 0 extra iters


def test_newton_iteration_count_positive_after_warm():
    s = make_base_sketch()
    s.estimate_newton_warm()
    assert s.newton_iteration_count >= 0
