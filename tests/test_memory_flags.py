"""Tests for memory_usage(flags) API using MemoryFlag bitflags."""

import pytest
from weighted_cardinality_estimation import MemoryFlag


def test_nothing_returns_zero(sketch) -> None:
    assert sketch.memory_usage(MemoryFlag.NOTHING) == 0


@pytest.mark.parametrize("flag", [MemoryFlag.TOTAL, MemoryFlag.REGISTERS])
def test_flag_positive(sketch, flag) -> None:
    assert sketch.memory_usage(flag) > 0


@pytest.mark.parametrize(("subset", "superset"), [
    (MemoryFlag.REGISTERS, MemoryFlag.ALL_WRITE),
    (MemoryFlag.FISHER_YATES_PERM_INIT, MemoryFlag.FISHER_YATES),
])
def test_monotonicity(sketch, subset, superset) -> None:
    assert sketch.memory_usage(superset) >= sketch.memory_usage(subset)


def test_total_geq_individual(sketch, memory_flag) -> None:
    assert sketch.memory_usage(MemoryFlag.TOTAL) >= sketch.memory_usage(memory_flag)


def test_total_exclude_flag(sketch, memory_flag) -> None:
    total = sketch.memory_usage(MemoryFlag.TOTAL)
    component = sketch.memory_usage(memory_flag)
    assert sketch.memory_usage(MemoryFlag.TOTAL | memory_flag) == total - component


@pytest.mark.parametrize(("flag_a", "flag_b"), [
    (MemoryFlag.REGISTERS, MemoryFlag.SEEDS),
    (MemoryFlag.FISHER_YATES_PERM_INIT, MemoryFlag.FISHER_YATES_NON_PERM_INIT),
    (MemoryFlag.REGISTERS, MemoryFlag.FISHER_YATES_PERM_INIT),
    (MemoryFlag.SEEDS, MemoryFlag.FISHER_YATES_NON_PERM_INIT),
])
def test_additivity_disjoint_flags(sketch, flag_a, flag_b) -> None:
    combined = sketch.memory_usage(flag_a | flag_b)
    separate = sketch.memory_usage(flag_a) + sketch.memory_usage(flag_b)
    assert combined == separate
