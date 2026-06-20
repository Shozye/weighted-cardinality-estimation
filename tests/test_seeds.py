"""Tests for Seeds(master_seed) via lazy splitmix64 derivation."""

import numpy as np
from weighted_cardinality_estimation import ExpSketch


def splitmix64(x: int) -> int:
    x = (x + 0x9E3779B97F4A7C15) & 0xFFFFFFFFFFFFFFFF
    x = ((x ^ (x >> 30)) * 0xBF58476D1CE4E5B9) & 0xFFFFFFFFFFFFFFFF
    x = ((x ^ (x >> 27)) * 0x94D049BB133111EB) & 0xFFFFFFFFFFFFFFFF
    return (x ^ (x >> 31)) & 0xFFFFFFFFFFFFFFFF


def derive_seeds(master_seed: int, count: int) -> list[int]:
    return [splitmix64(master_seed + i) & 0xFFFFFFFF for i in range(count)]


class TestSeedsDerived:
    def test_deterministic(self) -> None:
        s1 = ExpSketch(100, seed=42)
        s2 = ExpSketch(100, seed=42)
        s1.add("hello", 5.0)
        s2.add("hello", 5.0)
        assert s1.estimate() == s2.estimate()

    def test_all_distinct(self) -> None:
        seeds = derive_seeds(123, 400)
        assert len(set(seeds)) == 400

    def test_different_master_seeds_differ(self) -> None:
        s1 = derive_seeds(0, 100)
        s2 = derive_seeds(1, 100)
        assert s1 != s2

    def test_reasonable_estimates(self) -> None:
        sketch = ExpSketch(400, seed=7)
        for i in range(1000):
            sketch.add(f"elem_{i}", 1.0)
        est = sketch.estimate()
        assert 800 < est < 1200
