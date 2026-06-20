"""Tests for weighted_cardinality_estimation.stat module."""

import numpy as np
import pytest
from weighted_cardinality_estimation import stat


class TestMakeSeeds:
    def test_length(self) -> None:
        assert len(stat.make_seeds(10)) == 10

    def test_positive(self) -> None:
        seeds = stat.make_seeds(100)
        assert all(s > 0 for s in seeds)

    def test_deterministic(self) -> None:
        assert stat.make_seeds(10, seed=7) == stat.make_seeds(10, seed=7)

    def test_different_seeds_differ(self) -> None:
        assert stat.make_seeds(10, seed=1) != stat.make_seeds(10, seed=2)


class TestStream:
    def test_length(self) -> None:
        assert len(stat.elements_stream(50)) == 50

    def test_distinct(self) -> None:
        assert len(set(stat.elements_stream(100))) == 100


class TestWeightedStream:
    def test_length(self) -> None:
        elems, weights = stat.weighted_stream(50, 1000.0)
        assert len(elems) == 50
        assert len(weights) == 50

    def test_weights_sum(self) -> None:
        _, weights = stat.weighted_stream(100, 500.0, seed=0)
        assert abs(sum(weights) - 500.0) < 1e-8

    def test_all_positive(self) -> None:
        _, weights = stat.weighted_stream(100, 100.0, seed=0)
        assert all(w > 0 for w in weights)

    def test_deterministic(self) -> None:
        _, w1 = stat.weighted_stream(20, 10.0, seed=42)
        _, w2 = stat.weighted_stream(20, 10.0, seed=42)
        assert w1 == w2

    @pytest.mark.parametrize("dist_name", [
        "uniform", "exponential", "geometric", "constant", "pareto", "monotonic", "zipf",
    ])
    def test_string_alias(self, dist_name) -> None:
        _, weights = stat.weighted_stream(10, 100.0, dist=dist_name, seed=0)
        assert abs(sum(weights) - 100.0) < 1e-8

    def test_invalid_alias(self) -> None:
        with pytest.raises(ValueError, match="Unknown dist"):
            stat.weighted_stream(10, 100.0, dist="invalid")

    @pytest.mark.parametrize("dist", [
        stat.Uniform(), stat.Exponential(), stat.Geometric(),
        stat.Constant(), stat.Pareto(), stat.Monotonic(), stat.Zipf(),
    ])
    def test_all_distributions(self, dist) -> None:
        _, weights = stat.weighted_stream(20, 50.0, dist=dist, seed=0)
        assert abs(sum(weights) - 50.0) < 1e-8
        assert all(w > 0 for w in weights)


class TestJaccardStreams:
    @staticmethod
    def _true_weighted_jaccard(ea, wa, eb, wb):
        d_a = dict(zip(ea, wa, strict=False))
        d_b = dict(zip(eb, wb, strict=False))
        union = set(ea) | set(eb)
        s_min = sum(min(d_a.get(x, 0), d_b.get(x, 0)) for x in union)
        s_max = sum(max(d_a.get(x, 0), d_b.get(x, 0)) for x in union)
        return s_min / s_max if s_max > 0 else 0.0

    @pytest.mark.parametrize("target_j", [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0])
    @pytest.mark.parametrize("dist", ["uniform", "exponential"])
    @pytest.mark.parametrize("common", ["common_first", "common_last", "common_random"])
    def test_jaccard_accuracy(self, target_j, dist, common) -> None:
        (ea, wa), (eb, wb) = stat.jaccard_streams(
            500, 100.0, target_j, dist=dist, common=common, seed=0,
        )
        actual = self._true_weighted_jaccard(ea, wa, eb, wb)
        assert abs(actual - target_j) < 1e-6
        assert abs(sum(wa) - 100.0) < 1e-9
        assert abs(sum(wb) - 100.0) < 1e-9
        assert len(ea) == 500
        assert len(eb) == 500

    def test_deterministic_with_seed(self) -> None:
        r1 = stat.jaccard_streams(20, 100.0, 0.5, seed=7)
        r2 = stat.jaccard_streams(20, 100.0, 0.5, seed=7)
        assert r1 == r2

    def test_invalid_jaccard(self) -> None:
        with pytest.raises(ValueError):
            stat.jaccard_streams(10, 100.0, -0.1, seed=0)
        with pytest.raises(ValueError):
            stat.jaccard_streams(10, 100.0, 1.5, seed=0)


class TestRelativeError:
    def test_basic(self) -> None:
        assert stat.relative_error(1.1, 1.0) == pytest.approx(0.1)

    def test_symmetric(self) -> None:
        assert stat.relative_error(0.9, 1.0) == pytest.approx(0.1)

    def test_exact(self) -> None:
        assert stat.relative_error(5.0, 5.0) == 0.0


class TestSafeRelativeError:
    def test_finite(self) -> None:
        assert stat.safe_relative_error(1.1, 1.0) == pytest.approx(0.1)

    def test_nan(self) -> None:
        assert stat.safe_relative_error(float("nan"), 1.0) == 10.0

    def test_inf(self) -> None:
        assert stat.safe_relative_error(float("inf"), 1.0) == 10.0

    def test_non_positive(self) -> None:
        assert stat.safe_relative_error(-1.0, 1.0) == 10.0
        assert stat.safe_relative_error(0.0, 1.0) == 10.0

    def test_custom_cap(self) -> None:
        assert stat.safe_relative_error(float("nan"), 1.0, cap=5.0) == 5.0


class TestComputeRse:
    def test_zero_errors(self) -> None:
        assert stat.compute_rse([0.0, 0.0, 0.0]) == 0.0

    def test_known_value(self) -> None:
        # sqrt(mean([0.1^2, 0.2^2, 0.3^2])) = sqrt((0.01+0.04+0.09)/3)
        expected = float(np.sqrt(0.14 / 3))
        assert stat.compute_rse([0.1, 0.2, 0.3]) == pytest.approx(expected)

    def test_ignores_nan(self) -> None:
        result = stat.compute_rse([0.1, float("nan"), 0.1])
        assert np.isfinite(result)


class TestCompactVectorBytes:
    def test_exact(self) -> None:
        assert stat.compact_vector_bytes(8, 8) == 8

    def test_packing(self) -> None:
        assert stat.compact_vector_bytes(64, 5) == 40  # 64*5=320 bits = 40 bytes

    def test_rounding_up(self) -> None:
        assert stat.compact_vector_bytes(3, 5) == 2  # 15 bits -> 2 bytes


