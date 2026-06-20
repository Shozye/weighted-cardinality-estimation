"""Basic contract: universal properties every sketch must satisfy."""

import pytest
import weighted_cardinality_estimation as wce
from conftest import M
from weighted_cardinality_estimation import (
    CardinalitySketch,
    JaccardMixin,
    MergeableMixin,
    WeightedMixin,
)
from weighted_cardinality_estimation.stat import relative_error, elements_stream, weighted_stream

class TestBasic:
    def test_add_and_estimate_positive(self, sketch) -> None:
        sketch.add("element")
        assert sketch.estimate() > 0

    def test_duplicate_add_does_not_change_estimate(self, sketch) -> None:
        """Idempotency: adding the same element twice doesn't change the estimate."""
        sketch.add("element")
        estimate_before = sketch.estimate()
        sketch.add("element")
        assert sketch.estimate() == estimate_before

    def test_add_many(self, sketch) -> None:
        sketch.add_many(elements_stream(10))
        assert sketch.estimate() > 0

    def test_estimate_empty(self, sketch) -> None:
        est = sketch.estimate()
        assert isinstance(est, float)

    def test_add_many_empty_list(self, sketch) -> None:
        sketch.add_many([])

    @pytest.mark.parametrize("m", [1, 4], ids=["m=1", "m=4"])
    def test_small_m_does_not_crash(self, spec, m) -> None:
        """Small m values (including non-power-of-two) must not crash."""
        sketch = spec.factory(m)
        sketch.add("x")
        sketch.estimate()

    def test_empty_string_element(self, spec) -> None:
        """Empty string as element should work like any other element."""
        sketch = spec.factory(16)
        sketch.add("")
        assert sketch.estimate() > 0

    def test_monotonicity(self, spec) -> None:
        """Estimate must not decrease as more distinct elements are added."""
        if "WeightedHyperLogLog" in spec.name:
            pytest.skip(
                "WeightedHyperLogLog estimator is non-monotone: new buckets add unbounded "
                "-log(u) terms to the denominator, which can decrease the estimate",
            )
        sketch = spec.factory(M, 7)
        prev = 0.0
        for step, elem in enumerate(elements_stream(100, seed=0)):
            sketch.add(elem)
            est = sketch.estimate()
            assert est >= prev, (
                f"{spec.name}: estimate dropped from {prev} to {est} at step {step}"
            )
            prev = est

    def test_cardinality_accuracy(self, spec) -> None:
        """Median relative error over 3 independent runs must be within spec's tolerance."""
        errors = []
        for seed in range(3):
            sketch = spec.factory(400)
            sketch.add_many(elements_stream(500, seed=seed))
            errors.append(relative_error(sketch.estimate(), 500))
        errors.sort()
        med_err = errors[1]
        assert med_err <= spec.estimate_rel_error, (
            f"{spec.name}: median_rel_err={med_err:.2%} > {spec.estimate_rel_error:.0%}"
            f" (errors={[f'{e:.2%}' for e in errors]})"
        )


class TestAddMany:
    def test_equals_add_one_by_one(self, spec) -> None:
        """add_many must give identical result to sequential add calls (unweighted)."""
        elems = list(elements_stream(50, seed=0))
        s1 = spec.factory(M, 99)
        for e in elems:
            s1.add(e)
        s2 = spec.factory(M, 99)
        s2.add_many(elems)
        assert s1.estimate() == s2.estimate()

    def test_equals_add_one_by_one_weighted(self, weighted_spec) -> None:
        """add_many must give identical result to sequential add calls (weighted)."""
        elems, weights = weighted_stream(50, 1275.0, dist="monotonic", seed=0)
        s1 = weighted_spec.factory(M, 99)
        for e, w in zip(elems, weights, strict=False):
            s1.add(e, weight=w)
        s2 = weighted_spec.factory(M, 99)
        s2.add_many(elems, weights)
        assert s1.estimate() == s2.estimate()

    def test_add_many_weighted(self, weighted_spec) -> None:
        sketch = weighted_spec.factory(M)
        elems, weights = weighted_stream(10, 10.0, dist="constant", seed=0)
        sketch.add_many(elems, weights)
        assert sketch.estimate() > 0


class TestHierarchy:
    def test_is_cardinality_sketch(self, sketch) -> None:
        assert isinstance(sketch, CardinalitySketch)

    def test_weighted_mixin_iff_is_weighted(self, sketch, spec) -> None:
        assert isinstance(sketch, WeightedMixin) == spec.is_weighted

    def test_mergeable_mixin_iff_has_merge(self, sketch) -> None:
        assert isinstance(sketch, MergeableMixin) == hasattr(sketch, "merge")

    def test_jaccard_mixin_iff_has_jaccard_struct(self, sketch) -> None:
        assert isinstance(sketch, JaccardMixin) == hasattr(sketch, "jaccard_struct")


class TestRegisters:
    def test_get_registers_length(self, sketch) -> None:
        assert len(sketch.get_registers()) == M

    def test_get_registers_changes_after_add(self, sketch) -> None:
        regs_before = sketch.get_registers()
        sketch.add_many(elements_stream(10))
        assert sketch.get_registers() != regs_before


class TestWeightedEdgeCases:
    @pytest.mark.parametrize("weight", [0.0, -1.0, float("nan"), float("inf")])
    def test_invalid_weight_raises(self, weighted_spec, weight) -> None:
        sketch = weighted_spec.factory(16)
        with pytest.raises((ValueError, RuntimeError)):
            sketch.add("x", weight=weight)

    def test_add_many_mismatched_lengths(self, weighted_spec) -> None:
        sketch = weighted_spec.factory(16)
        with pytest.raises((ValueError, RuntimeError)):
            sketch.add_many(["a", "b", "c"], [1.0, 2.0])


_AMOUNT_BITS_ONE_SEEDS = list(range(1, 17))


@pytest.mark.parametrize("factory", [
    lambda: wce.QSketch(16, 42, 1),
    lambda: wce.QSketchDyn(16, 42, 1),
    lambda: wce.kQSketch(16, 42, 1, 2.0),
    lambda: wce.kQSketchRounding(16, 42, 1, 2.0),
    lambda: wce.kQSketchShifted(16, 42, 1, 2.0),
], ids=["QSketch", "QSketchDyn", "kQSketch",
        "kQSketchRounding", "kQSketchShifted"])
def test_amount_bits_one(factory) -> None:
    """amount_bits=1 should raise ValueError."""
    with pytest.raises((ValueError, RuntimeError)):
        factory()
