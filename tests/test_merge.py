"""Merge: merged sketch estimates the union cardinality."""

import pytest
import weighted_cardinality_estimation as wce
from conftest import M, make_sketches
from weighted_cardinality_estimation.stat import elements_stream, weighted_stream


def test_merge_equivalence(merge_spec) -> None:
    """merge(A,B) should give same estimate as a single sketch that saw both streams."""
    combined, sketch_a, sketch_b = make_sketches(merge_spec, M, seed=99, n=3)
    elems_a = list(elements_stream(50, seed=0))
    elems_b = list(elements_stream(50, seed=1))
    combined.add_many(elems_a)
    combined.add_many(elems_b)
    sketch_a.add_many(elems_a)
    sketch_b.add_many(elems_b)
    sketch_a.merge(sketch_b)
    assert sketch_a.estimate() == combined.estimate()


def test_merge_with_empty(merge_spec) -> None:
    """Merging with an empty sketch should not change estimate."""
    sketch, empty = make_sketches(merge_spec, M, seed=0)
    sketch.add_many(elements_stream(20))
    est_before = sketch.estimate()
    sketch.merge(empty)
    assert sketch.estimate() == pytest.approx(est_before, rel=1e-9)


def test_merge_size_mismatch(merge_spec) -> None:
    """Merging sketches of different sizes should raise."""
    small = merge_spec.factory(32)
    big = merge_spec.factory(M)
    with pytest.raises((ValueError, RuntimeError)):
        big.merge(small)


def test_merge_idempotent(merge_spec) -> None:
    """Merging a sketch with itself should not change estimate."""
    sketch = merge_spec.factory(M)
    sketch.add_many(elements_stream(30))
    est_before = sketch.estimate()
    sketch.merge(sketch)
    assert sketch.estimate() == pytest.approx(est_before, rel=1e-9)


def test_merge_weighted(weighted_merge_spec) -> None:
    """merge(A,B) gives same estimate as combined sketch that saw both weighted streams."""
    combined, sketch_a, sketch_b = make_sketches(weighted_merge_spec, M, seed=77, n=3)
    elems_a, weights_a = weighted_stream(30, total_weight=100.0, seed=10)
    elems_b, weights_b = weighted_stream(30, total_weight=200.0, seed=11)
    combined.add_many(elems_a, weights_a)
    combined.add_many(elems_b, weights_b)
    sketch_a.add_many(elems_a, weights_a)
    sketch_b.add_many(elems_b, weights_b)
    sketch_a.merge(sketch_b)
    assert sketch_a.estimate() == combined.estimate()


def test_merge_associative(merge_spec) -> None:
    """(A merged B) merged C == A merged (B merged C)."""
    elems_a = list(elements_stream(30, seed=20))
    elems_b = list(elements_stream(30, seed=21))
    elems_c = list(elements_stream(30, seed=22))

    a1, b1, c1 = make_sketches(merge_spec, M, seed=55, n=3)
    a1.add_many(elems_a); b1.add_many(elems_b); c1.add_many(elems_c)
    a1.merge(b1); a1.merge(c1)

    a2, b2, c2 = make_sketches(merge_spec, M, seed=55, n=3)
    a2.add_many(elems_a); b2.add_many(elems_b); c2.add_many(elems_c)
    b2.merge(c2); a2.merge(b2)

    assert a1.estimate() == a2.estimate()


@pytest.mark.parametrize("factory", [
    lambda q, p: wce.FastExpSketchCustomFloat(64, seed=42, exp_bits=q, mant_bits=p),
    lambda q, p: wce.WeightedHyperLogLogCustomFloat(64, seed=42, exp_bits=q, mant_bits=p),
], ids=["FastExpSketchCustomFloat", "WeightedHyperLogLogCustomFloat"])
def test_merge_format_mismatch(factory) -> None:
    """Merging custom float sketches with different formats should raise."""
    s1 = factory(5, 10)
    s2 = factory(8, 23)
    with pytest.raises(ValueError):
        s1.merge(s2)
