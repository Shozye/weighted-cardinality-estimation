import copy
import random
import pytest
from weighted_cardinality_estimation import BaseQSketch, FastExpSketch, ExpSketch, FastQSketch, QSketchDyn


SKETCH_CONSTRUCTORS = [
    pytest.param(ExpSketch, id="ExpSketch"),
    pytest.param(FastExpSketch, id="FastExpSketch"),
    pytest.param(lambda m, seeds: BaseQSketch(m, seeds, 8), id="BaseQSketch"),
    pytest.param(lambda m, seeds: FastQSketch(m, seeds, amount_bits=8), id="FastQSketch"),
    pytest.param(lambda m, seeds: QSketchDyn(m, seeds, amount_bits=8, g_seed=42), id="QSketchDyn"),
]

@pytest.mark.parametrize("sketch_cls", SKETCH_CONSTRUCTORS)
def test_unitary(sketch_cls):
    M=5
    seeds = [random.randint(1,10000000) for _ in range(M)]
    sketch = sketch_cls(M, seeds)
    sketch.add("I am just a simple element.", weight=1)

    estimate = sketch.estimate()
    assert estimate > 0

@pytest.mark.parametrize("sketch_cls", SKETCH_CONSTRUCTORS)
def test_estimate_adding_duplicate_does_not_change_estimation(sketch_cls):
    M=5
    seeds = [random.randint(1,10000000) for _ in range(M)]
    sketch = sketch_cls(M, seeds)
    sketch.add("I am just a simple element.", weight=1)

    estimate = sketch.estimate()
    sketch.add("I am just a simple element.", weight=1)

    assert estimate == sketch.estimate()

@pytest.mark.parametrize("sketch_cls", SKETCH_CONSTRUCTORS)
def test_copy_produces_same_estimate(sketch_cls):
    # here i just want to make some basic contract that it holds to ANY standard lol
    m = 5
    seeds = [1, 2, 3, 4, 5]
    original_sketch = sketch_cls(m, seeds)
    original_sketch.add("A single test element", weight=1.0)
    original_estimate = original_sketch.estimate()
    
    copied_sketch = copy.copy(original_sketch)
    copied_estimate = copied_sketch.estimate()

    assert original_estimate == copied_estimate
    assert original_estimate > 0
    assert original_sketch is not copied_sketch

@pytest.mark.parametrize("sketch_cls", SKETCH_CONSTRUCTORS)
def test_copy_has_identical_internal_state(sketch_cls):
    m = 5
    seeds = [1, 2, 3, 4, 5]
    original_sketch = sketch_cls(m, seeds)
    original_sketch.add("some element", weight=2.5)

    copied_sketch = copy.copy(original_sketch)

    # this test is essentialy, from testing point of view, pointless, but is good for dev.
    original_state = original_sketch.__getstate__()
    copied_state = copied_sketch.__getstate__()

    assert original_state == copied_state

@pytest.mark.parametrize("sketch_cls", SKETCH_CONSTRUCTORS)
def test_copy_different_memory_objects(sketch_cls):
    m = 5
    seeds = [1, 2, 3, 4, 5]
    original_sketch = sketch_cls(m, seeds)
    original_sketch.add("some element", weight=2.5)

    copied_sketch = copy.copy(original_sketch)
    # here i want to check if internal stuff is different so it is important to save it
    original_estimate = original_sketch.estimate() 
    copied_sketch.add("new element", weight=1)
    assert original_estimate == original_sketch.estimate()

@pytest.mark.parametrize("sketch_cls", SKETCH_CONSTRUCTORS)
def test_copy_independently_the_same_structures(sketch_cls):
    m = 5
    seeds = [1, 2, 3, 4, 5]
    original_sketch = sketch_cls(m, seeds)
    original_sketch.add("some element", weight=2.5)

    copied_sketch = copy.copy(original_sketch)
    # here I want to see if all internal stuff is the same so adding new element will cause same effect.
    copied_sketch.add("new element", weight=1)
    original_sketch.add("new element", weight=1)

    assert original_sketch.estimate() == copied_sketch.estimate()

@pytest.mark.parametrize("sketch_cls", SKETCH_CONSTRUCTORS)
def test_memory_usage_sanity_check(sketch_cls):
    m = 5
    seeds = [1, 2, 3, 4, 5]
    original_sketch = sketch_cls(m, seeds)
    original_sketch.add("some element", weight=2.5)

    total_memory = original_sketch.memory_usage_total()
    write_memory = original_sketch.memory_usage_write()
    estimate_memory = original_sketch.memory_usage_estimate()
    assert total_memory > write_memory
    assert write_memory >= estimate_memory
    assert estimate_memory > 0
