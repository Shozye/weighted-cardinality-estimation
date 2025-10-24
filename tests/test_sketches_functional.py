import random
import pytest
from tests.utils import _SKETCH_CLS_ALL, _SKETCH_CLS_JACC, _SKETCH_CLS_NO_SEEDS, _SKETCH_CLS_SEEDS, assert_error

M_SIZE = 400
AMOUNT_ELEMENTS = 1000
ELEMENTS_WEIGHT = 10.0
AMOUNT_TEST_RUNS = 50

def _make_params(sketch_clss: dict, allowed_error: float) -> list:
    return [
        pytest.param(sketch_cls, allowed_error, id=sketch_name) for sketch_name, sketch_cls in sketch_clss.items()
    ]

SKETCH_PARAMS = _make_params(_SKETCH_CLS_SEEDS, 0.1)
SKETCH_PARAMS_NO_SEEDS = _make_params(_SKETCH_CLS_NO_SEEDS, 0.2)
SKETCH_PARAMS_JACCARD = _make_params(_SKETCH_CLS_JACC, 0.7)

@pytest.mark.parametrize("sketch_cls, allowed_error", SKETCH_PARAMS)
def test_sketch_functional_accuracy(
        sketch_cls,
        allowed_error: float,
    ):
    # here I want to assert some level of error on the structures to make sure they are any good

    total_weight = AMOUNT_ELEMENTS * ELEMENTS_WEIGHT
    estimates = []
    for _ in range(AMOUNT_TEST_RUNS):
        seeds = [random.randint(1,10000000) for _ in range(M_SIZE)]
        s = sketch_cls(M_SIZE, seeds)
        elements = [f"e{i}" for i in range(AMOUNT_ELEMENTS)]
        weights = [ELEMENTS_WEIGHT] * AMOUNT_ELEMENTS
        s.add_many(elements, weights)
        estimates.append(s.estimate())
    
    average_estimate = sum(estimates)/len(estimates)
    assert_error(total_weight, average_estimate, allowed_error)

@pytest.mark.parametrize("sketch_cls, allowed_error", SKETCH_PARAMS_NO_SEEDS)
def test_sketch_no_seeds_functional_accuracy(
        sketch_cls,
        allowed_error: float,
    ):
    # here I want to assert some level of error on the structures to make sure they are any good

    total_weight = AMOUNT_ELEMENTS * ELEMENTS_WEIGHT

    s = sketch_cls(M_SIZE, [])
    elements = [f"e{i}" for i in range(AMOUNT_ELEMENTS)]
    weights = [ELEMENTS_WEIGHT] * AMOUNT_ELEMENTS
    s.add_many(elements, weights)
    
    estimate = s.estimate()
    assert_error(total_weight, estimate, allowed_error)


@pytest.mark.parametrize("sketch_cls, allowed_error", SKETCH_PARAMS_JACCARD)
def test_sketch_functional_jaccard(
        sketch_cls,
        allowed_error: float,
    ):
    POSSIBLE_JACCARDS = [x/20 for x in range(0, 21)] # 0, 0.05, 0.1, 0.15...

    TOTAL_ELEMS = 1000
    TOTAL_SET_A = [f"a{i}a{i}a{i}" for i in range(100_000)]
    TOTAL_SET_B = [f"b{i}b{i}b{i}" for i in range(100_000)]
    TOTAL_SET_COMMON = [f"c{i}c{i}c{i}" for i in range(100_000)]

    for possible_jaccard in POSSIBLE_JACCARDS:
        AMOUNT_COMMON =  int(TOTAL_ELEMS * possible_jaccard)
        AMOUNT_NOT_COMMON = TOTAL_ELEMS - AMOUNT_COMMON

        set_common = random.sample(TOTAL_SET_COMMON, AMOUNT_COMMON)
        set_a = random.sample(TOTAL_SET_A, AMOUNT_NOT_COMMON) + set_common
        set_b = random.sample(TOTAL_SET_B, AMOUNT_NOT_COMMON) + set_common

        seeds = [random.randint(1,10_000_000) for _ in range(M_SIZE)]
        sketch_a = sketch_cls(M_SIZE, seeds)
        sketch_b = sketch_cls(M_SIZE, seeds)

        weights = [ELEMENTS_WEIGHT] * AMOUNT_ELEMENTS
        sketch_a.add_many(set_a, weights)
        sketch_b.add_many(set_b, weights)
        estimated_jaccard = sketch_a.jaccard_struct(sketch_b)
        assert_error(possible_jaccard, estimated_jaccard, allowed_error, delta=0.05)
