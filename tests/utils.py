from typing import Callable

from weighted_cardinality_estimation import BaseLogExpSketch, BaseLogExpSketchJacc, BaseQSketch, ExpSketchFloat, FastExpSketch, ExpSketch, FastGMExpSketch, BaseLogExpSketch, FastLogExpSketch, FastQSketch, QSketchDyn, BaseShiftedLogExpSketch, FastShiftedLogExpSketch, QSketch

def assert_error(expected: float, actual: float, error: float, delta: float=0):
    min_range = (expected-delta) * (1-error)
    max_range = (expected+delta) * (1+error)
    assert min_range <= actual <= max_range, f"{min_range=} <= {actual=} <= {max_range=}"

_SKETCH_CLS_JACC_SEEDS = {
    "ExpSketch": lambda m, seeds: ExpSketch(m, seeds),
    "ExpSketchFloat": lambda m, seeds: ExpSketchFloat(m, seeds),
    "FastExpSketch": lambda m, seeds: FastExpSketch(m, seeds),
    "FastGMExpSketch": lambda m, seeds: FastGMExpSketch(m, seeds),
    "BaseLogExpSketchJacc": lambda m, seeds: BaseLogExpSketchJacc(m, seeds, amount_bits=8, logarithm_base=2, amount_bits_jaccard=8),
}

_SKETCH_CLS_NO_JACC_SEEDS = {
    "BaseQSketch": lambda m, seeds: BaseQSketch(m, seeds, 8),
    "FastQSketch": lambda m, seeds: FastQSketch(m, seeds, amount_bits=8),
    "QSketchDyn": lambda m, seeds: QSketchDyn(m, seeds, amount_bits=8, g_seed=42),
    "QSketch": lambda m, seeds: QSketch(m, seeds, amount_bits=8),
    "BaseLogExpSketch": lambda m, seeds: BaseLogExpSketch(m, seeds, amount_bits=8, logarithm_base=2),
    "FastLogExpSketch": lambda m, seeds: FastLogExpSketch(m, seeds, amount_bits=8, logarithm_base=2),
    "BaseShiftedLogExpSketch": lambda m, seeds: BaseShiftedLogExpSketch(m, seeds, amount_bits=8, logarithm_base=2),
    "FastShiftedLogExpSketch": lambda m, seeds: FastShiftedLogExpSketch(m, seeds, amount_bits=8, logarithm_base=2),
}

def _make_no_seeds_cls(sketch_clss: dict[str, Callable]) -> dict:
    ret = {}
    for sketch_name, sketch_cls in sketch_clss.items():
        ret[f"{sketch_name}_noseeds"] \
            = lambda m, seeds, _sketch_cls_lock=sketch_cls: _sketch_cls_lock(m, []) # this is busted
    return ret

_SKETCH_CLS_JACC_NO_SEEDS = _make_no_seeds_cls(_SKETCH_CLS_JACC_SEEDS)
_SKETCH_CLS_NO_JACC_NO_SEEDS = _make_no_seeds_cls(_SKETCH_CLS_NO_JACC_SEEDS)

_SKETCH_CLS_SEEDS = _SKETCH_CLS_JACC_SEEDS |_SKETCH_CLS_NO_JACC_SEEDS
_SKETCH_CLS_NO_SEEDS = _SKETCH_CLS_JACC_NO_SEEDS | _SKETCH_CLS_NO_JACC_NO_SEEDS
_SKETCH_CLS_JACC = _SKETCH_CLS_JACC_SEEDS | _SKETCH_CLS_JACC_NO_SEEDS
_SKETCH_CLS_ALL = _SKETCH_CLS_JACC_SEEDS | _SKETCH_CLS_JACC_NO_SEEDS |_SKETCH_CLS_NO_JACC_SEEDS | _SKETCH_CLS_NO_JACC_NO_SEEDS

