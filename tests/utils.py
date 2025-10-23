from weighted_cardinality_estimation import BaseLogExpSketch, BaseLogExpSketchJacc, BaseQSketch, FastExpSketch, ExpSketch, FastGMExpSketch, BaseLogExpSketch, FastLogExpSketch, FastQSketch, QSketchDyn, BaseShiftedLogExpSketch, FastShiftedLogExpSketch, QSketch

def assert_error(expected: float, actual: float, error: float):
    min_range = expected * (1-error)
    max_range = expected * (1+error)
    assert min_range < actual < max_range, f"{min_range=} < {actual=} < {max_range}"

_SKETCH_CONSTRUCTORS = {
    "ExpSketch": lambda m, seeds: ExpSketch(m, seeds),
    "FastExpSketch": lambda m, seeds: FastExpSketch(m, seeds),
    "FastGMExpSketch": lambda m, seeds: FastGMExpSketch(m, seeds),
    "BaseQSketch": lambda m, seeds: BaseQSketch(m, seeds, 8),
    "FastQSketch": lambda m, seeds: FastQSketch(m, seeds, amount_bits=8),
    "QSketchDyn": lambda m, seeds: QSketchDyn(m, seeds, amount_bits=8, g_seed=42),
    "QSketch": lambda m, seeds: QSketch(m, seeds, amount_bits=8),
    "BaseLogExpSketch": lambda m, seeds: BaseLogExpSketch(m, seeds, amount_bits=8, logarithm_base=2),
    "BaseLogExpSketchJacc": lambda m, seeds: BaseLogExpSketchJacc(m, seeds, amount_bits=8, logarithm_base=2, amount_bits_jaccard=8),
    "FastLogExpSketch": lambda m, seeds: FastLogExpSketch(m, seeds, amount_bits=8, logarithm_base=2),
    "BaseShiftedLogExpSketch": lambda m, seeds: BaseShiftedLogExpSketch(m, seeds, amount_bits=8, logarithm_base=2),
    "FastShiftedLogExpSketch": lambda m, seeds: FastShiftedLogExpSketch(m, seeds, amount_bits=8, logarithm_base=2),
}