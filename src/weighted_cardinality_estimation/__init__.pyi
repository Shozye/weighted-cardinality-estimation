# src/weighted_cardinality_estimation/__init__.pyi
from . import stat as stat


class Sketch:
    """Every sketch in Weighted Cardinality Estimation inherits it"""

    def add(self, x: str, weight: float = ...) -> None: 
        """Used to add a single element to the sketch"""
        ...
    def add_many(self, elems: list[str], weights: list[float]) -> None: 
        """Used to add a list of elements in batch. Uses C++ loop to increase speed"""
        ...
    def estimate(self) -> float: 
        """Main method to estimation using this Sketch. It is the best one."""
        ...
    def memory_usage_total(self) -> int: 
        """Method to get information about total memory usage during structure working."""
        ...
    def memory_usage_write(self) -> int: 
        """Method to get non-read-only memory usage. It assumes that seeds/permutation etc is 'given by system' """
        ...
    def memory_usage_estimate(self) -> int: 
        """Memory usage of stuff only needed to estimation. So it takes all structures etc used in estimate. """
        ...

class ExpSketch(Sketch):
    """Exponential sketch; Lemiesz Sketch. Exp on Double registers"""

    def __init__(self, m: int, seeds: list[int]) -> None: ...
    def jaccard_struct(self, other: "ExpSketch") -> float: ...

class ExpSketchFloat(Sketch):
    """The same as ExpSketch but uses float registers in the structure"""

    def __init__(self, m: int, seeds: list[int]) -> None: ...
    def jaccard_struct(self, other: "ExpSketchFloat") -> float: ...

class FastExpSketch(Sketch):
    """ExpSketch with a speed up framework made by dr. Lemiesz """

    def __init__(self, m: int, seeds: list[int]) -> None: ...
    def jaccard_struct(self, other: "FastExpSketch") -> float: ...

class FastGMExpSketch(Sketch):
    """Alternative to speed up framework. This one is made by other people."""

    def __init__(self, m: int, seeds: list[int]) -> None: ...
    def jaccard_struct(self, other: "FastGMExpSketch") -> float: ...

class BaseQSketch(Sketch):
    """This is a QSketch without speedup framework. It only has got the structure inside and the same estimator as in the QSketch Paper"""

    def __init__(self, m: int, seeds: list[int], amount_bits: int) -> None: ...

class FastQSketch(Sketch):
    """This is a BaseQSketch but sped up by FastExpSketch framework made by dr. Lemiesz"""

    def __init__(self, m: int, seeds: list[int], amount_bits: int) -> None: ...

class QSketch(Sketch):
    """This is a QSketch exactly as defined in the QSketch Paper"""

    def __init__(self, m: int, seeds: list[int], amount_bits: int) -> None: ...

class QSketchDyn(Sketch):
    """Sketch from the QSketch paper that uses histograms etc."""

    def __init__(self, m: int, seeds: list[int], amount_bits: int, g_seed: int = ...) -> None: ...

class BaseLogExpSketch(Sketch):
    """Exp sketch with logarithmic register quantisation. Doesn't have got any speed ups nor jaccard."""

    newton_iteration_count: int

    def __init__(self, m: int, seeds: list[int], amount_bits: int, logarithm_base: float) -> None: ...
    def estimate_direct(self) -> float: ...
    def estimate_newton_cold(self) -> float: ...
    def estimate_newton_warm(self) -> float: ...

class BaseLogExpSketchJacc(Sketch):
    """BaseLogExpSketch with H-Structure for Jaccard Estimation"""

    def __init__(self, m: int, seeds: list[int], amount_bits: int, logarithm_base: float, amount_bits_jaccard: int) -> None: ...
    def jaccard_struct(self, other: "BaseLogExpSketchJacc") -> float: ...

class FastLogExpSketch(Sketch):
    """Fast variant of BaseLogExpSketch. It uses 'FastExpSketch' framework"""

    def __init__(self, m: int, seeds: list[int], amount_bits: int, logarithm_base: float) -> None: ...

class BaseShiftedLogExpSketch(Sketch):
    """BaseLogExpSketch with shifted registers to save up on memory"""

    def __init__(self, m: int, seeds: list[int], amount_bits: int, logarithm_base: float) -> None: ...

class FastShiftedLogExpSketch(Sketch):
    """Combination of FastLogExpSketch and BaseShiftedLogExpSketch"""

    def __init__(self, m: int, seeds: list[int], amount_bits: int, logarithm_base: float) -> None: ...

class WeightedMinHash(Sketch):
    """Weighted MinHash based on Cohen, Katzir & Yehezkel (IPL 2015).
    Generalises the max-sketch cardinality estimator to weighted streams via
    the Beta(w,1) transform: h(x)^(1/w). Estimator: m*(m-1)/sum(1-M[k])."""

    def __init__(self, m: int, seeds: list[int]) -> None: ...
