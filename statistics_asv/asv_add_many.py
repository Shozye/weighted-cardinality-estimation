from weighted_cardinality_estimation.stat import elements_stream

from .common import IMPLS, get_seeds

# DRAFT: small values so the suite runs quickly.
SKETCH_SIZE = 64
AMOUNT_ELEMENTS = 500


class AddManySuite:
    param_names = ["sketch_type"]
    params = [list(IMPLS.keys())]

    def setup(self, impl_name: str):
        self.instance = IMPLS[impl_name](SKETCH_SIZE, get_seeds(SKETCH_SIZE))
        self.elems = elements_stream(AMOUNT_ELEMENTS, seed=0)
        self.weights = [1.0] * AMOUNT_ELEMENTS

    def time_add_many(self, impl_name: str):
        self.instance.add_many(self.elems, self.weights)

    time_add_many.rounds = 2  # type: ignore
    time_add_many.repeat = 3  # type: ignore
    time_add_many.warmup_time = 0.1  # type: ignore
