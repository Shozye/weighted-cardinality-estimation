from .common import IMPLS, get_seeds

# DRAFT: small values so the suite runs quickly.
SKETCH_SIZE = 64


class EstimateSuite:
    param_names = ["sketch_type"]
    params = [list(IMPLS.keys())]

    def setup(self, impl_name: str):
        self.instance = IMPLS[impl_name](SKETCH_SIZE, get_seeds(SKETCH_SIZE))
        self.instance.add("element", 1.0)

    def time_estimate(self, impl_name: str):
        self.instance.estimate()

    time_estimate.rounds = 2  # type: ignore
    time_estimate.repeat = 3  # type: ignore
    time_estimate.warmup_time = 0.1  # type: ignore
