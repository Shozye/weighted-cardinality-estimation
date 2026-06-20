from .common import IMPLS

# DRAFT: deterministic seeds so the tracked value is reproducible across commits.
SKETCH_SIZE = 64


class RegressionSuite:
    """Tracks the estimate value of each sketch with a fixed, deterministic seed.

    Purpose: catch if a code change accidentally alters the mathematical output.
    If this value changes between commits, something in the algorithm changed.
    """

    param_names = ["sketch_type"]
    params = [list(IMPLS.keys())]

    def setup(self, impl_name: str):
        seeds = list(range(1, SKETCH_SIZE + 1))
        self.instance = IMPLS[impl_name](SKETCH_SIZE, seeds)
        self.instance.add("fixed_element", 1.0)

    def track_estimate(self, impl_name: str) -> float:
        return self.instance.estimate()

    track_estimate.unit = "units"  # type: ignore
