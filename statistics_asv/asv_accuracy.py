import numpy as np

from weighted_cardinality_estimation.stat import elements_stream

from .common import IMPLS, get_seeds

# DRAFT: small values so the suite runs quickly.
# Bump STATISTICAL_RUNS to >=1000 before publishing results.
SKETCH_SIZE = 64
STREAM_SIZE = 200
STATISTICAL_RUNS = 50


class AccuracySuite:
    param_names = ["sketch_type"]
    params = [list(IMPLS.keys())]

    def setup_cache(self):
        elems = elements_stream(STREAM_SIZE, seed=0)
        weights = [1.0] * STREAM_SIZE
        true_cardinality = float(STREAM_SIZE)

        stats: dict[str, dict] = {}
        for impl_name in IMPLS:
            estimates = []
            for _ in range(STATISTICAL_RUNS):
                s = IMPLS[impl_name](SKETCH_SIZE, get_seeds(SKETCH_SIZE))
                s.add_many(elems, weights)
                estimates.append(s.estimate())

            arr = np.array(estimates)
            mean = float(np.mean(arr))
            mre = float(np.mean(np.abs(arr - true_cardinality) / true_cardinality))
            cv = float(np.std(arr) / mean) if mean != 0.0 else 0.0
            stats[impl_name] = {"mean": mean, "mre": mre, "cv": cv}

        return stats

    def track_mean_estimate(self, stats: dict, impl_name: str) -> float:
        return stats[impl_name]["mean"]

    track_mean_estimate.unit = "units"  # type: ignore

    def track_relative_error(self, stats: dict, impl_name: str) -> float:
        """Mean relative error in %"""
        return stats[impl_name]["mre"] * 100

    track_relative_error.unit = "%"  # type: ignore

    def track_coeff_of_variation(self, stats: dict, impl_name: str) -> float:
        """Coefficient of variation in %"""
        return stats[impl_name]["cv"] * 100

    track_coeff_of_variation.unit = "%"  # type: ignore
