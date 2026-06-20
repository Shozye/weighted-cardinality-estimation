"""Newton estimator tests for kQSketch variants."""

from weighted_cardinality_estimation.stat import relative_error, weighted_stream


class TestNewton:
    M = 200
    WEIGHT = 1e4

    def _make_sketch(self, newton_spec):
        s = newton_spec.factory(self.M, 42)
        elems, weights = weighted_stream(500, self.WEIGHT)
        s.add_many(elems, weights)
        return s

    def test_three_estimators_callable(self, newton_spec) -> None:
        s = self._make_sketch(newton_spec)
        assert s.estimate_direct() > 0
        assert s.estimate_newton_cold() > 0
        assert s.estimate_newton_warm() > 0

    def test_direct_and_warm_agree(self, newton_spec) -> None:
        s = self._make_sketch(newton_spec)
        assert relative_error(s.estimate_direct(), s.estimate_newton_warm()) < 0.02

    def test_cold_and_warm_same_ballpark(self, newton_spec) -> None:
        s = self._make_sketch(newton_spec)
        assert relative_error(s.estimate_newton_cold(), s.estimate_newton_warm()) < 0.05

    def test_warm_iterations_lte_cold_iterations(self, newton_spec) -> None:
        s = self._make_sketch(newton_spec)
        assert s.estimate_newton_warm_iterations() <= s.estimate_newton_cold_iterations()
