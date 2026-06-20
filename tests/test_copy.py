"""Copy and serialization: __getstate__/__setstate__, copy, pickle, deepcopy."""

import copy
import pickle


class TestCopy:
    def test_copy_same_estimate(self, sketch) -> None:
        sketch.add("elem")
        copied = copy.copy(sketch)
        assert copied.estimate() == sketch.estimate()
        assert copied is not sketch

    def test_copy_identical_state(self, sketch) -> None:
        sketch.add("elem")
        copied = copy.copy(sketch)
        assert sketch.__getstate__() == copied.__getstate__()

    def test_copy_independent_mutation(self, sketch) -> None:
        sketch.add("elem")
        estimate_before = sketch.estimate()
        copied = copy.copy(sketch)
        copied.add("new_elem")
        assert sketch.estimate() == estimate_before

    def test_copy_same_behavior_after_mutation(self, sketch) -> None:
        sketch.add("elem")
        copied = copy.copy(sketch)
        sketch.add("new_elem")
        copied.add("new_elem")
        assert sketch.estimate() == copied.estimate()

    def test_pickle_roundtrip(self, sketch) -> None:
        sketch.add("elem")
        data = pickle.dumps(sketch)
        restored = pickle.loads(data)
        assert restored.estimate() == sketch.estimate()
        assert restored.__getstate__() == sketch.__getstate__()

    def test_deepcopy_roundtrip(self, sketch) -> None:
        sketch.add("elem")
        deep = copy.deepcopy(sketch)
        assert deep.estimate() == sketch.estimate()
        assert deep.__getstate__() == sketch.__getstate__()
