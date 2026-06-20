import pickle
import zlib

from weighted_cardinality_estimation import MemoryFlag

from .common import IMPLS, get_seeds

# DRAFT: small values so the suite runs quickly.
SKETCH_SIZE = 64


class MemorySuite:
    param_names = ["sketch_type"]
    params = [list(IMPLS.keys())]

    def setup(self, impl_name: str):
        self.instance = IMPLS[impl_name](SKETCH_SIZE, get_seeds(SKETCH_SIZE))
        self.instance.add("element", 1.0)

    def track_total_memory(self, impl_name: str) -> int:
        return self.instance.memory_usage(MemoryFlag.TOTAL)

    track_total_memory.unit = "bytes"  # type: ignore

    def track_write_memory(self, impl_name: str) -> int:
        return self.instance.memory_usage(MemoryFlag.ALL_WRITE)

    track_write_memory.unit = "bytes"  # type: ignore

    def track_estimate_memory(self, impl_name: str) -> int:
        return self.instance.memory_usage(MemoryFlag.REGISTERS)

    track_estimate_memory.unit = "bytes"  # type: ignore

    def track_serialization_size(self, impl_name: str) -> int:
        return len(zlib.compress(pickle.dumps(self.instance)))

    track_serialization_size.unit = "bytes"  # type: ignore
