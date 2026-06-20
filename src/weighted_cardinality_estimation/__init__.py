from importlib import metadata as _md

from . import stat  # noqa: F401
from ._core import *  # noqa: F403  # pyright: ignore[reportMissingImports]
from ._core import (  # noqa: F401  # pyright: ignore[reportMissingImports]
    FastExpSketchCustomFloat,
    LogExpSketchFastNoShifted,
    LogExpSketchSlowNoShifted,
    MemoryFlag,
    QuantizationMode,
)

__version__ = _md.version(__name__)
