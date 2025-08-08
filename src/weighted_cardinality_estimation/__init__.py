from importlib import metadata as _md
from ._core import ExpSketch # type: ignore
from ._core import FastExpSketch # type: ignore

__all__ = ["ExpSketch", "FastExpSketch", "__version__"]
__version__ = _md.version(__name__)
