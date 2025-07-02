try:
    from ._version import __version__
except ImportError:  # pragma: no cover
    __version__ = "None"

from .FrequencySystemBuilder import FrequencySystemBuilder
from .TemporalSystemBuilder import TemporalSystemBuilder