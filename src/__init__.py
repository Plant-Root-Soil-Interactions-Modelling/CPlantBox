from importlib import resources as _resources

from . import plantbox as _core
from .plantbox import *


def data_path():
    """Return the installed model parameter data directory."""
    return str(_resources.files(__package__) / "modelparameter")
