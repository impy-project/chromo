from impy.version import __version__
from impy import models
from impy import kinematics
from impy import constants
import os
from impy.util import AliveInstanceWarning, FixedStableListWarning

debug_level = int(os.environ.get("DEBUG", "0"))

__all__ = [
    "models",
    "kinematics",
    "constants",
    "debug_level",
    "AliveInstanceWarning",
    "FixedStableListWarning",
    "__version__",
]
