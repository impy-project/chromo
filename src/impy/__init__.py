from impy.version import __version__
from impy import models
from impy import kinematics
from impy import constants
import os

debug_level = int(os.environ.get("DEBUG", "0"))

__all__ = ["models", "kinematics", "constants", "debug_level", "__version__"]
