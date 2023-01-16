from impy import models
from impy import kinematics
from impy import constants
import os
from importlib.metadata import version

__version__ = version("impy")

debug_level = int(os.environ.get("DEBUG", "0"))

__all__ = ["models", "kinematics", "constants", "debug_level", "__version__"]
