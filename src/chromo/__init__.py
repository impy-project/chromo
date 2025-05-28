from chromo import models
from chromo import kinematics
from chromo import constants

import os
from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("chromo")
except PackageNotFoundError:
    __version__ = "0.0.0"

debug_level = int(os.environ.get("DEBUG", "0"))

__all__ = ["models", "kinematics", "constants", "debug_level", "__version__"]
