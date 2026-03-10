import os
from importlib.metadata import version

from chromo import constants, kinematics, models

__version__ = version("chromo")

debug_level = int(os.environ.get("DEBUG", "0"))

__all__ = ["__version__", "constants", "debug_level", "kinematics", "models"]
