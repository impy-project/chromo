from impy.version import __version__
from impy import config as _config

impy_config = _config.__dict__
# these must be after impy_config for now
# FIXME code should not depend on import order (HD)
from impy import models  # noqa
from impy import kinematics  # noqa
from impy import constants  # noqa

__all__ = ["models", "kinematics", "constants", "__version__"]
