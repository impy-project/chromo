from .version import __version__
from pathlib import Path
from . import config as _config
from particletools.tables import PYTHIAParticleData

impy_config = _config.__dict__
base_path = Path(__file__).parent.absolute()

pdata = PYTHIAParticleData(
    cache_file=open(base_path / impy_config["pdata_cachefile"], "wb")
)

# these must be after pdata for now
# FIXME code should not depend on import order (HD)
from . import models  # noqa
from . import kinematics  # noqa
from . import constants  # noqa

__all__ = ["models", "kinematics", "constants", "__version__"]
