from .version import __version__  # noqa
from os.path import join, abspath, dirname
from . import config as _config

impy_config = _config.__dict__

from particletools.tables import PYTHIAParticleData

base_path = abspath(dirname(__file__))

pdata = PYTHIAParticleData(
    cache_file=open(join(base_path, impy_config["pdata_cachefile"]), "wb")
)


import impy.models as models  # noqa
import impy.kinematics as kinematics  # noqa
import impy.constants as constants  # noqa
