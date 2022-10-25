from .version import __version__  # noqa
from os.path import join, abspath, dirname
from . import config as _config

impy_config = _config.__dict__

from particletools.tables import PYTHIAParticleData

base_path = abspath(dirname(__file__))


# This is not nice, but the paths in the config should become absolute
# in case impy is used outside of the folder
for dpmmod in ["dpmjetIII", "phojet"]:
    for version_key in impy_config[dpmmod]["param_file"]:
        impy_config[dpmmod]["param_file"][version_key] = join(
            base_path, impy_config[dpmmod]["param_file"][version_key]
        )
        impy_config[dpmmod]["dat_dir"][version_key] = join(
            base_path, impy_config[dpmmod]["dat_dir"][version_key]
        )
        if dpmmod == "phojet":
            continue
        impy_config[dpmmod]["evap_file"][version_key] = join(
            base_path, impy_config[dpmmod]["evap_file"][version_key]
        )

impy_config["epos"]["datdir"] = join(base_path, impy_config["epos"]["datdir"])

pdata = PYTHIAParticleData(
    cache_file=open(join(base_path, impy_config["pdata_cachefile"]), "wb")
)

import pathlib
from .util import PathFromZip

base_url = "https://github.com/afedynitch/MCEq/releases/download/"
release_tag = "builds_on_azure/"
zip_fname = "zipped_files.zip"
remote_zip_file = base_url + release_tag + zip_fname


_impy_base_dir = pathlib.Path(base_path)
# Creating instance of a class for lazy downloading
_get_impy_data_path = PathFromZip(
    local_zip=_impy_base_dir / "impy_iamdata.zip",
    remote_zip=remote_zip_file,
    cache_downloaded=_impy_base_dir / "cache_downloaded_files.cch",
)

# Checks the full_path (relative to _impy_base_dir)
def _check_impy_data_path(full_path):
    return _get_impy_data_path(
        pathlib.Path(full_path).relative_to(_impy_base_dir),
        _impy_base_dir,
    ).as_posix()


import impy.models as models  # noqa
import impy.kinematics as kinematics  # noqa
import impy.constants as constants  # noqa
