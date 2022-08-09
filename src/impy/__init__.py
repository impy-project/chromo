import sys
from os.path import join, abspath, dirname
import yaml
from particletools.tables import PYTHIAParticleData


base_path = abspath(dirname(__file__))
sys.path.append(join(base_path, "lib"))


# Globals
impy_config = yaml.load(
    open(join(base_path, "impy_config.yaml")), Loader=yaml.FullLoader
)

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


import impy.models as models  # noqa
import impy.kinematics as kinematics  # noqa
import impy.constants as constants  # noqa
