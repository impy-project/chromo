import sys
from pathlib import Path
import warnings
import numpy
from packaging.version import parse as parse_version

# This script substitute patched f2py files from numpy==1.26
# Works only for numpy >= 1.23
if (parse_version(numpy.__version__) >= parse_version("1.23.0")) and (
    parse_version(numpy.__version__) < parse_version("1.26.1")
):
    # import _download_file
    sys.path.append(str(Path(__file__).parent / "src"))
    from chromo.util import _download_file  # noqa

    # Files to substitute
    files_to_patch = ["f2py/crackfortran.py", "f2py/auxfuncs.py"]
    numpy_dir = Path(numpy.__file__).parent
    # Commit of numpy to copy from
    # Change the url to the url of official release when it will be ready!!!
    base_url = "https://raw.githubusercontent.com/numpy/numpy/1e6a322a9514c0050f4cc656aead1aeffba0b1b5/numpy"

    warnings.warn(f"Files in numpy: {files_to_patch} are substituted")
    for fp in files_to_patch:
        _download_file(numpy_dir / fp, f"{base_url}/{fp}")
elif (parse_version(numpy.__version__) < parse_version("1.23.0")) and (
    parse_version(numpy.__version__) >= parse_version("1.21.0")
):
    warnings.warn(
        f"numpy=={numpy.__version__} f2py generates corrupted files!!! "
        f"Chromo will not work."
    )
