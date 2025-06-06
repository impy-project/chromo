from pathlib import Path
import warnings
import numpy
import urllib.request
from packaging.version import parse as parse_version


# Version from util without rich
def _download_file(outfile, url):
    """Download a file from 'url' to 'outfile'"""
    fname = Path(url).name
    try:
        response = urllib.request.urlopen(url)
    except BaseException:
        raise ConnectionError(
            f"_download_file: probably something wrong with url = '{url}'"
        )
    total_size = response.getheader("content-length")

    min_blocksize = 4096
    if total_size:
        total_size = int(total_size)
        blocksize = max(min_blocksize, total_size // 100)
    else:
        blocksize = min_blocksize

    wrote = 0
    with open(outfile, "wb") as f:
        chunk = True
        while chunk:
            chunk = response.read(blocksize)
            f.write(chunk)
            nchunk = len(chunk)
            wrote += nchunk
            if total_size:
                print(
                    f"Downloading {fname}: {wrote/total_size*100:.0f}% "
                    f"done ({wrote/(1024*1024):.0f} Mb) \r",
                    end=" ",
                )
            else:
                print(
                    f"Downloading {fname}: {wrote/(1024*1024):.0f} Mb downloaded \r",
                    end=" ",
                )
    print()
    if total_size and wrote != total_size:
        raise ConnectionError(f"{fname} has not been downloaded")


# This script substitute patched f2py files from numpy==1.26
# Works only for numpy >= 1.23
if (parse_version(numpy.__version__) >= parse_version("1.23.0")) and (
    parse_version(numpy.__version__) < parse_version("1.26.1")
):
    # Files to substitute
    files_to_patch = ["f2py/crackfortran.py", "f2py/auxfuncs.py"]
    numpy_dir = Path(numpy.__file__).parent
    # Commit of numpy to copy from
    base_url = "https://raw.githubusercontent.com/numpy/numpy/maintenance/1.26.x/numpy"

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
