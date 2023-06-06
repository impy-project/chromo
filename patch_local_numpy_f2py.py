import sys
import zipfile
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory
import os
import numpy

fn = (Path() / "numpy-1.19.5.zip").absolute()

if not Path(fn).exists():
    sys.path.append("src/chromo")

    from chromo.util import _download_file  # noqa

    url = "https://github.com/numpy/numpy/releases/download/v1.19.5/numpy-1.19.5.zip"

    _download_file(fn, url)

with TemporaryDirectory() as d:
    os.chdir(d)
    with zipfile.ZipFile(fn) as f:
        for x in f.infolist():
            if "numpy/f2py" in x.filename:
                f.extract(x)

    numpy_dir = Path(numpy.__file__).parent
    shutil.rmtree(numpy_dir / "f2py")
    shutil.move("numpy-1.19.5/numpy/f2py", numpy_dir)
