import subprocess as subp
import zipfile
from pathlib import Path
from tempfile import TemporaryDirectory
import shutil
import numpy

with TemporaryDirectory() as d:
    subp.run(["pip", "download", "numpy~=1.19.0"], cwd=d)
    d = Path(d)
    p = list(d.glob("numpy-1.19*"))[0]
    with zipfile.ZipFile(p) as f:
        f.extractall(d)
    p = list(d.glob("numpy*/numpy/f2py"))[0]
    p2 = Path(numpy.__file__).parent / "f2py"
    shutil.rmtree(p2)
    shutil.copytree(p, p2)
