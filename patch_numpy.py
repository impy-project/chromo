import subprocess as subp
import zipfile
from pathlib import Path
from tempfile import TemporaryDirectory
import shutil
import numpy

with TemporaryDirectory() as d:
    subp.run(["pip", "download", "numpy~=1.19.0"], cwd=d)
    p = list(Path(d).glob("numpy-1.19*"))[0]
    with zipfile.ZipFile(p) as f:
        p = [x for x in f.namelist() if x.endswith("numpy/f2py/")][0]
        p = f.extract(p)
    p2 = Path(numpy.__file__).parent / "f2py"
    shutil.rmtree(p2)
    shutil.copytree(p, p2)
