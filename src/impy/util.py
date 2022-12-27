"""Utility module for auxiliary methods and classes."""

import warnings
import inspect
import platform
from pathlib import Path
import urllib.request
import zipfile
import shutil
import numpy as np
from typing import Sequence
from particle import Particle, PDGID, ParticleNotFound, InvalidParticle
from impy.constants import MeV, nucleon_mass


def energy2momentum(E, m):
    return np.sqrt((E + m) * (E - m))


def elab2ecm(elab, m1, m2):
    return np.sqrt(2.0 * elab * m2 + m2**2 + m1**2)


def mass(pdgid):
    m = Particle.from_pdgid(pdgid).mass
    if m is None:
        a = pdg2AZ(pdgid)[0]
        if a == 0:
            raise ValueError(f"cannot get mass for {pdgid}")
        return a * nucleon_mass
    return m * MeV


def pdg2name(pdgid):
    try:
        return Particle.from_pdgid(pdgid).name
    except ParticleNotFound:
        return f"Unknown({int(pdgid)})"
    except InvalidParticle:
        return f"Invalid({int(pdgid)})"


def is_AZ(arg):
    if not isinstance(arg, Sequence):
        return False
    if len(arg) != 2:
        return False
    for x in arg:
        if not isinstance(x, int):
            return False
    return True


def pdg2AZ(pdgid):
    """Returns mass number :math:`A`, charge :math:`Z` and neutron
    number :math:`N` of ``pdgid``.

    Note::

        PDG ID for nuclei is coded according to 10LZZZAAAI. For iron-52 it is 1000260520.

    Args:
        pdgid (int): PDG ID of nucleus/mass group
    Returns:
        (int, int): (A, Z) tuple
    """
    p = PDGID(pdgid)
    if p.is_nucleus:
        return p.A, p.Z
    elif pdgid == 2112:
        return 1, 0
    elif pdgid == 2212:
        return 1, 1
    return 0, 0


def AZ2pdg(A, Z):
    """Conversion of nucleus with mass A and charge Z
    to PDG nuclear code"""
    # 10LZZZAAAI
    pdg_id = 1000000000
    pdg_id += 10000 * Z
    pdg_id += 10 * A
    return pdg_id


def fortran_chars(array_ref, char_seq):
    """Helper to set fortran character arrays with python strings"""
    info(10, "Setting fortran array with", char_seq)
    # Reset
    import numpy as np

    len_arr = int(str(array_ref.dtype)[2:])
    len_seq = len(char_seq)
    return np.array(
        [c for c in char_seq + (len_arr - len_seq) * " "], dtype="S" + str(len_arr)
    )


def caller_name(skip=2):
    """Get a name of a caller in the format module.class.method

    `skip` specifies how many levels of stack to skip while getting caller
    name. skip=1 means "who calls me", skip=2 "who calls my caller" etc.
    An empty string is returned if skipped levels exceed stack height.abs

    From https://gist.github.com/techtonik/2151727
    """

    stack = inspect.stack()
    start = 0 + skip

    if len(stack) < start + 1:
        return ""

    parentframe = stack[start][0]

    name = []

    module = inspect.getmodule(parentframe)
    # `modname` can be None when frame is executed directly in console
    if module:
        name.append(module.__name__ + ".")

    # detect classname
    if "self" in parentframe.f_locals:
        # I don't know any way to detect call from the object method
        # there seems to be no way to detect static method call - it will
        # be just a function call

        name.append(parentframe.f_locals["self"].__class__.__name__ + "::")

    codename = parentframe.f_code.co_name
    if codename != "<module>":  # top level usually
        name.append(codename + "(): ")  # function or a method

    del parentframe
    return "".join(name)


def info(min_dbg_level, *args):
    """Print to console if min_dbg_level <= impy.debug_level.

    The function determines automatically the name of caller and appends
    the message to it. Message can be a tuple of strings or objects
    which can be converted to string using `str()`.

    Args:
        min_dbg_level (int): Minimum debug level in config for printing
        message (tuple): Any argument or list of arguments that casts to str
    """
    import impy

    if min_dbg_level <= impy.debug_level:
        print(caller_name(), *args)


def _download_file(outfile, url):
    """Download a file from 'url' to 'outfile'"""
    from rich.progress import (
        Progress,
        TextColumn,
        BarColumn,
        SpinnerColumn,
        DownloadColumn,
        TimeRemainingColumn,
    )

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
    with Progress(
        TextColumn("[progress.description]{task.description}"),
        BarColumn() if total_size else SpinnerColumn(),
        DownloadColumn(),
        TimeRemainingColumn(),
        transient=True,
    ) as bar:
        task_id = bar.add_task(f"Downloading {fname}", total=total_size)

        with open(outfile, "wb") as f:
            chunk = True
            while chunk:
                chunk = response.read(blocksize)
                f.write(chunk)
                nchunk = len(chunk)
                wrote += nchunk
                bar.advance(task_id, nchunk)

    if total_size and wrote != total_size:
        raise ConnectionError(f"{fname} has not been downloaded")


# Function to check and download dababase files on github
def _cached_data_dir(url):
    """Checks for existence of version file
    "model_name_vxxx.zip". Downloads and unpacks
    zip file from url in case the file is not found

    Args:
        url (str): url for zip file
    """

    base_dir = Path(__file__).parent.absolute() / "iamdata"
    base_dir.mkdir(parents=True, exist_ok=True)

    vname = Path(url).stem
    model_dir = base_dir / vname.split("_v")[0]
    version_file = model_dir / vname
    if not version_file.exists():
        zip_file = base_dir / Path(url).name
        temp_dir = Path(model_dir.parent / f".{model_dir.name}")
        shutil.rmtree(temp_dir, ignore_errors=True)
        if model_dir.exists():
            shutil.move(model_dir, temp_dir)
        _download_file(zip_file, url)
        if zipfile.is_zipfile(zip_file):
            with zipfile.ZipFile(zip_file, "r") as zf:
                zf.extractall(base_dir.as_posix())
            zip_file.unlink()
            shutil.rmtree(temp_dir, ignore_errors=True)

        version_glob = vname.split("_v")[0]
        for vfile in model_dir.glob(f"{version_glob}_v*"):
            vfile.unlink
        with open(version_file, "w") as vf:
            vf.write(url)
    return str(model_dir) + "/"


class TaggedFloat:
    """Floating point type that is distinct from an ordinary float.

    TaggedFloat is a base class to inherit from. We use it to declare
    that certain numbers have special meaning, e.g.::

        class KineticEnergy(TaggedFloat):
            pass

    can be used to declare a float-like class which stores kinetic energies,
    as opposed to another float-like class which stores the total energy or
    the momentum. Functions with an energy parameter can react differently
    on a KineticEnergy or TotalEnergy.

    A tagged float only allows arithmetic with its own type or
    with type-less numbers (int, float).
    """

    __slots__ = "_value"

    def __init__(self, val):
        self._value = float(val)

    def __repr__(self):
        return f"{self.__class__.__name__}({self._value!r})"

    def __float__(self):
        return self._value

    @classmethod
    def _reduce(cls, other):
        if isinstance(other, (int, float)):
            return other
        elif isinstance(other, cls):
            return other._value
        raise ValueError("{other!r} is not a number or {cls.__name__}")

    def __eq__(self, other):
        return self.__class__ is other.__class__ and self._value == other._value

    def __ne__(self, val):
        return not self == val

    def __mul__(self, val):
        return self.__class__(self._value * self._reduce(val))

    def __rmul__(self, val):
        return self * self._reduce(val)

    def __add__(self, val):
        return self.__class__(self._value + self._reduce(val))

    def __radd__(self, val):
        return self + self._reduce(val)

    def __truediv__(self, val):
        return self.__class__(self._value / self._reduce(val))

    def __rtruediv__(self, val):
        return self.__class__(self._reduce(val) / self._value)

    def __sub__(self, val):
        return self.__class__(self._value - self._reduce(val))

    def __rsub__(self, val):
        return self.__class__(self._reduce(val) - self._value)


# from Python-3.9 onwards, classmethod can be combined
# with property to replace this, which can then be removed
class classproperty:
    def __init__(self, f):
        self.f = f

    def __get__(self, obj, owner):
        return self.f(owner)


def _select_parents(mask, parents):
    # This algorithm is slow in pure Python and should be
    # speed up by compiling the logic.

    # attach parentless particles to beam particles,
    # unless those are also removed
    fallback = (0, 0)
    if mask[0] and mask[1]:
        fallback = (1, 2)

    n = len(parents)
    indices = np.arange(n)[mask] + 1
    result = parents[mask]
    mapping = {old: i + 1 for i, old in enumerate(indices)}

    n = len(result)
    for i in range(n):
        a = result[i, 0]
        if a == 0:
            continue
        p = mapping.get(a, -1)
        if p == -1:
            a, b = fallback
            result[i, 0] = a
            result[i, 1] = b
        elif p != a:
            q = 0
            b = result[i, 1]
            if b > 0:
                q = mapping.get(b, 0)
            result[i, 0] = p
            result[i, 1] = q
    return result


def select_parents(arg, parents):
    if parents is None:
        return None

    n = len(parents)

    if isinstance(arg, np.ndarray) and arg.dtype is bool:
        mask = arg
    else:
        mask = np.zeros(n, dtype=bool)
        mask[arg] = True

    with warnings.catch_warnings():
        # suppress numba safety warning that we can ignore
        warnings.simplefilter("ignore")
        return _select_parents(mask, parents)


try:
    # accelerate with numba if numba is available
    import numba as nb

    _select_parents = nb.njit(_select_parents)

except ModuleNotFoundError:
    pass


def tolerant_string_match(a, b):
    """
    Return True if all characters in appear also in b in same order.

    This algorithm is slow and should only be used were speed does
    not matter.
    """
    last = 0
    for c in a:
        i = b.find(c)
        if i == -1:
            return False
        if i < last:
            return False
        last = i
    return True


def get_all_models(skip=None):
    from impy import models
    from impy.common import MCRun

    if skip is None:
        skip = []
    if platform.system() == "Windows":
        skip = list(skip) + [models.UrQMD34, models.Pythia8]

    result = []
    for key in dir(models):
        obj = getattr(models, key)
        if skip and obj in skip:
            continue
        try:
            if issubclass(obj, MCRun):  # fails if obj is not a class
                result.append(obj)
        except TypeError:
            pass

    return result


class Singleton(type):
    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]


def naneq(a, b):
    """
    Return True if a == b or if a and b are both NaN.

    Parameters
    ----------
    a : float
        First float.
    b : float
        Second float.
    """
    return a == b or (np.isnan(a) and np.isnan(b))


def fortran_array_insert(array, size, value):
    if value in array[:size]:
        return
    if len(array) == size:
        raise RuntimeError("array is full")
    array[size] = value
    size += 1


def fortran_array_remove(array, size, value):
    for i, val in enumerate(array[:size]):
        if val == value:
            size -= 1
            for j in range(i, size):
                array[j] = array[j + 1]
            break
