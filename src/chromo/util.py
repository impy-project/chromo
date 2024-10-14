"""Utility module for auxiliary methods and classes."""

import warnings
import inspect
import platform
import dataclasses
import copy
import math
import zipfile
import shutil
from pathlib import Path
from enum import Enum
from typing import Sequence, Set, Tuple, Collection, Union
import urllib.request
import numpy as np
from particle import Particle, PDGID, ParticleNotFound, InvalidParticle
from chromo.constants import MeV, nucleon_mass, sec2cm

EventFrame = Enum("EventFrame", ["CENTER_OF_MASS", "FIXED_TARGET", "GENERIC"])


@dataclasses.dataclass(init=False)
class CompositeTarget:
    """Definition of composite targets made of multiple (atomic) nuclei.

    Examples of such composite targets are Air, CO_2, HCl, C_2H_60.
    """

    label: str
    components: Tuple[PDGID]
    fractions: np.ndarray

    def __init__(
        self, components: Collection[Tuple[Union[str, int], float]], label: str = ""
    ):
        """
        Parameters
        ----------
        components : collection of (str|int, float)
            The components of the targets. Each component is given by a string or PDGID
            that identifies the element, and its relative amount in the material.
            Amounts do not have to add up to 1, fractions are computed automatically.
        label : str, optional
            Give the target a name. This is purely cosmetic.
        """

        if len(components) == 0:
            raise ValueError("components cannot be empty")
        fractions = np.empty(len(components))
        c = []
        for i, (particle, amount) in enumerate(components):
            fractions[i] = amount
            p = process_particle(particle)
            if not p.is_nucleus:
                raise ValueError(f"component {particle} is not a nucleus")
            c.append(p)
        self.label = label
        self.components = tuple(c)
        self.fractions = fractions / np.sum(fractions)
        self.fractions.flags["WRITEABLE"] = False

    def copy(self):
        new_target = CompositeTarget([("N", 1)])
        for field in dataclasses.fields(CompositeTarget):
            setattr(new_target, field.name, copy.copy(getattr(self, field.name)))
        return new_target

    def __eq__(self, other):
        if not isinstance(other, CompositeTarget):
            return False
        return (
            (self.label == other.label)
            & (self.components == other.components)
            & (np.allclose(self.fractions, other.fractions))
        )

    @property
    def Z(self):
        """Return maximum charge number."""
        # needed for compatibility with PDGID interface and for dpmjet initialization
        return max(p.Z for p in self.components)

    @property
    def A(self):
        """Return maximum number of nucleons."""
        # needed for compatibility with PDGID interface and for dpmjet initialization
        return max(p.A for p in self.components)

    @property
    def is_nucleus(self):
        return True

    @property
    def is_hadron(self):
        return False

    def __int__(self):
        """Return PDGID for heaviest of elements."""
        return int(max((c.A, c) for c in self.components)[1])

    # this allows us to use CompositeTarget in Set[int].__contains__
    def __hash__(self):
        return self.__int__().__hash__()

    def average_mass(self):
        return sum(
            f * p.A * nucleon_mass for (f, p) in zip(self.fractions, self.components)
        )

    def __abs__(self):
        return abs(self.__int__())

    def __repr__(self):
        components = [
            (pdg2name(c), float(amount))
            for (c, amount) in zip(self.components, self.fractions)
        ]
        args = f"{components}"
        if self.label:
            args += f", label={self.label!r}"
        return f"CompositeTarget({args})"


def is_real_nucleus(pdgid: Union[int, PDGID, CompositeTarget]) -> bool:
    """
    Return True if pdgid is a nucleus with A > 1.

    PDGID.is_nucleus is True also for proton and neutrons,
    which is correct in some sense, but often we want to
    handle only nuclei with A > 1 in the interface.

    Also works for CompositeTarget.
    """
    if not isinstance(pdgid, PDGID):
        pdgid = PDGID(pdgid)
    return pdgid.A and pdgid.A > 1


def energy2momentum(E, m):
    """
    Compute the momentum of a particle given its energy and mass.

    Numerically more stable way to compute E^2 - m^2.

    Args:
        E (float): The energy of the particle.
        m (float): The mass of the particle.

    Returns:
        float: The momentum of the particle.
    """
    return np.sqrt((E + m) * (E - m))


def momentum2energy(p, m):
    # only a minor trick can be used here, add in order
    # of increasing scale
    a, b = (p, m) if p < m else (m, p)
    return np.sqrt(a**2 + b**2)


def elab2ecm(elab, m1, m2):
    # ecm^2 = s = ((p1^2 + m1^2)^0.5 + (p2^2 + m2^2)^0.5)^2 - (p1 + p2)^2
    #   with   elab = (p1^2 + m1^2)^0.5,    p2 = 0
    #       = (elab + m2)^2 - p1^2
    #       = (elab + m2)^2 - (elab^2 - m1^2)
    #       = elab^2 + 2 elab m2 + m2^2 - elab^2 + m1^2
    #       = 2 elab m2 + m1^2 + m2^2
    # sum in order of increasing size to improve numerical accuracy
    return np.sqrt(m1**2 + m2**2 + 2.0 * elab * m2)


def ecm2elab(ecm, m1, m2):
    """
    Calculates the center-of-mass energy (ECM) in the lab frame given
    the masses of two particles (m1, m2).

    Args:
        ecm (float): The center-of-mass energy.
        m1 (float): The mass of particle 1.
        m2 (float): The mass of particle 2.

    Returns:
        float: The center-of-mass energy in the lab frame.
    """
    return 0.5 * (ecm**2 - m1**2 - m2**2) / m2


def mass(pdgid):
    """
    Returns the mass of a particle with the given PDG ID in MeV/c^2.

    Args:
        pdgid (int): The PDG ID of the particle.

    Returns:
        float: The mass of the particle in MeV/c^2.

    Raises:
        ValueError: If the mass cannot be determined for the given PDG ID.
    """
    m = Particle.from_pdgid(pdgid).mass
    if m is None:
        a = pdg2AZ(pdgid)[0]
        if a == 0:
            raise ValueError(f"cannot get mass for {pdgid}")
        return a * nucleon_mass
    return m * MeV


def _make_name2pdg_db():
    all_particles = Particle.findall()
    db = {p.name: p.pdgid for p in all_particles}
    db.update({p.programmatic_name: p.pdgid for p in all_particles})
    db["p"] = PDGID(2212)
    db["n"] = PDGID(2112)
    db["p~"] = -db["p"]
    db["n~"] = -db["n"]
    db.update(
        H=db["p"],
        H1=db["p"],
        He=db["He4"],
        C=db["C12"],
        N=db["N14"],
        O=db["O16"],
        Ne=db["Ne20"],
        Ar=db["Ar40"],
        Xe=db["Xe131"],
        Pb=db["Pb206"],
        photon=db["gamma"],
        proton=db["p"],
        neutron=db["n"],
        antiproton=-db["p"],
        antineutron=-db["n"],
        pbar=-db["p"],
        nbar=-db["n"],
        p_bar=-db["p"],
        n_bar=-db["n"],
    )
    return db


_name2pdg_db = _make_name2pdg_db()


def name2pdg(name: str):
    """
    Given a particle name, returns the corresponding PDG code.

    Args:
        name (str): The name of the particle.

    Returns:
        int: The PDG code of the particle.
    """
    return _name2pdg_db[name]


def pdg2name(pdgid):
    """
    Given a particle's PDG ID, returns its name.

    Args:
        pdgid (int): The particle's PDG ID.

    Returns:
        str: The particle's name, or "Unknown" or "Invalid" if the
        PDG ID is not recognized.
    """
    try:
        return Particle.from_pdgid(pdgid).name
    except ParticleNotFound:
        return f"Unknown({int(pdgid)})"
    except InvalidParticle:
        return f"Invalid({int(pdgid)})"


def is_AZ(arg):
    """
    Check if the input is a tuple of mass and charge number.

    Args:
        arg : The input to check.

    Returns:
        bool: True if the input is a sequence of two integers, False otherwise.
    """
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
    return PDGID(pdg_id)


def process_particle(x):
    """
    Process a particle specification and return a PDGID object.

    Args:
        x: A particle specification. Can be an integer, string,
        PDGID object, or CompositeTarget object.

    Returns:
        A PDGID object representing the particle.

    Raises:
        ValueError: If the particle specification is not recognized.
    """
    if isinstance(x, (PDGID, CompositeTarget)):
        return x
    if isinstance(x, int):
        return PDGID(x)
    if isinstance(x, str):
        try:
            return PDGID(name2pdg(x))
        except KeyError:
            raise ValueError(f"particle with name {x} not recognized")
    if is_AZ(x):
        return PDGID(AZ2pdg(*x))
    raise ValueError(f"{x} is not a valid particle specification")


def fortran_chars(array_ref, char_seq):
    """Helper to set fortran character arrays with python strings"""
    info(10, "Setting fortran array with", char_seq)
    len_arr = int(str(array_ref.dtype)[2:])
    len_seq = len(char_seq)
    return char_seq + (len_arr - len_seq) * " "


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
    """Print to console if min_dbg_level <= chromo.debug_level.

    The function determines automatically the name of caller and appends
    the message to it. Message can be a tuple of strings or objects
    which can be converted to string using `str()`.

    Args:
        min_dbg_level (int): Minimum debug level in config for printing
        message (tuple): Any argument or list of arguments that casts to str
    """
    import chromo

    if min_dbg_level <= chromo.debug_level:
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
    ) as probar:
        task_id = probar.add_task(f"Downloading {fname}", total=total_size)

        with open(outfile, "wb") as f:
            chunk = True
            while chunk:
                chunk = response.read(blocksize)
                f.write(chunk)
                nchunk = len(chunk)
                wrote += nchunk
                probar.advance(task_id, nchunk)

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
            shutil.move(str(model_dir), str(temp_dir))
        _download_file(zip_file, url)
        if zipfile.is_zipfile(zip_file):
            with zipfile.ZipFile(zip_file, "r") as zf:
                zf.extractall(base_dir.as_posix())
            zip_file.unlink()
            shutil.rmtree(temp_dir, ignore_errors=True)

        version_glob = vname.split("_v")[0]
        for vfile in model_dir.glob(f"{version_glob}_v*"):
            vfile.unlink()
        with open(version_file, "w", encoding="utf-8") as vf:
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


def _select_mothers(mask, mothers):
    # This algorithm is slow in pure Python and should be
    # speed up by compiling the logic.

    # attach parentless particles to beam particles,
    # unless those are also removed
    fallback = (-1, -1)
    if mask[0] and mask[1]:
        fallback = (0, 1)

    n = len(mothers)
    indices = np.arange(n)[mask]
    result = mothers[mask]
    mapping = {old: i for i, old in enumerate(indices)}

    n = len(result)
    for i in range(n):
        a = result[i, 0]
        if a == -1:
            continue
        p = mapping.get(a, -1)
        if p == -1:
            a, b = fallback
            result[i, 0] = a
            result[i, 1] = b
        elif p != a:
            q = -1
            b = result[i, 1]
            if b > -1:
                q = mapping.get(b, -1)
            result[i, 0] = p
            result[i, 1] = q
    return result


def select_mothers(arg, mothers):
    if mothers is None:
        return None

    n = len(mothers)

    if isinstance(arg, np.ndarray) and arg.dtype is bool:
        mask = arg
    else:
        mask = np.zeros(n, dtype=bool)
        mask[arg] = True

    with warnings.catch_warnings():
        # suppress numba safety warning that we can ignore
        warnings.simplefilter("ignore")
        return _select_mothers(mask, mothers)


try:
    # accelerate with numba if numba is available
    import numba as nb

    _select_mothers = nb.njit(_select_mothers)

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
    from chromo import models
    from chromo.common import MCRun

    if skip is None:
        skip = []
    if platform.system() == "Windows":
        skip = list(skip) + [
            models.UrQMD34,
            models.EposLHC,
            models.Pythia8,
            models.DpmjetIII191,
            models.DpmjetIII193,
            models.Phojet191,
            models.Phojet193,
        ]

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


def naneq(a, b, rtol=None):
    """
    Return True if a == b or if a and b are both NaN.

    Parameters
    ----------
    a : float
        First float.
    b : float
        Second float.
    """
    if rtol is not None:
        return np.isclose(a, b, rtol=rtol) or (np.isnan(a) and np.isnan(b))
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


class Nuclei:
    """
    Class to specify ranges of nuclei supported by a model.

    It acts like a set and can be combined with sets via operator |.

    The default is to accept any nucleus.
    """

    def __init__(
        self, *, a_min: int = 1, a_max: int = 1000, z_min: int = 0, z_max: int = 1000
    ):
        self._a_min = a_min
        self._a_max = a_max
        self._z_min = z_min
        self._z_max = z_max
        self._other = set()

    def __contains__(self, pdgid: Union[int, CompositeTarget]):
        if pdgid in self._other:
            return True
        if not isinstance(pdgid, PDGID):
            pdgid = PDGID(pdgid)
        if pdgid.A is None:
            return False
        return (
            self._a_min <= pdgid.A <= self._a_max
            and self._z_min <= pdgid.Z <= self._z_max
        )

    def __repr__(self):
        s = (
            f"Nuclei(a_min={self._a_min}, a_max={self._a_max}, "
            f"z_min={self._z_min}, z_max={self._z_max})"
        )
        if self._other:
            s += f" | {self._other}"
        return s

    def __ior__(self, other: Set[PDGID]):
        self._other |= other
        return self

    def __or__(self, other: Set[PDGID]):
        from copy import deepcopy

        result = deepcopy(self)
        result |= other
        return result

    def __ror__(self, other: Set[PDGID]):
        return self.__or__(other)


def unique_sorted_pids(ids):
    """np.arrays of unique ids sorted by abs value
    with negative value first, e.g.
    -11, 11, -12, 12, ...
    """
    uids = np.unique(np.fromiter(ids, dtype=np.int64))
    return uids[np.argsort(2 * np.abs(uids) - (uids < 0))]


def select_long_lived(tau=0, mm=False):
    """
    Returns unstable particles that are stable
    for `tau` sec (or mm if mm=True).
    By default returns all unstable particles excluding nuclei.
    """
    if not mm:
        tau = tau * sec2cm * 1e1  # in mm

    long_lived = []
    for p in Particle.findall():
        pid = int(p.pdgid)
        ctau = p.ctau
        if (
            (ctau is not None)
            and (not math.isinf(ctau))
            and (not math.isnan(ctau))
            and (abs(pid) < 1000000000)  # exclude nuclei
            and (ctau > tau)
        ):
            long_lived.append(pid)

    return long_lived
