"""Utility module for auxiliary methods and classes."""

import copy
import dataclasses
import importlib.util
import inspect
import math
import pkgutil
import shutil
import urllib.request
import warnings
import zipfile
from collections.abc import Collection, Sequence
from enum import Enum
from pathlib import Path
from typing import Union

import numpy as np
from particle import PDGID, InvalidParticle, Particle, ParticleNotFound

import chromo
from chromo.constants import MeV, nucleon_mass, sec2cm

EventFrame = Enum("EventFrame", ["CENTER_OF_MASS", "FIXED_TARGET", "GENERIC"])


@dataclasses.dataclass(init=False)
class CompositeTarget:
    """Definition of composite targets made of multiple (atomic) nuclei.

    Examples of such composite targets are Air, CO_2, HCl, C_2H_60.
    """

    label: str
    components: tuple[PDGID]
    fractions: np.ndarray

    def __init__(
        self, components: Collection[tuple[Union[str, int], float]], label: str = ""
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
                msg = f"component {particle} is not a nucleus"
                raise ValueError(msg)
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
            msg = f"cannot get mass for {pdgid}"
            raise ValueError(msg)
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
    if pdgid == 2112:
        return 1, 0
    if pdgid == 2212:
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
            msg = f"particle with name {x} not recognized"
            raise ValueError(msg)
    if is_AZ(x):
        return PDGID(AZ2pdg(*x))
    msg = f"{x} is not a valid particle specification"
    raise ValueError(msg)


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
        print(caller_name(), *args)  # noqa: T201


def _download_file(outfile, url):
    """Download a file from 'url' to 'outfile'"""
    from rich.progress import (
        BarColumn,
        DownloadColumn,
        Progress,
        SpinnerColumn,
        TextColumn,
        TimeRemainingColumn,
    )

    fname = Path(url).name
    try:
        response = urllib.request.urlopen(url)
    except BaseException:
        msg = f"_download_file: probably something wrong with url = '{url}'"
        raise ConnectionError(msg)
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
        msg = f"{fname} has not been downloaded"
        raise ConnectionError(msg)


# Helper function to extract zip files
def _extract_zip(zip_file_path: Path, destination_dir: Path):
    """Extracts a zip file to a destination directory."""
    if not zip_file_path.is_file():
        msg = f"Zip file not found: {zip_file_path}"
        raise OSError(msg)
    if not zipfile.is_zipfile(zip_file_path):
        msg = f"File {zip_file_path} is not a valid zip file."
        raise OSError(msg)
    with zipfile.ZipFile(zip_file_path, "r") as zf:
        zf.extractall(destination_dir)
    info(1, f"Successfully extracted {zip_file_path.name} to {destination_dir}")


# REFACTORED _cached_data_dir function
def _cached_data_dir(url: str) -> str:
    """
    Ensures model data is available in src/chromo/iamdata/model_name.

    1. Checks for a version file (e.g., model_name_vXXX) in the target directory.
       If present, assumes data is correct and returns the path.
    2. If version file not present, prepares the target model directory by clearing
       any existing content and recreating it.
    3. Tries to obtain the source ZIP:
        - If "CI" env var is set, copies from ~/.cache/chromo/zip_filename.
        - Else, downloads from the given URL.
       The ZIP is temporarily placed in src/chromo/iamdata/.
    4. Extracts the ZIP to src/chromo/iamdata/model_name/.
    5. Cleans up the temporary ZIP and any old version files.
    6. Creates the new version file.

    Args:
        url: URL for the model data zip file.

    Returns:
        Path to the model data directory (e.g., src/chromo/iamdata/model_name/).

    Raises:
        ConnectionError: If the zip file cannot be obtained (download/copy failed).
        OSError: For file operation issues (e.g., extraction, cleanup).
    """
    iamdata_dir = Path(__file__).parent.absolute() / "iamdata"
    iamdata_dir.mkdir(parents=True, exist_ok=True)

    zip_filename = Path(url).name
    vname_stem = Path(url).stem  # e.g., qgsjet_v002 (used for version file)
    model_name = vname_stem.split("_v")[0]  # e.g., qgsjet (used for dir name)

    model_dir = iamdata_dir / model_name  # Target directory for extracted data
    version_file = (
        model_dir / vname_stem
    )  # Expected version file, e.g., iamdata/qgsjet/qgsjet_v002

    if version_file.exists():
        info(
            2,
            f"Version file {version_file.name} found for {model_name}. Using existing data.",
        )
        return str(model_dir) + "/"

    # Prepare model_dir: remove if exists, then recreate for a clean slate
    if model_dir.exists():
        info(1, f"Removing existing model directory: {model_dir}")
        try:
            shutil.rmtree(model_dir)
        except OSError as e:
            warnings.warn(
                f"Could not remove existing model directory {model_dir}: {e}. "
                "Proceeding, but extraction might fail or be incomplete."
            )
    try:
        model_dir.mkdir(parents=True, exist_ok=True)
    except OSError as e:
        msg = f"Could not create model directory {model_dir}: {e}"
        raise OSError(msg) from e

    # temp_zip_path is where the zip will be before extraction (inside iamdata_dir)
    temp_zip_path = iamdata_dir / zip_filename
    zip_obtained = False

    try:
        # Attempt to get ZIP from CI cache or download
        ci_cache_file = Path.home() / ".cache" / "chromo" / zip_filename
        info(1, f"CI: Checking cache: {ci_cache_file}")
        if ci_cache_file.is_file():
            try:
                shutil.copy2(ci_cache_file, temp_zip_path)
                info(1, f"CI: Copied {zip_filename} from cache to {temp_zip_path}")
                zip_obtained = True
            except Exception as e:
                warnings.warn(
                    f"CI: Failed to copy {zip_filename} from cache: {e}. Will try download."
                )
        else:
            info(1, f"CI: Cache miss for {zip_filename}. Will try download.")

        if not zip_obtained:
            info(1, f"Downloading {url} to {temp_zip_path}")
            _download_file(temp_zip_path, url)  # _download_file is an existing function
            zip_obtained = True

        if not zip_obtained or not temp_zip_path.is_file():
            # This case should ideally be caught by _download_file raising an error
            msg = f"Failed to obtain zip file: {zip_filename}"
            raise ConnectionError(msg)  # noqa: TRY301

        # Extract the zip
        info(1, f"Extracting {temp_zip_path.name} to {model_dir}")
        _extract_zip(temp_zip_path, model_dir.parent)

        # Post-extraction: Clean up old version files and create new one
        for old_vfile in model_dir.glob(f"{model_name}_v*"):
            if old_vfile.name != vname_stem:
                try:
                    old_vfile.unlink()
                    info(5, f"Removed old version artifact: {old_vfile.name}")
                except OSError as e:
                    warnings.warn(
                        f"Could not remove old version file {old_vfile.name}: {e}"
                    )

        try:
            with open(version_file, "w", encoding="utf-8") as vf:
                vf.write(url)  # Store the URL as the content of the version file
            info(1, f"Created version file: {version_file.name} for {model_name}")
        except OSError as e:
            # If version file creation fails, the data is there but might be re-processed.
            # This is a critical step for recognizing the data next time.
            msg = f"Failed to create version file {version_file.name} for {model_name}: {e}"
            raise OSError(msg) from e

        return str(model_dir) + "/"

    except (ConnectionError, OSError) as e:
        # If any crucial step failed, log it and re-raise.
        # The model_dir might be in an incomplete state.
        # On next run, version_file check will fail, and it will try to rebuild.
        warnings.warn(f"Error during data caching for {model_name} from {url}: {e}")
        # Clean up model_dir to ensure a fresh start next time if it exists
        if model_dir.exists():
            shutil.rmtree(model_dir, ignore_errors=True)
        raise  # Re-raise the caught specific exception (ConnectionError or OSError)
    except Exception as e:
        # Catch any other unexpected errors
        warnings.warn(f"Unexpected error during data caching for {model_name}: {e}")
        if model_dir.exists():
            shutil.rmtree(model_dir, ignore_errors=True)
        # Wrap in a RuntimeError for unexpected issues
        msg = f"Unexpected error processing {model_name}: {e}"
        raise RuntimeError(msg) from e
    finally:
        # Always try to clean up the temporary zip file from iamdata_dir
        if temp_zip_path.exists():
            try:
                temp_zip_path.unlink()
                info(5, f"Cleaned up temporary zip: {temp_zip_path.name}")
            except OSError as e:
                warnings.warn(
                    f"Could not remove temporary zip {temp_zip_path.name}: {e}"
                )


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
        if isinstance(other, cls):
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
    """
    Select mothers using the best available implementation.
    Uses numba-accelerated loop if available, otherwise numpythonic version.
    """
    try:
        import numba as nb

        return nb.njit(_select_mothers_loop)(mask, mothers)
    except ModuleNotFoundError:
        return _select_mothers_numpy(mask, mothers)


def _select_mothers_loop(mask, mothers):
    # Original for-loop version (numba-compatible)
    fallback = (-1, -1)
    if mask[0] and mask[1]:
        fallback = (0, 1)
    n = len(mothers)
    indices = np.arange(n)[mask]
    result = mothers[mask].copy()
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


def _select_mothers_numpy(mask, mothers):
    # Numpythonic version (not numba-compatible)
    fallback = (-1, -1)
    if mask[0] and mask[1]:
        fallback = (0, 1)
    n = len(mothers)
    indices = np.arange(n)[mask]
    result = mothers[mask].copy()
    mapping = {old: i for i, old in enumerate(indices)}
    a = result[:, 0]
    b = result[:, 1]
    p = np.array([mapping.get(x, -1) for x in a])
    q = np.array([mapping.get(x, -1) if x > -1 else -1 for x in b])
    mask_a_valid = a != -1
    fallback_rows = mask_a_valid & (p == -1)
    result[fallback_rows, 0] = fallback[0]
    result[fallback_rows, 1] = fallback[1]
    update_rows = mask_a_valid & (p != -1) & (p != a)
    result[update_rows, 0] = p[update_rows]
    result[update_rows, 1] = q[update_rows]
    return result


def select_mothers(arg, mothers):
    """
    Select mothers using the best available implementation.
    Uses numba-accelerated loop if available, otherwise numpythonic version.
    Handles boolean mask or index array for selection, and ignores numba warnings.
    """
    if mothers is None:
        return None
    n = len(mothers)
    if isinstance(arg, np.ndarray) and arg.dtype is bool:
        mask = arg
    else:
        mask = np.zeros(n, dtype=bool)
        mask[arg] = True
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # ignore numba warnings if any
        return _select_mothers(mask, mothers)


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


def get_available_binary_modules():
    """
    Get all available binary modules in chromo.models.

    This function lists all importable binary modules in the chromo.models
    package and returns a list of their names.
    """

    # List all importable binary modules in chromo.models,
    # excluding those starting with '_'
    available_binaries = []
    for _, name, _ in pkgutil.iter_modules(chromo.models.__path__):
        if not name.startswith("_"):
            continue
        spec = importlib.util.find_spec(f"chromo.models.{name}")
        if spec and spec.origin and not spec.origin.endswith(".py"):
            available_binaries.append(name)

    return available_binaries


def get_all_models(only_names=False):
    """
    Get all available models in chromo.models.

    This function lists all importable binary modules in the chromo.models
    package. It then checks which model classes can be used and returns
    a list of classes of class names in `only_names` attribute is set.
    """
    import inspect

    import chromo.models
    from chromo.common import MCRun

    # List all importable binary modules in chromo.models,
    # excluding those starting with '_'
    available_binaries = get_available_binary_modules()

    active_classes = []
    active_class_names = []
    for key in dir(chromo.models):
        obj = getattr(chromo.models, key)
        if not inspect.isclass(obj):
            continue
        if issubclass(obj, MCRun) and hasattr(obj, "_library_name"):
            library_name = obj._library_name
            if library_name in available_binaries:
                active_classes.append(obj)
                active_class_names.append(obj.__name__)

    return active_class_names if only_names else active_classes


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

    def __ior__(self, other: set[PDGID]):
        self._other |= other
        return self

    def __or__(self, other: set[PDGID]):
        from copy import deepcopy

        result = deepcopy(self)
        result |= other
        return result

    def __ror__(self, other: set[PDGID]):
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
