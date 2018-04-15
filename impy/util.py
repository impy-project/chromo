"""Utility module for auxiliary methods and classes."""

import inspect
# This dependency might be overkill for just reading a few
# variables. Should be changed at some point.
import scipy.constants as spc


#===============================================================================
# Standard stable particles for for fast air shower cascade calculation
#===============================================================================
# Particles having an anti-partner
standard_particles = [11, 13, 15, 211, 321, 2212, 2112, 3122, 411, 421, 431]
standard_particles += [-pid for pid in standard_particles]

#unflavored particles
standard_particles = tuple(standard_particles + [111, 130, 310, 221, 223, 333])

def convert_to_namedtuple(dictionary, name='GenericNamedTuple'):
    """Converts a dictionary to a named tuple."""
    from collections import namedtuple
    return namedtuple(name, dictionary.keys())(**dictionary)


# The stuff below is from an astrophysical code that I write,
# and most are not needed for particle physics. But since
# Hans added some unit module here we go.

# Default units in impy are ***cm, s, GeV***
# Define here all constants and unit conversions and use
# throughout the code. Don't write c=2.99.. whatever.
# Write clearly the units returned by a function.
# Convert them if not standard unit.
# Accept only arguments in the units above.

UNITS_AND_CONVERSIONS_DEF = dict(
    c=1e2 * spc.c,
    cm2Mpc=1. / (spc.parsec * spc.mega * 1e2),
    Mpc2cm=spc.mega * spc.parsec * 1e2,
    m_proton=spc.physical_constants['proton mass energy equivalent in MeV'][0]
    * 1e-3,
    m_electron=spc.physical_constants[
        'electron mass energy equivalent in MeV'][0] * 1e-3,
    r_electron=spc.physical_constants['classical electron radius'][0] * 1e2,
    fine_structure=spc.fine_structure,
    GeV2erg=1. / 624.15,
    erg2GeV=624.15,
    km2cm=1e5,
    yr2sec=spc.year,
    Gyr2sec=spc.giga * spc.year,
    cm2sec=1e-2 / spc.c,
    sec2cm=spc.c * 1e2)

# This is the immutable unit object to be imported throughout the code
IMPY_UNITS = convert_to_namedtuple(UNITS_AND_CONVERSIONS_DEF, "ImpyUnits")


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
        return ''

    parentframe = stack[start][0]

    name = []

    if config["print_module"]:
        module = inspect.getmodule(parentframe)
        # `modname` can be None when frame is executed directly in console
        if module:
            name.append(module.__name__ + '.')

    # detect classname
    if 'self' in parentframe.f_locals:
        # I don't know any way to detect call from the object method
        # there seems to be no way to detect static method call - it will
        # be just a function call

        name.append(parentframe.f_locals['self'].__class__.__name__ + '::')

    codename = parentframe.f_code.co_name
    if codename != '<module>':  # top level usually
        name.append(codename + '(): ')  # function or a method

    del parentframe
    return "".join(name)


def info(min_dbg_level, *message):
    """Print to console if `min_debug_level <= config["debug_level"]`

    The fuction determines automatically the name of caller and appends
    the message to it. Message can be a tuple of strings or objects
    which can be converted to string using `str()`.

    Args:
        min_dbg_level (int): Minimum debug level in config for printing
        message (tuple): Any argument or list of arguments that casts to str
    """

    if min_dbg_level <= config["debug_level"]:
        message = [str(m) for m in message]
        print(caller_name() + " ".join(message))
