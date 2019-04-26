"""Utility module for auxiliary methods and classes."""

from __future__ import print_function
import inspect

from impy.common import impy_config

# Global debug flags that would be nice to have in some sort
# of config or other ideas?
print_module = True

# Standard stable particles for for fast air shower cascade calculation
# Particles with an anti-partner
standard_particles = [
    11, 13, 15, 211, 321, 2212, 2112, 3122, 411, 421, 431]
standard_particles += [-pid for pid in standard_particles]
#unflavored particles
standard_particles = tuple(
    standard_particles + [111, 130, 310, 221, 223, 333])


def getAZN(pdgid):
    """Returns mass number :math:`A`, charge :math:`Z` and neutron
    number :math:`N` of ``pdgid``.

    Note::
    
        PDG ID for nuclei is coded according to 10LZZZAAAI. For iron-52 it is 1000260520.

    Args:
        pdgid (int): PDG ID of nucleus/mass group
    Returns:
        (int,int,int): (Z,A) tuple
    """
    Z, A = 1, 1
    if pdgid < 2000:
        return 0,0,0
    elif pdgid == 2112:
        return 1,0,1
    elif pdgid == 2212:
        return 1,1,0
    elif pdgid > 1000000000:
        A = pdgid % 1000 / 10
        Z = pdgid % 1000000 / 10000
        return A, Z, A - Z
    else:
        return 1, 0, 0

def AZ2pdg(A, Z):
    """Conversion of nucleus with mass A and chage Z
    to PDG nuclear code"""
    # 10LZZZAAAI
    pdg_id = 1000000000
    pdg_id += 10*A
    pdg_id += 10000*Z
    return pdg_id


def clear_and_set_fortran_chars(array_ref, char_seq):
    """Helper to set fortran character arrays with python strings"""
    info(10, 'Setting fortran array with', char_seq)
    # Reset
    array_ref.data[:] = len(array_ref.data) * ' '
    array_ref.data[:len(char_seq)] = char_seq

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

    if print_module:
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

    # Would prefer here a global debug
    # if min_dbg_level <= config["debug_level"]:
    if min_dbg_level <= impy_config['debug_level']:
        message = [str(m) for m in message]
        print(caller_name() + " ".join(message))
