import os
import sys
sys.path.append(os.dirname(__file__) + "/..")
sys.path.append("../ParticleDataTool")
import numpy as np
from common import EventKinematics
from ParticleDataTool import PYTHIAParticleData
from collections import Sequence, namedtuple
Nucleus = namedtuple("Nucleus", "A Z")


persistent_obj = None

def epos(version, projectile, target, nevents, **kwargs):
    __process_version(version, {"lhc": "eposlhc"})
    from epos import EPOSMCRun

    run = None

def sibyll(version, projectile, target, nevents, **kwargs):
    __process_version(version, {"2.3c", "sib23c"})
    from sibyll import SibyllMCRun
    id1, p1, id2, p2 = __process_particles(projectile, target)
    evtkin = EventKinematics(p1pdg=id1, p2pdg=id2, pab )
    run = SibyllMCRun(evtkin, )
    if persistent_obj is None:
        # initial
    while run.is_going_on():
        # somehow get event
        yield event

__pdg_db__ = PYTHIAParticleData()


def __interpret_particle_id(arg):
    if type(arg) is str:
        return __pdg_db__.pdg_id(arg)
    if type(arg) is Nucleus:
        return arg
    if type(arg) is int:
        return arg
    raise ValueError("argument not understood")


def __process_version(version, version_map):
    key = version.lower()
    if key not in version_map:
        raise ValueError(version + " not recognized")
    return __import__(version_map[key])


def __process_particles(projectile, target):
    id1, p1 = projectile
    id2, p2 = target
    id1 = __interpret_particle_id(id1)
    id2 = __interpret_particle_id(id2)
    return id1, p1, id2, p2
