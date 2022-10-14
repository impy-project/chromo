# This dependency might be overkill for just reading a few
# variables. Should be changed at some point.
from scipy.constants import speed_of_light as _c_SI
from particle import literals as lp

c = 1e2 * _c_SI
cm2sec = 1e-2 / _c_SI
sec2cm = _c_SI * 1e2
eV = 1e-9
keV = 1e-6
MeV = 1e-3
GeV = 1.0
TeV = 1e3
PeV = 1e6
EeV = 1e9

quarks_and_diquarks_and_gluons = (
    1,
    2,
    3,
    4,
    5,
    6,
    21,
    1103,
    2101,
    2103,
    2203,
    3101,
    3103,
    3201,
    3203,
    3303,
    4101,
    4103,
    4201,
    4203,
    4301,
    4303,
    4403,
    5101,
    5103,
    5201,
    5203,
    5301,
    5303,
    5401,
    5403,
    5503,
)

nucleon_mass = 0.5 * (lp.proton.mass + lp.neutron.mass) / 1e3
