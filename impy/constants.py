# This dependency might be overkill for just reading a few
# variables. Should be changed at some point.
from scipy.constants import speed_of_light as _c_SI

c = 1e2 * _c_SI
cm2sec = 1e-2 / _c_SI
sec2cm = _c_SI * 1e2
eV = 1e-9
keV = 1e-6
MeV = 1e-3
GeV = 1.
TeV = 1e3
PeV = 1e6
EeV = 1e9
