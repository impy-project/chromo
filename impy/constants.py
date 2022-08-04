# This dependency might be overkill for just reading a few
# variables. Should be changed at some point.
import scipy.constants as spc

c = 1e2 * spc.c
cm2sec = 1e-2 / spc.c
sec2cm = spc.c * 1e2
eV = 1e-9
keV = 1e-6
MeV = 1e-3
GeV = 1.0
TeV = 1e3
PeV = 1e6
EeV = 1e9
