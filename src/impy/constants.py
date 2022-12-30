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
millibarn = 1.0
microbarn = 1e-3 * millibarn

nucleon_mass = 0.5 * (lp.proton.mass + lp.neutron.mass) * MeV

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

# standard long-lived particles with life-time > 30 ps
long_lived = (
    13,
    2112,
    130,
    211,
    310,
    321,
    -3334,
    -3322,
    -3312,
    -3222,
    -3122,
    -3112,
    3112,
    3122,
    3222,
    3312,
    3322,
    3334,
    -321,
    -211,
    -2112,
    -13,
)

# only positive PDGIDs!
standard_projectiles = {
    p.pdgid for p in (lp.p, lp.n, lp.pi_plus, lp.K_plus, lp.K_S_0, lp.K_L_0)
}

# only positive PDGIDs!
em_particles = {p.pdgid for p in (lp.photon, lp.e_minus)}

# # Standard stable particles for for fast air shower cascade calculation
# # Particles with an anti-partner
# standard_particles = [11, 13, 15, 211, 321, 2212, 2112, 3122, 411, 421, 431]
# standard_particles += [-pid for pid in standard_particles]
# # unflavored particles
# standard_particles = tuple(standard_particles + [111, 130, 310, 221, 223, 333])

# Default definition of final state particle
# All particles with proper lifetime shorter than this
# will decay
tau_stable = 30e-12  # 30 ps, typical value at LHC


# Air composition for special cross section functions
# (source https://en.wikipedia.org/wiki/Atmosphere_of_Earth)
air_composition = {
    1000070140: 0.78084,  # nitrogen
    1000080160: 0.20946,  # oxygen
    1000180400: 0.00934,  # argon
}
