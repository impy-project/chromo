"""DPMJET-rs: DPMJET-III / PHOJET physics written anew in Rust.

The generator lives in the DPMJET-rs project (crate ``dpm``), which
installs the Python module ``_dpmjet_rs`` (``tools/install_python.sh``
there); it is not built with chromo.
"""

import numpy as np
from particle import literals as lp

from chromo.common import CrossSectionData, MCEvent, MCRun
from chromo.constants import GeV
from chromo.kinematics import EventFrame
from chromo.util import Nuclei, is_real_nucleus

# Hadron beams with DPMJET parameters (on protons and neutrons); K0_L and
# K0_S collide as K0 or anti-K0 with equal probability. Anti-particles are
# checked per pair in _check_kinematics.
_hadrons = {
    p.pdgid
    for p in (
        lp.p,
        lp.n,
        lp.pi_plus,
        lp.pi_0,
        lp.K_plus,
        lp.K_0,
        lp.K_S_0,
        lp.K_L_0,
        lp.Lambda,
        lp.Sigma_minus,
    )
}

# nucleon codes of chromo's nuclear PDG IDs
_nucleons = {1000010010: 2212, 1000000010: 2112}


def _beam(p):
    """(pdg, A, Z) of a chromo beam particle."""
    if is_real_nucleus(p):
        return (2212, int(p.A), int(p.Z))
    pdg = int(p)
    return (_nucleons.get(pdg, pdg), 1, 0)


class DpmjetRSEvent(MCEvent):
    """Event of ``DpmjetRS``: the final hadrons and the residual nuclei."""

    def __init__(self, generator):
        super().__init__(generator)
        self._prepend_initial_beam()
        self.mothers[2:] = [0, 1]
        self.daughters[:] = -1

    def _get_charge(self, npart):
        return self._lib.ichg[:npart]

    def _get_impact_parameter(self):
        return self._lib.b

    def _get_n_wounded(self):
        return self._lib.nwounded


class DpmjetRS(MCRun):
    """DPMJET-III / PHOJET (DPMJET 19.4 physics) written anew in Rust.

    Hadron-hadron, hadron-nucleus and nucleus-nucleus collisions in the Dual
    Parton Model: two-channel eikonal with soft and hard pomerons,
    diffraction, minijets, Glauber geometry, chain ends of wounded nucleons,
    Lund string fragmentation with PYTHIA 6.4 physics, PYTHIA's decays of
    unstable hadrons before the formation-zone intranuclear cascade (as
    DPMJET, following ``set_stable``), the cascade and the nuclear
    potential. Below 4-6 GeV/c hadron projectiles make their first
    collision with HADRES, a low-energy hadron-nucleon model (resonances,
    two-body channels, strings), and cascade with it. Validated statistically
    against DPMJET 19.4 above that window. Photon beams
    and the evaporation of residual nuclei are not implemented: residual
    nuclei carry the baryon number and charge the hadrons do not. Decays
    of charm and bottom hadrons are left to chromo's decay handler. Uses
    the module ``_dpmjet_rs`` of the DPMJET-rs project, with its own random
    generator seeded from chromo's; tables are built at the first event of
    a beam pair at a new highest energy (about a second) and kept.
    """

    _name = "DPMJET-rs"
    _version = "19.4-rs.2"
    _library_name = "_dpmjet_rs"
    _event_class = DpmjetRSEvent
    _frame = EventFrame.CENTER_OF_MASS
    _projectiles = _hadrons | Nuclei(a_max=280)
    _targets = Nuclei()
    # hadron projectiles below p_lab = 4-6 GeV/c take the HADRES path
    # (one collision of the low-energy hadron-nucleon model HADRES, then the
    # cascade with it; DPMJET's HADRIN regime written anew); nucleus
    # projectiles keep the DPM path there
    _ecm_min = 2 * GeV
    _restartable = True
    _restore_beam_and_history = False
    randomize_azimuth = True

    def __init__(self, evt_kin, *, seed=None):
        super().__init__(seed)
        self._lib = self._lib.Generator(int(self._rng.integers(0, 2**63 - 1)))
        self.kinematics = evt_kin
        self._set_final_state_particles()
        self._activate_decay_handler(
            on=True, seed=int(self._rng.integers(0, 900_000_000))
        )

    def _check_kinematics(self, kin):
        super()._check_kinematics(kin)
        targets = getattr(kin.p2, "components", [kin.p2])
        for t in targets:
            if not self._lib.supports(_beam(kin.p1), _beam(t)):
                msg = f"{self.pyname}: no parameters for {int(kin.p1)} on {int(t)}"
                raise ValueError(msg)

    def _set_kinematics(self, kin):
        self._lib.set_kinematics(_beam(kin.p1), _beam(kin.p2), kin.ecm)

    def _cross_section(self, kin=None, max_info=False):
        kin = self.kinematics if kin is None else kin
        xs = self._lib.cross_section(_beam(kin.p1), _beam(kin.p2), kin.ecm)
        if "prod" in xs:
            return CrossSectionData(prod=xs["prod"])
        sd = xs["sd_proj"] + xs["sd_targ"]
        return CrossSectionData(
            total=xs["total"],
            inelastic=xs["inelastic"],
            elastic=xs["elastic"],
            prod=xs["inelastic"],
            diffractive_xb=xs["sd_proj"],
            diffractive_ax=xs["sd_targ"],
            diffractive_xx=xs["dd"],
            diffractive_axb=xs["cd"],
            diffractive_sum=sd + xs["dd"] + xs["cd"],
            b_elastic=xs["b_elastic"],
        )

    def _set_stable(self, pdgid, stable):
        # decays before the intranuclear cascade, as DPMJET's MDCY; species
        # with partonic decays (charm, bottom) are left to the decay handler,
        # which must then follow direct calls of _set_stable too
        self._lib.set_stable(pdgid, stable)
        if self._decay_handler is None:
            return
        fsp = np.asarray(self._final_state_particles, dtype=np.int64)
        if stable:
            fsp = np.unique(np.append(fsp, pdgid))
        else:
            fsp = fsp[fsp != pdgid]
        self._final_state_particles = fsp
        self._sync_decay_handler()

    def _generate(self):
        self._lib.event()
        return True

    def set_cascade(self, on=True):
        """Switch the intranuclear cascade on (default) or off."""
        self._lib.set_cascade(on)

    def set_hadres(self, on=True):
        """HADRES path for hadron projectiles below ``hadres.primary_window``
        (default on): off, every event takes the DPM path."""
        self._lib.set_hadres(on)

    def set_low_energy(self, on=True):
        """Former name of :meth:`set_hadres`."""
        self.set_hadres(on)

    def set_decays(self, on=True):
        """Decays before the cascade on (default) or off: with off, the
        cascade sees the primary hadrons and chromo's decay handler decays
        them afterwards."""
        self._lib.set_decays(on)

    def stats(self):
        """Counters of retries and rejections in the elementary collisions."""
        return self._lib.stats()

    def reset_stats(self):
        self._lib.reset_stats()
