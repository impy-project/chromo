import numpy as np
from impy.common import MCRun, MCEvent, CrossSectionData
from particle import literals as lp
from impy.kinematics import EventFrame
from impy.constants import tau_stable, sec2cm


class PYTHIA6Event(MCEvent):
    """Wrapper class around HEPEVT particle stack."""

    def _charge_init(self, npart):
        k = self._lib.pyjets.k[:npart, 1]
        # TODO accelerate by implementing this loop in Fortran
        return np.fromiter((self._lib.pychge(ki) / 3 for ki in k), np.double)


class Pythia6(MCRun):
    """Implements all abstract attributes of MCRun for the
    EPOS-LHC series of event generators."""

    _name = "Pythia"
    _version = "6.428"
    _library_name = "_pythia6"
    _event_class = PYTHIA6Event
    _frame = EventFrame.CENTER_OF_MASS
    _projectiles = {
        p.pdgid: code
        for (p, code) in (
            (lp.p, "p"),
            (lp.n, "n"),
            (lp.K_plus, "K+"),
            (lp.pi_plus, "pi+"),
        )
    }
    _targets = _projectiles

    def __init__(self, evt_kin, seed=None, new_mpi=False):
        super().__init__(seed)

        # setup logging
        lun = 6  # stdout
        self._lib.pydat1.mstu[10] = lun

        if new_mpi:
            assert False, "new_mpi=True is currently broken"
            # Pythia output:
            #   Error: you did not link PDFLIB correctly.
            #   Dummy routine PDFSET in PYTHIA file called instead.

            # Latest Pythia 6 is tune 383
            self._lib.pytune(383)
            self._event_call = self._lib.pyevnw
        else:
            self._event_call = self._lib.pyevnt

        # self.mstp[51]

        # Setup pythia processes (set to custom mode)
        self._lib.pysubs.msel = 0

        # Enable minimum bias processes incl diffraction, low-pt
        # but no elastic (see p227 of hep-ph/0603175)
        for isub in [11, 12, 13, 28, 53, 68, 92, 93, 94, 95, 96]:
            self._lib.pysubs.msub[isub - 1] = 1

        self.kinematics = evt_kin

        # Set default stable
        self._set_final_state_particles()
        # Set PYTHIA decay flags to follow all changes to MDCY
        self._lib.pydat1.mstj[21 - 1] = 1
        self._lib.pydat1.mstj[22 - 1] = 2
        self._lib.pydat1.parj[70] = tau_stable * sec2cm * 10.0  # mm

    def _cross_section(self):
        s = self._lib.pyint7.sigt[0, 0]
        c = CrossSectionData(
            total=s[0],
            elastic=s[1],
            inelastic=s[0] - s[1],
            diffractive_xb=s[2],
            diffractive_ax=s[3],
            diffractive_xx=s[4],
            diffractive_axb=0,
        )
        return c

    def _set_kinematics(self, kin):
        codes = []
        for pdg in (kin.p1, kin.p2):
            c = self._projectiles[abs(pdg)]
            if abs(pdg) != pdg:
                last = c[-1]
                if last == "+":
                    c = c[:-1] + "-"
                elif last == "-":
                    c = c[:-1] + "+"
                else:
                    c += "bar"
            codes.append(c)
        self._lib.pyinit("CMS", *codes, kin.ecm)

    def _set_stable(self, pdgid, stable):
        kc = self._lib.pycomp(pdgid)
        self._lib.pydat3.mdcy[kc - 1, 0] = 0 if stable else 1

    def _generate(self):
        self._event_call()
        self._lib.pyhepc(1)
        return True
