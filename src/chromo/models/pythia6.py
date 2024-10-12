import numpy as np
from chromo.common import MCRun, MCEvent, CrossSectionData
from particle import literals as lp
from chromo.kinematics import EventFrame
from chromo.constants import standard_projectiles


class PYTHIA6Event(MCEvent):
    """Wrapper class around HEPEVT particle stack."""

    def _charge_init(self, npart):
        k = self._lib.pyjets.k[:npart, 1]
        # TODO accelerate by implementing this loop in Fortran
        return np.fromiter((self._lib.pychge(ki) / 3 for ki in k), np.double)

    def _repair_initial_beam(self):
        self.status[0:2] = 4


class Pythia6(MCRun):
    """Implements all abstract attributes of MCRun for the
    EPOS-LHC series of event generators."""

    _name = "Pythia"
    _version = "6.428"
    _library_name = "_pythia6"
    _event_class = PYTHIA6Event
    _frame = None
    _projectiles = standard_projectiles
    _targets = standard_projectiles

    def __init__(self, evt_kin, *, seed=None, new_mpi=False):
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

        self._set_final_state_particles()

    def _cross_section(self, kin=None, max_info=False):
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
            code = {
                lp.proton.pdgid: "p",
                lp.neutron.pdgid: "n",
                lp.K_plus.pdgid: "K+",
                lp.pi_plus.pdgid: "pi+",
                lp.K_S_0.pdgid: "K_S0",
                lp.K_L_0.pdgid: "K_L0",
            }[abs(pdg)]
            if abs(pdg) != pdg:
                last = code[-1]
                if last == "+":
                    code = code[:-1] + "-"
                elif last == "-":
                    code = code[:-1] + "+"
                else:
                    code += "bar"
            codes.append(code)
        if kin.frame == EventFrame.FIXED_TARGET:
            self._frame = EventFrame.FIXED_TARGET
            self._lib.pyinit("FIXT", *codes, kin.plab)
        else:
            self._frame = EventFrame.CENTER_OF_MASS
            self._lib.pyinit("CMS", *codes, kin.ecm)

    def _set_stable(self, pdgid, stable):
        kc = self._lib.pycomp(pdgid)
        self._lib.pydat3.mdcy[kc - 1, 0] = 0 if stable else 1

    def _generate(self):
        self._event_call()
        self._lib.pyhepc(1)
        return True
