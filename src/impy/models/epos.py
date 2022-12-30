import numpy as np
from impy.kinematics import EventFrame
from impy.common import MCEvent, MCRun, CrossSectionData
from impy.util import (
    _cached_data_dir,
    fortran_array_insert,
    fortran_array_remove,
    Nuclei,
)
from impy.constants import standard_projectiles


class EPOSEvent(MCEvent):
    """Wrapper class around EPOS particle stack."""

    def __init__(self, generator):
        super().__init__(generator)
        # EPOS sets parents of beam particles to (-1, -1).
        # We change it to (0, 0)
        self.parents[self.status == 4] = 0

    def _charge_init(self, npart):
        return self._lib.charge_vect(self._lib.hepevt.idhep[:npart])

    def _get_impact_parameter(self):
        # return self._lib.nuc3.bimp
        return float(self._lib.cevt.bimevt)

    def _get_n_wounded(self):
        return int(self._lib.cevt.npjevt), int(self._lib.cevt.ntgevt)


class EposLHC(MCRun):
    """Implements all abstract attributes of MCRun for the
    EPOS-LHC series of event generators."""

    _name = "EPOS"
    _version = "LHC"
    _library_name = "_eposlhc"
    _event_class = EPOSEvent
    _frame = None
    _projectiles = standard_projectiles | Nuclei()
    _data_url = (
        "https://github.com/impy-project/impy"
        + "/releases/download/zipped_data_v1.0/epos_v001.zip"
    )

    def __init__(self, evt_kin, *, seed=None):
        import impy

        super().__init__(seed)

        self._lib.aaset(0)
        datdir = _cached_data_dir(self._data_url)

        lun = 6  # stdout
        self._lib.initepos(
            float(self._seed),
            evt_kin.ecm,
            datdir,
            len(datdir),
            impy.debug_level,
            lun,
        )

        self._set_final_state_particles()
        self._lib.charge_vect = np.vectorize(self._lib.getcharge, otypes=[np.float32])
        self.kinematics = evt_kin

    def _cross_section(self, kin=None):
        total, inel, el, dd, sd, _ = self._lib.xsection()
        return CrossSectionData(
            total=total,
            inelastic=inel,
            elastic=el,
            diffractive_xb=sd / 2,  # this is an approximation
            diffractive_ax=sd / 2,  # this is an approximation
            diffractive_xx=dd,
            diffractive_axb=0,
        )

    def _set_kinematics(self, kin):
        if self._frame == EventFrame.FIXED_TARGET:
            iframe = 2
            self._frame = kin.frame
        else:
            iframe = 1
            self._frame = EventFrame.CENTER_OF_MASS

        self._lib.initeposevt(kin.ecm, -1.0, kin.p1, kin.p2, iframe)

    def _set_stable(self, pdgid, stable):
        # EPOS decays all unstable particles by default. It uses a nodcy common block
        # to prevent decay of particles. The common block contains the array
        # nody and the length nrnody. The array holds EPOS particle ids of
        # particles that should not be decayed.

        idx = self._lib.idtrafo("pdg", "nxs", pdgid)

        c = self._lib.nodcy  # common block
        if stable:
            fortran_array_insert(c.nody, c.nrnody, idx)
        else:
            fortran_array_remove(c.nody, c.nrnody, idx)

    def _get_stable(self):
        result = []
        c = self._lib.nodcy  # common block
        for i in range(c.nrnody):
            pid = self._lib.idtrafo("nxs", "pdg", c.nody[i])
            result.append(pid)
        return result

    def _generate(self):
        self._lib.aepos(-1)
        self._lib.afinal()
        self._lib.hepmcstore()
        return True
