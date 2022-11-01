import numpy as np

from impy import base_path, impy_config
from impy.common import MCEvent, MCRun
from impy.util import info, _cached_data_dir


class EPOSEvent(MCEvent):
    """Wrapper class around EPOS particle stack."""

    def __init__(self, generator):
        super().__init__(generator)
        # EPOS sets parents of beam particles to (-1, -1).
        # We change it to (0, 0)
        self.parents[self.status == 4] = 0

    def _charge_init(self, npart):
        return self._lib.charge_vect(self._lib.hepevt.idhep[:npart])

    # Nuclear collision parameters
    @property
    def impact_parameter(self):
        """Returns impact parameter for nuclear collisions."""
        # return self._lib.nuc3.bimp
        return self._lib.cevt.bimevt

    @property
    def n_wounded_A(self):
        """Number of wounded nucleons side A"""
        return self._lib.cevt.npjevt

    @property
    def n_wounded_B(self):
        """Number of wounded nucleons (target) side B"""
        return self._lib.cevt.ntgevt

    @property
    def n_wounded(self):
        """Number of total wounded nucleons"""
        return self.n_wounded_A + self.n_wounded_B

    @property
    def n_spectator_A(self):
        """Number of spectator nucleons side A"""
        return self._lib.cevt.npnevt + self._lib.cevt.nppevt

    @property
    def n_spectator_B(self):
        """Number of spectator nucleons (target) side B"""
        return self._lib.cevt.ntnevt + self._lib.cevt.ntpevt


class EposLHC(MCRun):
    """Implements all abstract attributes of MCRun for the
    EPOS-LHC series of event generators."""

    _name = "EPOS"
    _version = "LHC"
    _event_class = EPOSEvent
    _library_name = "_eposlhc"
    _output_frame = "center-of-mass"

    def __init__(self, event_kinematics, seed=None, logfname=None):
        from os import path

        super().__init__(seed, logfname)

        epos_conf = impy_config["epos"]

        info(1, "First initialization")
        self._lib.aaset(0)

        if impy_config["user_frame"] == "center-of-mass":
            iframe = 1
            self._output_frame = "center-of-mass"
        elif impy_config["user_frame"] == "laboratory":
            iframe = 2
            self._output_frame = "laboratory"

        k = event_kinematics

        _cached_data_dir(
            "https://github.com/impy-project/impy/releases/download"
            "/zipped_data_v1.0/epos_v001.zip"
        )
        datdir = path.join(base_path, epos_conf["datdir"])
        self._lib.initializeepos(
            float(self._seed),
            k.ecm,
            datdir,
            len(datdir),
            iframe,
            k.p1pdg,
            k.p2pdg,
            impy_config["epos"]["debug_level"],
            self._lun,
        )

        # Set default stable
        self._set_final_state_particles()
        self._lib.charge_vect = np.vectorize(self._lib.getcharge, otypes=[np.float32])
        self.event_kinematics = event_kinematics

    def _sigma_inel(self, evt_kin):
        with self._temporary_evt_kin(evt_kin):
            return self._lib.xsection()[1]

    def _set_event_kinematics(self, k):
        info(5, "Setting event kinematics")
        self._lib.initeposevt(k.ecm, -1.0, k.p1pdg, k.p2pdg)

    def _attach_log(self, fname=None):
        """Routes the output to a file or the stdout."""
        fname = impy_config["output_log"] if fname is None else fname
        if fname == "stdout":
            lun = 6
            info(5, "Output is routed to stdout.")
        else:
            lun = self._attach_fortran_logfile(fname)
            info(5, "Output is routed to", fname, "via LUN", lun)

        self._lun = lun

    def _set_stable(self, pdgid, stable):
        if stable:
            self._lib.setstable(pdgid)
        else:
            self._lib.setunstable(pdgid)

    def _generate_event(self):
        self._lib.aepos(-1)
        self._lib.afinal()
        self._lib.hepmcstore()
        return False
