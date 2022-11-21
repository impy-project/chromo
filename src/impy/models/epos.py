import numpy as np

from impy import impy_config
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
    _data_url = (
        "https://github.com/impy-project/impy"
        + "/releases/download/zipped_data_v1.0/epos_v001.zip"
    )

    def __init__(self, event_kinematics, seed=None, logfname=None):
        super().__init__(seed, logfname)

        self._lib.aaset(0)

        if impy_config["user_frame"] == "center-of-mass":
            iframe = 1
            self._output_frame = "center-of-mass"
        elif impy_config["user_frame"] == "laboratory":
            iframe = 2
            self._output_frame = "laboratory"

        k = event_kinematics

        datdir = _cached_data_dir(self._data_url)
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
        # EPOS decays all unstable particles by default. It uses a nodcy common block
        # to prevent decay of particles. The common block contains the array
        # nody and the length nrnody. The array holds EPOS particle ids of
        # particles that should not be decayed.

        idx = self._lib.idtrafo("pdg", "nxs", pdgid)

        c = self._lib.nodcy  # common block
        if stable:
            # append to nodcy.nody array
            if c.nrnody == len(c.nody):
                raise RuntimeError(
                    f"maximum number of stable particles reached ({c.nrnody})"
                )
            c.nody[c.nrnody] = idx
            c.nrnody += 1
        else:
            # find and remove entry from nodcy.nody array
            for i in range(c.nrnody):
                if c.nody[i] == idx:
                    c.nrnody -= 1
                    for j in range(i, c.nrnody):
                        c.nody[j] = c.nody[j + 1]
                    break

    def _get_stable(self):
        result = []
        c = self._lib.nodcy  # common block
        for i in range(c.nrnody):
            pid = self._lib.idtrafo("nxs", "pdg", c.nody[i])
            result.append(pid)
        return result

    def _generate_event(self):
        self._lib.aepos(-1)
        self._lib.afinal()
        self._lib.hepmcstore()
        return False
