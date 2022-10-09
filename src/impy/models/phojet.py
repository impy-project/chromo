"""
Created on 15.05.2012

@author: afedynitch
"""
from impy.common import MCRun, MCEvent
from impy import impy_config
from impy.util import info, fortran_chars


class PhojetEvent(MCEvent):
    """Wrapper class around (pure) PHOJET particle stack."""

    _hepevt = "poevt1"

    def _charge_init(self, npart):
        return self._lib.poevt2.icolor[0, :npart] / 3

    def _gen_cut_info(self):
        """Init variables tracking the number of soft and hard cuts"""

        self._mpi = self._lib.podebg.kspom + self._lib.podebg.khpom
        self._kspom = self._lib.podebg.kspom
        self._khpom = self._lib.podebg.khpom
        self._ksoft = self._lib.podebg.ksoft
        self._khard = self._lib.podebg.khard

    @property
    def mpi(self):
        """Total number of cuts"""
        return self._lib.podebg.kspom + self._lib.podebg.khpom

    @property
    def kspom(self):
        """Total number of soft cuts"""
        return self._lib.podebg.kspom

    @property
    def khpom(self):
        """Total number of hard cuts"""
        return self._lib.podebg.khpom

    @property
    def ksoft(self):
        """Total number of realized soft cuts"""
        return self._lib.podebg.ksoft

    @property
    def khard(self):
        """Total number of realized hard cuts"""
        return self._lib.podebg.khard

    # def elastic_t(self):
    #     """Squared momentum transfer t for elastic interaction.

    #     This only makes sense if the interaction is indeed elastic
    #     and only 4 particles are on the stack. The initial 2 and
    #     the final 2. Handle with care!!
    #     """
    #     return np.sum(
    #         np.square(self._lib.poevt1.phep[0:4, 0] -
    #                   self._lib.poevt1.phep[0:4, 5]))


class PHOJETRun(MCRun):
    """Implements all abstract attributes of MCRun for the
    PHOJET series of event generators.

    PHOJET is part of DPMJET and is run via the same fortran library.
    The results for pp and hadron-p `should` be the same.
    """

    _name = "PhoJet"
    _event_class = PhojetEvent
    _output_frame = "center-of-mass"

    def __init__(self, event_kinematics, seed="random", logfname=None):
        from impy.constants import c

        super().__init__(seed, logfname)

        # Detect what kind of PHOJET interface is attached. If PHOJET
        # is run through DPMJET, initial init needs -2 else -1
        init_flag = -2 if "dpmjetIII" in self._lib.__name__ else -1

        pho_conf = impy_config["phojet"]
        # Set the dpmjpar.dat file
        if hasattr(self._lib, "pomdls") and hasattr(self._lib.pomdls, "parfn"):
            pfile = pho_conf["param_file"][self.version]
            info(3, "PHOJET parameter file at", pfile)
            self._lib.pomdls.parfn = fortran_chars(self._lib.pomdls.parfn, pfile)

        # Set the data directory for the other files
        if hasattr(self._lib, "poinou") and hasattr(self._lib.poinou, "datdir"):
            pfile = pho_conf["dat_dir"][self.version]
            info(3, "PHOJET data dir is at", pfile)
            self._lib.poinou.datdir = fortran_chars(self._lib.poinou.datdir, pfile)
            self._lib.poinou.lendir = len(pfile)

        # Set debug level of the generator
        for i in range(self._lib.podebg.ideb.size):
            self._lib.podebg.ideb[i] = pho_conf["debug_level"]

        # Initialize PHOJET's parameters
        if self._lib.pho_init(init_flag, self._lun):
            raise Exception("PHOJET unable to initialize or set LUN")

        process_switch = self._lib.poprcs.ipron
        # non-diffractive elastic scattering (1 - on, 0 - off)
        process_switch[0, 0] = 1
        # elastic scattering
        process_switch[1, 0] = 0
        # quasi-elastic scattering (for incoming photons only)
        process_switch[2, 0] = 1
        # central diffration (double-pomeron scattering)
        process_switch[3, 0] = 1
        # particle 1 single diff. dissociation
        process_switch[4, 0] = 1
        # particle 2 single diff. dissociation
        process_switch[5, 0] = 1
        # double diff. dissociation
        process_switch[6, 0] = 1
        # direct photon interaction (for incoming photons only)
        process_switch[7, 0] = 1

        self._set_event_kinematics(event_kinematics)
        if self._lib.pho_event(-1, self.p1, self.p2)[1]:
            raise Exception(
                "PHOJET failed to initialize with the current", "event kinematics"
            )

        # if self.def_settings:
        #     print self.class_name + \
        #         "::init_generator(): Using default settings:", \
        #         self.def_settings.__class__.__name__
        #     self.def_settings.enable()

        self._set_final_state_particles()
        # Set PYTHIA decay flags to follow all changes to MDCY
        self._lib.pydat1.mstj[21 - 1] = 1
        self._lib.pydat1.mstj[22 - 1] = 2
        # Set ctau threshold in PYTHIA for the default stable list
        self._lib.pydat1.parj[70] = impy_config["tau_stable"] * c * 1e-3  # mm

    def sigma_inel(self):
        """Inelastic cross section according to current
        event setup (energy, projectile, target)"""

        info(3, "PHOJET workaround for cross-section", "(need to generate dummy event)")
        self._lib.pho_event(1, self._curr_event_kin.p1pdg, self._curr_event_kin.p2pdg)
        return self._lib.powght.siggen[3]

    def sigma_inel_air(self, **kwargs):
        """PHOJET does not support nuclear targets."""
        raise Exception("PHOJET does not support nuclear targets.")

    def _set_stable(self, pdgid, stable):
        if abs(pdgid) == 2212:
            return
        kc = self._lib.pycomp(pdgid)
        self._lib.pydat3.mdcy[kc - 1, 0] = not stable

    def _set_event_kinematics(self, event_kinematics):
        """Set new combination of energy, momentum, projectile
        and target combination for next event."""
        info(5, "Setting event kinematics")

        self._curr_event_kin = event_kinematics
        k = event_kinematics

        if k.p1_is_nucleus or k.p2_is_nucleus:
            raise Exception("PHOJET does not support nuclei.")

        self._lib.pho_setpar(1, k.p1pdg, 0, 0.0)
        self._lib.pho_setpar(2, k.p2pdg, 0, 0.0)
        self.p1, self.p2 = k.beam_as_4vec

    def _attach_log(self, fname=None):
        """Routes the output to a file or the stdout."""
        fname = impy_config["output_log"] if fname is None else fname
        if fname == "stdout":
            lun = 6
            info(5, "Output is routed to stdout.")
        else:
            lun = self._attach_fortran_logfile(fname)
            info(5, "Output is routed to", fname, "via LUN", lun)

        if hasattr(self._lib, "dtflka"):
            self._lib.dtflka.lout = lun
            self._lib.dtflka.lpri = 50
        elif hasattr(self._lib, "dtiont"):
            self._lib.dtiont.lout = lun
        elif hasattr(self._lib, "poinou"):
            self._lib.poinou.lo = lun
        else:
            raise Exception(
                "Unknown PHOJET (DPMJET) version, IO common block not detected."
            )

        self._lib.pydat1.mstu[10] = lun

        # Save lun for initialization of PHOJET
        self._lun = lun

    def _generate_event(self):
        return self._lib.pho_event(1, self.p1, self.p2)[1]


class Phojet112(PHOJETRun):
    _version = "1.12-35"
    _library_name = "_phojet112"


class Phojet191(PHOJETRun):
    _version = "19.1"
    _library_name = "_phojet191"


class Phojet193(PHOJETRun):
    _version = "19.3"
    _library_name = "_phojet193"
