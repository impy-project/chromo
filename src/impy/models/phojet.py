"""
Created on 15.05.2012

@author: afedynitch
"""
import numpy as np
from impy.common import MCRun, MCEvent
from impy import impy_config
from impy.util import info, fortran_chars


class PhojetEvent(MCEvent):
    """Wrapper class around (pure) PHOJET particle stack."""

    def __init__(self, lib, event_kinematics, event_frame):
        # HEPEVT (style) common block
        evt = lib.poevt1

        # Save selector for implementation of on-demand properties
        px, py, pz, en, m = evt.phep
        vx, vy, vz, vt = evt.vhep

        MCEvent.__init__(
            self,
            lib=lib,
            event_kinematics=event_kinematics,
            event_frame=event_frame,
            nevent=evt.nevhep,
            npart=evt.nhep,
            p_ids=evt.idhep,
            status=evt.isthep,
            px=px,
            py=py,
            pz=pz,
            en=en,
            m=m,
            vx=vx,
            vy=vy,
            vz=vz,
            vt=vt,
            pem_arr=evt.phep,
            vt_arr=evt.vhep,
        )

    def filter_final_state(self):
        self.selection = np.where(self.status == 1)
        self._apply_slicing()

    def filter_final_state_charged(self):
        self.selection = np.where((self.status == 1) & (self.charge != 0))
        self._apply_slicing()

    @property
    def _charge_init(self):
        return self.lib.poevt2.icolor[0, self.selection] // 3

    def _gen_cut_info(self):
        """Init variables tracking the number of soft and hard cuts"""

        self._mpi = self.lib.podebg.kspom + self.lib.podebg.khpom
        self._kspom = self.lib.podebg.kspom
        self._khpom = self.lib.podebg.khpom
        self._ksoft = self.lib.podebg.ksoft
        self._khard = self.lib.podebg.khard

    @property
    def mpi(self):
        """Total number of cuts"""
        return self.lib.podebg.kspom + self.lib.podebg.khpom

    @property
    def kspom(self):
        """Total number of soft cuts"""
        return self.lib.podebg.kspom

    @property
    def khpom(self):
        """Total number of hard cuts"""
        return self.lib.podebg.khpom

    @property
    def ksoft(self):
        """Total number of realized soft cuts"""
        return self.lib.podebg.ksoft

    @property
    def khard(self):
        """Total number of realized hard cuts"""
        return self.lib.podebg.khard

    @property
    def parents(self):
        MCEvent.parents(self)
        return self.lib.poevt1.jmohep

    @property
    def children(self):
        MCEvent.children(self)
        return self.lib.poevt1.jdahep

    # def elastic_t(self):
    #     """Squared momentum transfer t for elastic interaction.

    #     This only makes sense if the interaction is indeed elastic
    #     and only 4 particles are on the stack. The initial 2 and
    #     the final 2. Handle with care!!
    #     """
    #     return np.sum(
    #         np.square(self.lib.poevt1.phep[0:4, 0] -
    #                   self.lib.poevt1.phep[0:4, 5]))


class PHOJETRun(MCRun):
    """Implements all abstract attributes of MCRun for the
    PHOJET series of event generators.

    PHOJET is part of DPMJET and is run via the same fortran library.
    The results for pp and hadron-p `should` be the same.
    """

    def sigma_inel(self, *args, **kwargs):
        """Inelastic cross section according to current
        event setup (energy, projectile, target)"""

        info(3, "PHOJET workaround for cross-section", "(need to generate dummy event)")
        self.lib.pho_event(1, self._curr_event_kin.p1pdg, self._curr_event_kin.p2pdg)
        return self.lib.powght.siggen[3]

    def sigma_inel_air(self, **kwargs):
        """PHOJET does not support nuclear targets."""
        raise Exception("PHOJET does not support nuclear targets.")

    def set_stable(self, pdgid, stable=True):
        if abs(pdgid) == 2212:
            return
        kc = self.lib.pycomp(pdgid)
        if stable:
            self.lib.pydat3.mdcy[kc - 1, 0] = 0
            info(5, "defining", pdgid, "as stable particle")
        else:
            self.lib.pydat3.mdcy[kc - 1, 0] = 1
            info(5, "forcing decay of", pdgid)

    def _set_event_kinematics(self, event_kinematics):
        """Set new combination of energy, momentum, projectile
        and target combination for next event."""
        info(5, "Setting event kinematics")

        self._curr_event_kin = event_kinematics
        k = event_kinematics

        if k.p1_is_nucleus or k.p2_is_nucleus:
            raise Exception("PHOJET does not support nuclei.")

        self.lib.pho_setpar(1, k.p1pdg, 0, 0.0)
        self.lib.pho_setpar(2, k.p2pdg, 0, 0.0)
        self.p1, self.p2 = k.beam_as_4vec

    def attach_log(self, fname=None):
        """Routes the output to a file or the stdout."""
        fname = impy_config["output_log"] if fname is None else fname
        if fname == "stdout":
            lun = 6
            info(5, "Output is routed to stdout.")
        else:
            lun = self._attach_fortran_logfile(fname)
            info(5, "Output is routed to", fname, "via LUN", lun)

        if hasattr(self.lib, "dtflka"):
            self.lib.dtflka.lout = lun
            self.lib.dtflka.lpri = 50
        elif hasattr(self.lib, "dtiont"):
            self.lib.dtiont.lout = lun
        else:
            raise Exception(
                "Unknown PHOJET (DPMJET) version, IO common block not detected."
            )

        self.lib.pydat1.mstu[10] = lun

        # Save lun for initialization of PHOJET
        self._lun = lun

    def init_generator(self, event_kinematics, seed="random", logfname=None):
        from impy.constants import c
        from random import randint

        self._abort_if_already_initialized()

        if seed == "random":
            seed = randint(1000000, 10000000)
        else:
            seed = int(seed)
        info(3, "Using seed:", seed)

        # Detect what kind of PHOJET interface is attached. If PHOJET
        # is run through DPMJET, initial init needs -2 else -1
        init_flag = -2 if "dpmjetIII" in self.lib.__name__ else -1

        pho_conf = impy_config["phojet"]
        # Set the dpmjpar.dat file
        if hasattr(self.lib, "pomdls") and hasattr(self.lib.pomdls, "parfn"):
            pfile = pho_conf["param_file"][self.version]
            info(3, "PHOJET parameter file at", pfile)
            self.lib.pomdls.parfn = fortran_chars(self.lib.pomdls.parfn, pfile)

        # Set the data directory for the other files
        if hasattr(self.lib, "poinou") and hasattr(self.lib.poinou, "datdir"):
            pfile = pho_conf["dat_dir"][self.version]
            info(3, "PHOJET data dir is at", pfile)
            self.lib.poinou.datdir = fortran_chars(self.lib.poinou.datdir, pfile)
            self.lib.poinou.lendir = len(pfile)

        # Set debug level of the generator
        for i in range(self.lib.podebg.ideb.size):
            self.lib.podebg.ideb[i] = pho_conf["debug_level"]

        self.attach_log(fname=logfname)

        # Initialize PHOJET's parameters
        if self.lib.pho_init(init_flag, self._lun):
            raise Exception("PHOJET unable to initialize or set LUN")

        process_switch = self.lib.poprcs.ipron
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
        if self.lib.pho_event(-1, self.p1, self.p2)[1]:
            raise Exception(
                "PHOJET failed to initialize with the current", "event kinematics"
            )

        # Set seed of random number generator
        sseed = str(seed)
        n1, n2, n3, n4 = (
            int(sseed[0:2]),
            int(sseed[2:4]),
            int(sseed[4:6]),
            int(sseed[6:]),
        )
        self.lib.dt_rndmst(n1, n2, n3, n4)

        # if self.def_settings:
        #     print self.class_name + \
        #         "::init_generator(): Using default settings:", \
        #         self.def_settings.__class__.__name__
        #     self.def_settings.enable()

        self._define_default_fs_particles()
        # Set PYTHIA decay flags to follow all changes to MDCY
        self.lib.pydat1.mstj[21 - 1] = 1
        self.lib.pydat1.mstj[22 - 1] = 2
        # Set ctau threshold in PYTHIA for the default stable list
        self.lib.pydat1.parj[70] = impy_config["tau_stable"] * c * 1e-3  # mm

    def generate_event(self):
        return self.lib.pho_event(1, self.p1, self.p2)[1]


class Phojet112(PHOJETRun):
    def __init__(self, event_kinematics, seed="random", logfname=None):
        from impy.definitions import interaction_model_by_tag as models_dict

        interaction_model_def = models_dict["PHOJET112"]
        super().__init__(interaction_model_def)
        self.init_generator(event_kinematics, seed, logfname)


class Phojet191(PHOJETRun):
    def __init__(self, event_kinematics, seed="random", logfname=None):
        from impy.definitions import interaction_model_by_tag as models_dict

        interaction_model_def = models_dict["PHOJET191"]
        super().__init__(interaction_model_def)
        self.init_generator(event_kinematics, seed, logfname)
