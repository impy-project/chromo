"""
Created on 17.03.2014

@author: afedynitch
"""

# Some notes on UrQMD. The model works quite well now. It is slower than anything else,
# except EPOS maybe. It describes quite okaish the fixed target results and some LHC
# results. What is a strange feature is that decays of pions, kaons or neutrinos are not
# supported in the model and can not be disabled by flags.

# The current settings are taken from CORSIKA and they are optimized for speed aparently.
# The license of UrQMD is quite restrictive, they won't probably permit distributing it.

import numpy as np
from impy.common import MCRun, MCEvent, impy_config
from impy.util import info


class UrQMDEvent(MCEvent):
    """Wrapper class around EPOS particle stack."""

    def _charge_init(self, npart):
        return self._lib.uqchg.ichg[:npart]

    # Nuclear collision parameters
    @property
    def impact_parameter(self):
        """Returns impact parameter for nuclear collisions."""
        return self._lib.rsys.bimp


class UrQMD34(MCRun):
    """Implements all abstract attributes of MCRun for the
    UrQMD series of event generators.

    The version included here is UrQMD 3.4. The manual and further
    references can be accessed on this webpage https://urqmd.org/.
    """

    _name = "UrQMD"
    _version = "3.4"
    _library_name = "_urqmd34"
    _event_class = UrQMDEvent
    _output_frame = "center-of-mass"

    def __init__(self, event_kinematics, seed="random", logfname=None):
        from particletools.tables import UrQMDParticleTable

        self._stab = UrQMDParticleTable()

        super().__init__(seed, logfname)

        info(1, "First initialization")
        self._lib.urqini(self._lun, impy_config["urqmd"]["debug_level"])

        # Set default stable
        self._set_final_state_particles()

        self._lib.inputs.nevents = 1
        self._lib.rsys.bmin = 0
        # Use bdb weighting for impact parameter selection
        self._lib.options.ctoption[5 - 1] = 1

        # Disable elastic collision
        self._lib.options.ctoption[7 - 1] = 1

        # Change CTParams and/or CTOptions if needed
        if "CTParams" in impy_config["urqmd"]:
            for ctp in impy_config["urqmd"]["CTParams"]:
                self._lib.options.ctparams[ctp[0]] = ctp[1]
                info(5, "CTParams[{}] changed to {}".format(ctp[0], ctp[1]))
        if "CTOptions" in impy_config["urqmd"]:
            for cto in impy_config["urqmd"]["CTOptions"]:
                self._lib.options.ctoptions[cto[0]] = cto[1]
                info(5, "CTOptions[{}] changed to {}".format(cto[0], cto[1]))

        # Time evolution (only relevant for QMD?)
        caltim = impy_config["urqmd"]["caltim"]
        outtim = impy_config["urqmd"]["outtim"]
        # Don't do time evolution in QMD (as in CORSIKA)
        # This sets timestep to final time -> all steps are = 1
        self._lib.pots.dtimestep = outtim
        self._lib.sys.nsteps = int(0.01 + caltim / self._lib.pots.dtimestep)
        self._lib.inputs.outsteps = int(0.01 + caltim / self._lib.pots.dtimestep)
        self._set_event_kinematics(event_kinematics)

    def sigma_inel(self):
        """Inelastic cross section according to current
        event setup (energy, projectile, target)"""
        return self._lib.ptsigtot()

    def _set_event_kinematics(self, event_kinematics):
        """Set new combination of energy, momentum, projectile
        and target combination for next event."""
        k = event_kinematics
        self._curr_event_kin = k

        if not k.p1_is_nucleus:
            # Special projectile
            self._lib.inputs.prspflg = 1
            self._lib.sys.ap = 1
            self._lib.inputs.spityp[0] = self._stab.pdg2modid[k.p1pdg][0]
            self._lib.inputs.spiso3[0] = self._stab.pdg2modid[k.p1pdg][1]
        else:
            self._lib.inputs.prspflg = 0
            self._lib.sys.ap = k.A1
            self._lib.sys.zp = k.Z1

        if not k.p2_is_nucleus:
            # Special projectile
            self._lib.inputs.trspflg = 1
            self._lib.sys.at = 1
            self._lib.inputs.spityp[1] = self._stab.pdg2modid[k.p2pdg][0]
            self._lib.inputs.spiso3[1] = self._stab.pdg2modid[k.p2pdg][1]
        else:
            self._lib.inputs.trspflg = 0
            self._lib.sys.at = k.A2
            self._lib.sys.zt = k.Z2

        # Set impact parameter (to be revisited)
        # AF: Does this work for pions or kaons??
        self._lib.rsys.bdist = (
            self._lib.nucrad(self._lib.sys.ap)
            + self._lib.nucrad(self._lib.sys.at)
            + 2 * self._lib.options.ctparam[30 - 1]
        )

        # Output only correct in lab frame, don't try using the
        # "equal speed" frame.
        self._lib.input2.pbeam = k.plab
        self._lib.inputs.srtflag = 2

        # Unclear what the effect of the equation of state is but
        # should be required for very low energy.
        if k.plab > 4.9:
            self._lib.inputs.eos = 0
            # This is the fast method as set default in urqinit.f
            self._lib.options.ctoption[24 - 1] = 2
        else:
            self._lib.inputs.eos = 1
            # A hard sphere potential is required for equation of state
            self._lib.options.ctoption[24 - 1] = 0

        info(5, "Setting event kinematics")

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
        stable_ids = self._lib.stables.stabvec
        nstab = self._lib.stables.nstable
        try:
            uid = self._stab.pdg2modid[pdgid][0]
        except KeyError as e:
            e.msg += f"\nParticle {pdgid} unknown to UrQMD"
            info(5, e.msg)
            raise
        if stable:
            if uid not in stable_ids[:nstab]:
                self._lib.stables.stabvec[nstab] = uid
                self._lib.stables.nstable += 1
        else:
            if uid in stable_ids[:nstab]:
                self._lib.stables.stabvec[: nstab - 1] = np.array(
                    [pid for pid in stable_ids[:nstab] if pid != uid]
                )
                self._lib.stables.nstable -= 1

    def _generate_event(self):
        # If new event, initialize projectile and target
        if self._lib.options.ctoption[40 - 1] == 0:
            if self._lib.inputs.prspflg == 0:
                self._lib.cascinit(self._lib.sys.zp, self._lib.sys.ap, 1)
            if self._lib.inputs.trspflg == 0:
                self._lib.cascinit(self._lib.sys.zt, self._lib.sys.at, 2)

        self._lib.urqmd(0)
        # Convert URQMD event to HEPEVT
        self._lib.chepevt()
        return False
