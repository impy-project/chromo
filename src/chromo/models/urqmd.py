# Some notes on UrQMD. The model works quite well now. It is slower than anything else,
# except EPOS maybe. It describes quite okaish the fixed target results and some LHC
# results. What is a strange feature is that decays of pions, kaons or neutrinos are not
# supported in the model and can not be disabled by flags.

# The current settings are taken from CORSIKA and they are optimized for speed aparently.
# The license of UrQMD is quite restrictive, they won't probably permit distributing it.

from chromo.common import MCRun, MCEvent, CrossSectionData
from chromo.util import (
    info,
    fortran_array_insert,
    fortran_array_remove,
    Nuclei,
)
from chromo.kinematics import EventFrame
from chromo.constants import standard_projectiles, GeV
import warnings


class UrQMDEvent(MCEvent):
    """Wrapper class around EPOS particle stack."""

    def _charge_init(self, npart):
        return self._lib.uqchg.ichg[:npart]

    def _get_impact_parameter(self):
        return self._lib.rsys.bimp

    def _history_zero_indexing(self):
        # Urqmd produces wrong history
        self.mothers[:] = [-1, -1]
        self.daughters[:] = [-1, -1]

    def _repair_initial_beam(self):
        self._prepend_initial_beam()
        # Repair history
        self.mothers[(self.mothers == [1, 1]).all(axis=1)] = [0, 1]
        # Set [i, i] to [i, -1]
        condition = self.mothers[:, 0] == self.mothers[:, 1]
        self.mothers[condition, 1] = -1
        # No daughters
        self.daughters[:] = [-1, -1]


_urqmd_unstable_pids = set(
    [
        111,
        113,
        211,
        -211,
        213,
        -213,
        221,
        223,
        313,
        -313,
        321,
        -321,
        323,
        -323,
        333,
        411,
        -411,
        413,
        -413,
        421,
        -421,
        431,
        -431,
        441,
        443,
        1114,
        -1114,
        2112,
        -2112,
        2114,
        -2114,
        2214,
        -2214,
        2224,
        -2224,
        3112,
        -3112,
        3114,
        -3114,
        3122,
        -3122,
        3212,
        -3212,
        3214,
        -3214,
        3222,
        -3222,
        3224,
        -3224,
        3312,
        -3312,
        3314,
        -3314,
        3322,
        -3322,
        3324,
        -3324,
        3334,
        -3334,
        10421,
        -10421,
    ]
)


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
    _frame = EventFrame.FIXED_TARGET
    _unstable_pids = _urqmd_unstable_pids
    _projectiles = standard_projectiles | Nuclei()
    _targets = Nuclei()
    _ecm_min = 2 * GeV

    def __init__(
        self,
        evt_kin,
        *,
        seed=None,
        caltim=200,
        outtim=200,
        ct_params=None,
        ct_options=None,
    ):
        import chromo

        self._pdg2modid = {
            22: (100, 0),
            111: (101, 0),
            211: (101, 2),
            -211: (101, -2),
            2212: (1, 1),
            2112: (1, -1),
            221: (102, 0),
            223: (103, 0),
            213: (104, 2),
            -213: (104, -2),
            113: (104, 0),
            321: (106, 1),
            -321: (-106, -1),
            311: (106, -1),
            -311: (-106, 1),
            441: (107, 0),
            323: (108, 2),
            -323: (108, -2),
            313: (108, 0),
            -313: (-108, 0),
            333: (109, 0),
            3222: (40, 2),
            3212: (40, 0),
            3112: (40, -2),
            3322: (49, 0),
            3312: (49, -1),
            3122: (27, 0),
            2224: (17, 4),
            2214: (17, 2),
            2114: (17, 0),
            1114: (17, -2),
            3224: (41, 2),
            3214: (41, 0),
            3114: (41, -2),
            3324: (50, 0),
            3314: (50, -1),
            3334: (55, 0),
            411: (133, 2),
            -411: (133, -2),
            421: (133, 0),
            -421: (-133, 0),
            431: (138, 1),
            -431: (138, -1),
            433: (139, 1),
            -433: (139, -1),
            413: (134, 1),
            -413: (134, -1),
            10421: (134, 0),
            -10421: (-134, 0),
            443: (135, 0),
            -2212: (-1, 1),
            -2112: (-1, -1),
            -3222: (-40, 2),
            -3212: (-40, 0),
            -3112: (-40, -2),
            -3322: (-49, 0),
            -3312: (-49, -1),
            -3122: (-27, 0),
            -2224: (-17, 4),
            -2214: (-17, 2),
            -2114: (-17, 0),
            -1114: (-17, -2),
            -3224: (-41, 2),
            -3214: (-41, 0),
            -3114: (-41, -2),
            -3324: (-50, 0),
            -3314: (-50, -1),
            -3334: (-55, 0),
        }

        super().__init__(seed)

        # logging
        lun = 6  # stdout
        self._lib.urqini(lun, chromo.debug_level)

        self._lib.inputs.nevents = 1
        self._lib.rsys.bmin = 0
        # Use bdb weighting for impact parameter selection
        self._lib.options.ctoption[5 - 1] = 1

        # Disable elastic collision
        self._lib.options.ctoption[7 - 1] = 1

        # Change CTParams and/or CTOptions if needed
        if ct_params:
            for k, v in ct_params:
                self._lib.options.ctparams[k] = v
                info(5, f"CTParams[{k}] changed to {v}")
        if ct_options:
            for k, v in ct_options:
                self._lib.options.ctoptions[k] = v
                info(5, f"CTOptions[{k}] changed to {v}")

        # Don't do time evolution in QMD (as in CORSIKA)
        # This sets timestep to final time -> all steps are = 1
        self._lib.pots.dtimestep = outtim
        self._lib.sys.nsteps = int(0.01 + caltim / self._lib.pots.dtimestep)
        self._lib.inputs.outsteps = int(0.01 + caltim / self._lib.pots.dtimestep)
        self.kinematics = evt_kin

        self._set_final_state_particles()

    def _cross_section(self, kin=None, max_info=False):
        tot = self._lib.ptsigtot()
        return CrossSectionData(total=tot)

    def _set_kinematics(self, kin):
        for i, (p, x) in enumerate(zip((kin.p1, kin.p2), "pt")):
            if p.is_nucleus:
                setattr(self._lib.inputs, f"{x}rspflg", 0)
                setattr(self._lib.sys, f"a{x}", p.A)
                setattr(self._lib.sys, f"z{x}", p.Z)
            else:
                # Special projectile
                setattr(self._lib.inputs, f"{x}rspflg", 1)
                setattr(self._lib.sys, f"a{x}", 1)
                self._lib.inputs.spityp[i] = self._pdg2modid[p][0]
                self._lib.inputs.spiso3[i] = self._pdg2modid[p][1]

        # Set impact parameter (to be revisited)
        # AF: Does this work for pions or kaons??
        self._lib.rsys.bdist = (
            self._lib.nucrad(self._lib.sys.ap)
            + self._lib.nucrad(self._lib.sys.at)
            + 2 * self._lib.options.ctparam[30 - 1]
        )

        # Output only correct in lab frame, don't try using the
        # "equal speed" frame.
        self._lib.input2.pbeam = kin.plab
        self._lib.inputs.srtflag = 2

        # Unclear what the effect of the equation of state is but
        # should be required for very low energy.
        if kin.plab > 4.9:
            self._lib.inputs.eos = 0
            # This is the fast method as set default in urqinit.f
            self._lib.options.ctoption[24 - 1] = 2
        else:
            self._lib.inputs.eos = 1
            # A hard sphere potential is required for equation of state
            self._lib.options.ctoption[24 - 1] = 0

    def _set_stable(self, pdgid, stable):
        if pdgid not in self._unstable_pids:
            return

        try:
            uid = self._pdg2modid[pdgid][0]
        except KeyError:
            warnings.warn(f"{pdgid} unknown to UrQMD", RuntimeWarning)
            return

        # FIXME changing stabvec has no effect
        s = self._lib.stables
        if stable:
            fortran_array_insert(s.stabvec, s.nstable, uid)
        else:
            fortran_array_remove(s.stabvec, s.nstable, uid)

    def _generate(self):
        # If new event, initialize projectile and target
        if self._lib.options.ctoption[40 - 1] == 0:
            if self._lib.inputs.prspflg == 0:
                self._lib.cascinit(self._lib.sys.zp, self._lib.sys.ap, 1)
            if self._lib.inputs.trspflg == 0:
                self._lib.cascinit(self._lib.sys.zt, self._lib.sys.at, 2)

        self._lib.urqmd(0)
        # Convert URQMD event to HEPEVT
        self._lib.chepevt()
        return True
