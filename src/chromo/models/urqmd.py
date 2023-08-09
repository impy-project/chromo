# Some notes on UrQMD. The model works quite well now. It is slower than anything else,
# except EPOS maybe. It describes quite okaish the fixed target results and some LHC
# results. What is a strange feature is that decays of pions, kaons or neutrinos are not
# supported in the model and can not be disabled by flags.

# The current settings are taken from CORSIKA and they are optimized for speed aparently.
# The license of UrQMD is quite restrictive, they won't probably permit distributing it.

from chromo.common import MCRun, MCEvent, CrossSectionData
from chromo.util import info, fortran_array_insert, fortran_array_remove, Nuclei
from chromo.kinematics import EventFrame
from chromo.constants import standard_projectiles, GeV, all_decaying_pids
import warnings


class UrQMDEvent(MCEvent):
    """Wrapper class around EPOS particle stack."""

    def _charge_init(self, npart):
        return self._lib.uqchg.ichg[:npart]

    def _get_impact_parameter(self):
        return self._lib.rsys.bimp


def get_urqmd_decaying_pids():
    unknown_pids_urqmd = [
        6,
        -6,
        13,
        -13,
        15,
        -15,
        23,
        24,
        -24,
        25,
        115,
        117,
        119,
        130,
        215,
        -215,
        217,
        -217,
        219,
        -219,
        225,
        227,
        229,
        310,
        315,
        -315,
        317,
        -317,
        319,
        -319,
        325,
        -325,
        327,
        -327,
        329,
        -329,
        331,
        335,
        337,
        415,
        -415,
        423,
        -423,
        425,
        -425,
        435,
        -435,
        445,
        511,
        -511,
        515,
        -515,
        521,
        -521,
        525,
        -525,
        531,
        -531,
        535,
        -535,
        541,
        -541,
        551,
        553,
        1112,
        -1112,
        1116,
        -1116,
        1118,
        -1118,
        1212,
        -1212,
        1214,
        -1214,
        1216,
        -1216,
        1218,
        -1218,
        2116,
        -2116,
        2118,
        -2118,
        2122,
        -2122,
        2124,
        -2124,
        2126,
        -2126,
        2128,
        -2128,
        2216,
        -2216,
        2218,
        -2218,
        2222,
        -2222,
        2226,
        -2226,
        2228,
        -2228,
        3116,
        -3116,
        3118,
        -3118,
        3124,
        -3124,
        3126,
        -3126,
        3128,
        -3128,
        3216,
        -3216,
        3218,
        -3218,
        3226,
        -3226,
        3228,
        -3228,
        4112,
        -4112,
        4114,
        -4114,
        4122,
        -4122,
        4132,
        -4132,
        4212,
        -4212,
        4214,
        -4214,
        4222,
        -4222,
        4224,
        -4224,
        4232,
        -4232,
        4314,
        -4314,
        4324,
        -4324,
        4332,
        -4332,
        5112,
        -5112,
        5114,
        -5114,
        5122,
        -5122,
        5132,
        -5132,
        5222,
        -5222,
        5224,
        -5224,
        5232,
        -5232,
        5332,
        -5332,
        10111,
        10113,
        10115,
        10211,
        -10211,
        10213,
        -10213,
        10215,
        -10215,
        10221,
        10223,
        10225,
        10311,
        -10311,
        10313,
        -10313,
        10315,
        -10315,
        10321,
        -10321,
        10323,
        -10323,
        10325,
        -10325,
        10331,
        10333,
        10335,
        10411,
        -10411,
        10413,
        -10413,
        10423,
        -10423,
        10431,
        -10431,
        10433,
        -10433,
        10441,
        10443,
    ]

    decaying_pids = []
    for pid in all_decaying_pids:
        if (abs(pid) not in unknown_pids_urqmd) and (abs(pid) < unknown_pids_urqmd[-1]):
            decaying_pids.append(pid)

    return decaying_pids


urqmd_decaying_pids = get_urqmd_decaying_pids()


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
    _decaying_pids = urqmd_decaying_pids
    _projectiles = standard_projectiles | Nuclei()
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

    def _cross_section(self, kin=None):
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
        if pdgid not in self._decaying_pids:
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
