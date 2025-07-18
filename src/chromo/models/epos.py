import numpy as np

from chromo.common import CrossSectionData, MCEvent, MCRun
from chromo.constants import standard_projectiles
from chromo.kinematics import EventFrame, GeV
from chromo.util import (
    Nuclei,
    _cached_data_dir,
    fortran_array_insert,
    fortran_array_remove,
    select_long_lived,
)


def _epos_unstable_pids():
    unknown_pids_epos = [
        117,
        119,
        217,
        219,
        227,
        229,
        317,
        319,
        327,
        329,
        337,
        1218,
        2128,
        3118,
        3128,
        3218,
        3228,
        10115,
        10215,
        10225,
        10315,
        10325,
        10335,
        13314,
        13324,
        14122,
        20315,
        20325,
        23126,
        30113,
        30213,
        30223,
        30313,
        30323,
        30443,
        100111,
    ]

    # From all unstable particles take that
    # are known to epos
    return [
        pid
        for pid in select_long_lived()
        if (abs(pid) not in unknown_pids_epos) and (abs(pid) < unknown_pids_epos[-1])
    ]


class EPOSEvent(MCEvent):
    """Wrapper class around EPOS particle stack."""

    def _get_charge(self, npart):
        return self._lib.charge_vect(self._lib.hepevt.idhep[:npart])

    def _get_impact_parameter(self):
        # return self._lib.nuc3.bimp
        return float(self._lib.cevt.bimevt)

    def _get_n_wounded(self):
        return int(self._lib.cevt.npjevt), int(self._lib.cevt.ntgevt)

    def _history_zero_indexing(self):
        # Beam particles are [-1, -1], but indexing starts from 1
        # adjust to [0, 0]
        condition = (self.mothers == [-1, -1]).all(axis=1)
        self.mothers[condition] = [0, 0]

        # Adjust to zero-base indexing
        self.mothers = self.mothers - 1

        # Set [i, i] to [i, -1]
        condition = self.mothers[:, 0] == self.mothers[:, 1]
        self.mothers[condition, 1] = -1

        self.daughters = self.daughters - 1
        condition = self.daughters[:, 0] == self.daughters[:, 1]
        self.daughters[condition, 1] = -1

    @staticmethod
    def _beam_mother_daughters_fix(event, ind):
        """
        Set correct mothers and daughters of Epos's beam particles
        after prepending projectile/target

        ind = 0 for projectile
        ind = 1 for target
        """
        # Search among epos's beam particles the particles with the
        # same energy as projectile/target
        is_parent = np.isclose(event.en, event.en[ind], rtol=1e-2) & (event.status == 4)
        # Attach them to the parent
        event.mothers[is_parent] = [ind, -1]
        # Assign daughters for projectile/target
        daughters = np.where(is_parent)[0]
        event.daughters[ind] = [np.min(daughters), np.max(daughters)]

    def _repair_initial_beam(self):
        """
        Attach spectators/wounded nucleons to corresponding nucleus
        if projectile/target is nucleus
        """
        beam = self.kin._get_beam_data(self._generator_frame)
        is_nucleus = np.abs(beam["pid"]) > 1000000000

        # Don't do anything if it's not a nucleus
        if not (np.any(is_nucleus)):
            return

        bstatus = 555  # status of prepended particles
        beam["status"][:] = bstatus
        for field, beam_field in beam.items():
            event_field = getattr(self, field)
            if np.all(is_nucleus):
                res = np.concatenate((beam_field, event_field))
            # projectile
            elif is_nucleus[0]:
                res = np.concatenate((beam_field[0:1], event_field))
            # target
            elif is_nucleus[1]:
                res = np.concatenate(
                    (event_field[0:1], beam_field[1:], event_field[1:])
                )
            setattr(self, field, res)

        shift = 0
        # if projectile is nucleus
        if is_nucleus[0]:
            self._beam_mother_daughters_fix(self, 0)
            shift += 1

        # if target is nucleus
        if is_nucleus[1]:
            self._beam_mother_daughters_fix(self, 1)
            shift += 1

        # Shift all daughters to a number of prepended particles excluding
        # prepended particles
        self.daughters[
            (self.status[:, np.newaxis] != bstatus) & (self.daughters > -1)
        ] += shift
        # Reset statuses of prepended particles
        self.status[self.status == bstatus] = 4
        # Shift all mothers of all particles excluding beam particles
        if is_nucleus[1]:
            # Attach to the projectile only.
            # We cannot attacht to, for example, [0, 2] (it occurs sometimes),
            # because we insert nucleus between 0 and 1
            # If we try we would get an error from hepmc
            self.mothers[(self.status != 4) & (self.mothers[:, 0] == 0)] = [0, -1]

        self.mothers[(self.status[:, np.newaxis] != 4) & (self.mothers > 0)] += shift


class EposLHC(MCRun):
    """Implements all abstract attributes of MCRun for the
    EPOS-LHC series of event generators."""

    _name = "EPOS"
    _version = "LHC"
    _library_name = "_eposlhc"
    _event_class = EPOSEvent
    _frame = None
    _projectiles = standard_projectiles | Nuclei()
    _unstable_pids = set(_epos_unstable_pids())
    _hadronic_rescattering = False
    _init_protection_E_A1_A2 = None
    _data_url = (
        "https://github.com/impy-project/chromo"
        "/releases/download/zipped_data_v1.0/eposlhc_v001.zip"
    )
    _ecm_min = 6 * GeV

    def __init__(self, evt_kin, *, seed=None, isigma=1):
        """
        isigma=1: cross-section is calculated by a numerical method
                  which is valid only for h-p or h-A (h being pion, kaon or nucleon)
                  but not A-B (nucleus-nucleus)
                  (not good for ionudi=2)

        isigma=0: same as isigma=1 but do not print the cross section on screen

        isigma=2: all the nuclear cross-sections are calculated by AA pseudo simulations
                  but it takes several minutes to compute
        """
        import chromo

        super().__init__(seed)

        self._lib.aaset(0)
        datdir = _cached_data_dir(self._data_url)

        lun = 6  # stdout
        self._lib.initepos(
            evt_kin.ecm,
            self._hadronic_rescattering,  # No effect on EPOS-LHC
            datdir,
            len(datdir),
            chromo.debug_level,
            lun,
            isigma,
        )

        self._set_final_state_particles()
        self._lib.charge_vect = np.vectorize(self._lib.getcharge, otypes=[np.float32])
        self._init_protection_E_A1_A2 = (evt_kin.ecm, evt_kin.p1.A, evt_kin.p2.A)
        self.kinematics = evt_kin

    def _cross_section(self, kin=None, max_info=False):
        total, inel, el, dd, sd, _ = self._lib.xsection()
        return CrossSectionData(
            total=total,
            inelastic=inel,
            prod=total
            - el
            - sd
            - dd,  # TODO: We need to ask Tanguy if this is correct to get sigprod
            elastic=el,
            diffractive_xb=sd / 2,  # this is an approximation
            diffractive_ax=sd / 2,  # this is an approximation
            diffractive_xx=dd,
            diffractive_axb=0,
        )

    def _check_init_protection(self, kin):
        Ecm_max, A1_max, A2_max = self._init_protection_E_A1_A2
        explanation = """Initialize Epos models with the highest energy and heaviest
            masses expected in the simulation, e.g.:

            kinematics = CenterOfMass(maximal_energy, (238, 92), (238, 92))
            epos = EposLHC(kinematics)

            The energy or mass can be changed during the simulation, e.g.
            epos.kinematics = CenterOfMass(100, (4, 2), (16, 8))
            for event in epos(100):
                ...
            """
        if kin.ecm > Ecm_max:
            message = "Requested energy exceeds the initialization energy: > "
            message += f"{kin.ecm} GeV > {Ecm_max} GeV.\n" + explanation
            raise ValueError(message)
        if kin.p1.A is not None and kin.p1.A > A1_max:
            message = (
                "Requested projectile mass number exceeds the initialization mass: > "
            )
            message += f"{kin.p1.A} > {A1_max}.\n" + explanation
            raise ValueError(message)
        if kin.p2.A is not None and kin.p2.A > A2_max:
            message = "Requested target mass number exceeds the initialization mass: > "
            message += f"{kin.p2.A} > {A2_max}.\n" + explanation
            raise ValueError(message)

    def _set_kinematics(self, kin):
        if kin.frame == EventFrame.FIXED_TARGET:
            iframe = 2
            self._frame = EventFrame.FIXED_TARGET
        else:
            iframe = 1
            self._frame = EventFrame.CENTER_OF_MASS
        self._check_init_protection(kin)
        self._lib.initeposevt(kin.ecm, -1.0, int(kin.p1), int(kin.p2), iframe)

    def _set_stable(self, pdgid, stable):
        # EPOS decays all unstable particles by default. It uses a nodcy common block
        # to prevent decay of particles. The common block contains the array
        # nody and the length nrnody. The array holds EPOS particle ids of
        # particles that should not be decayed.

        if pdgid not in self._unstable_pids:
            return

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

    def print_native_event(self, nparticles=100):
        self._lib.alist("EposLHC listing&", 1, nparticles)


class EposLHCR(EposLHC):
    """CRMC 2.2.0 version of EPOS LHC-R. No explicit hadronic rescattering."""

    _name = "EPOS"
    _version = "LHC-R"
    _library_name = "_eposlhcr"
    _event_class = EPOSEvent
    _data_url = (
        "https://github.com/impy-project/chromo"
        "/releases/download/zipped_data_v1.0/eposlhcr_v001.zip"
    )
    _hadronic_rescattering = False

    def _generate(self):
        self._lib.aepos(-1)
        self._lib.afinal()
        self._lib.hepmcstore(-1)
        return True


class EposLHCRHadrRescattering(EposLHCR):
    """CRMC 2.2.0 version of EPOS LHC-R."""

    _version = "LHC-R hadr. rescattering"
    _hadronic_rescattering = True
