import warnings

from particle import PDGID, Particle

from chromo.common import CrossSectionData, MCEvent, MCRun
from chromo.constants import GeV, standard_projectiles
from chromo.kinematics import EventFrame
from chromo.util import (
    Nuclei,
    _cached_data_dir,
    fortran_chars,
    info,
    pdg2AZ,
)

# The list below are all known particles to DPMJET, which can be used as
# projectiles. To generate this list, run the following script:
# ```python
# from particle import Particle
# projectile_list = Particle.findall(
#     lambda p: not (
#         abs(p.pdgid) < 11
#         or abs(p.pdgid) in [21, 22]
#         or (22 < abs(p.pdgid) < 111)
#         or abs(p.pdgid) > 5000
#         or p.pdgid.has_bottom
#         or (p.mass is None)
#         or (p.pdgid.is_lepton)
#         or dpm._lib.idt_icihad(p.pdgid) == 0
#     ),
#     particle=True,
# )
# print(set([int(p.pdgid) for p in projectile_list]))
# ```
# fmt: off
dpmjet_extended_projectiles = {
    130, 3334, 4232, 3212, 3214, 3216, 4112, 3218, 3222, 3224,
    4122, 411, 413, 2212, 421, 2214, 423, 3112, 4132, 3114,
    431, 2224, 433, 3122, 310, 311, 313, 441, 443, 2112,
    321, 2114, 323, 331, 333, 211, 213, 1114, 221, 223,
    111, 3312, 113, 3314, 4212, 3322, 3324, 4222
}
dpmjet_extended_projectiles = {Particle.from_pdgid(p).pdgid for p in dpmjet_extended_projectiles}
# fmt: on


class DpmjetIIIEvent(MCEvent):
    """Wrapper class around DPMJET-III HEPEVT-style particle stack."""

    _hepevt = "dtevt1"
    _phep = "phkk"
    _vhep = "vhkk"
    _nevhep = "nevhkk"
    _nhep = "nhkk"
    _idhep = "idhkk"
    _isthep = "isthkk"
    _jmohep = "jmohkk"
    _jdahep = "jdahkk"

    def _get_charge(self, npart):
        return self._lib.dtpart.iich[self._lib.dtevt2.idbam[:npart] - 1]

    def _get_impact_parameter(self):
        return self._lib.dtglcp.bimpac

    def _get_n_wounded(self):
        return self._lib.dtglcp.nwasam, self._lib.dtglcp.nwbsam

    def _repair_initial_beam(self):
        beam = self.kin._get_beam_data(self._generator_frame)
        for field in ["pid", "status", "charge", "px", "py", "pz", "en", "m"]:
            event_field = getattr(self, field)
            event_field[0:2] = beam[field]

    def _prepare_for_hepmc(self):
        model, version = self.generator
        warnings.warn(
            f"{model}-{version}: only part of the history available in HepMC3 event",
            RuntimeWarning,
        )
        mask = (
            (self.status == 1)
            | (self.status == 2)
            | (self.status == 4)
            | (self.pid == 99999)
        )
        return self[mask]

    # Unfortunately not that simple since this is bounced through
    # entire code as argument not in COMMON
    # @property
    # def n_inel_NN_interactions(self):
    #     """Number of inelastic nucleon-nucleon interactions"""
    #     return self._lib.dtglcp.nwtsum


# =========================================================================
# DpmjetIIIMCRun
# =========================================================================
class DpmjetIIIRun(MCRun):
    """Implements all abstract attributes of MCRun for the
    DPMJET-III series of event generators.

    It should work identically for the new 'dpmjet3' module and the legacy
    dpmjet307. No special constructor is necessary and everything is
    handled by the default constructor of the base class.

    Notes
    -----
    Initialize at the highest energy of the simulation: the PHOJET
    hadron-nucleon tables (``DT_INIT``) end there, and higher-energy
    requests raise ``ValueError``.

    For cross-section tabulation use a fresh instance per point: the
    Glauber sigma drifts over many kinematics switches (19.3, p-air
    100 GeV: sigma_prod 278 mb on first query, ~258 mb late in a loop).
    """

    _name = "DPMJET-III"
    _event_class = DpmjetIIIEvent
    _frame = None
    # Photon projectiles on nuclear targets are enabled in DpmjetIII307.
    _projectiles = dpmjet_extended_projectiles | Nuclei(a_max=280)
    _targets = Nuclei()
    _param_file_name = "dpmjpar.dat"
    _evap_file_name = "dpmjet.dat"
    _data_url = (
        "https://github.com/impy-project/chromo"
        "/releases/download/zipped_data_v1.0/dpm3191_v001.zip"
    )
    _ecm_min = 1 * GeV
    _max_A1 = 0
    _max_A2 = 0
    _max_plab = 0.0

    def __init__(self, evt_kin, *, seed=None):
        import chromo

        super().__init__(seed)

        data_dir = _cached_data_dir(self._data_url)
        # Set the dpmjpar.dat file
        if hasattr(self._lib, "pomdls") and hasattr(self._lib.pomdls, "parfn"):
            pfile = data_dir + self._param_file_name
            info(3, "DPMJET parameter file at", pfile)
            self._lib.pomdls.parfn = fortran_chars(self._lib.pomdls.parfn, pfile)

        # Set the data directory for the other files
        if hasattr(self._lib, "poinou") and hasattr(self._lib.poinou, "datdir"):
            pfile = data_dir
            info(3, "DPMJET data dir is at", pfile)
            self._lib.poinou.datdir = fortran_chars(self._lib.poinou.datdir, pfile)
            self._lib.poinou.lendir = len(pfile)
        # TODO: Rename the common block to chromo
        if hasattr(self._lib, "dtchro"):
            evap_file = data_dir + self._evap_file_name
            info(3, "DPMJET evap file at", evap_file)
            self._lib.dtchro.fnevap = fortran_chars(self._lib.dtchro.fnevap, evap_file)

        # Setup logging
        lun = 6  # stdout
        if hasattr(self._lib, "dtflka"):
            self._lib.dtflka.lout = lun
            self._lib.dtflka.lpri = 5 if chromo.debug_level else 1
        elif hasattr(self._lib, "dtiont"):
            self._lib.dtiont.lout = lun
        else:
            assert False, "Unknown DPMJET version, IO common block not detected"
        self._lib.pydat1.mstu[10] = lun

        self.kinematics = evt_kin

        # Relax momentum and energy conservation checks at very high energies
        if evt_kin.ecm > 5e4:
            # Relative allowed deviation
            self._lib.pomdls.parmdl[74] = 0.05
            # Absolute allowed deviation
            self._lib.pomdls.parmdl[75] = 0.05

        # Prevent DPMJET from overwriting decay settings
        self._lib.dtfrpa.ovwtdc = False
        # Tell PHOJET to not overwrite decay settings
        self._lib.pomdls.iswmdl[6 - 1] = 4
        # Recover the decay settings due to how DPMJET works
        self._lib.pydat1.mstj[21 - 1] = 1
        self._lib.pydat1.mstj[22 - 1] = 1

        self._set_final_state_particles()

    def _run_glauber(self, kin, photon_x, prod_only):
        """Run Glauber calculation for nuclear cross sections.

        Parameters
        ----------
        prod_only : bool
            If True, only compute production cross section (fast).
            If False, compute all components including total/elastic (slow).
        """
        self._lib.dtglgp.lprod = prod_only
        self._lib.dt_xsglau(
            kin.p1.A or 1,
            kin.p2.A or 1,
            (
                self._lib.idt_icihad(2212)
                if (kin.p1.A and kin.p1.A > 1)
                else self._lib.idt_icihad(kin.p1)
            ),
            photon_x,
            kin.virt_p1,
            kin.ecm,
            1,
            1,
            1,
        )

    def _cross_section(self, kin=None, photon_x=0, max_info=False):
        kin = self.kinematics if kin is None else kin
        # we override to set precision
        if (
            (kin.p1.is_nucleus and kin.p1.A > 1) or (kin.p2.is_nucleus and kin.p2.A > 1)
        ) and max_info:
            assert kin.p2.A >= 1, "DPMJET requires nucleons or nuclei on side 2."
            self._run_glauber(kin, photon_x, prod_only=False)
            glxs = self._lib.dtglxs

            def _generate():
                raise RuntimeError(
                    "Do not generate events with DPMJET after calculations "
                    "of nuclear cross sections."
                )

            self._generate = _generate
            return CrossSectionData(
                total=glxs.xstot[0, 0, 0],
                elastic=glxs.xsela[0, 0, 0],
                inelastic=glxs.xstot[0, 0, 0] - glxs.xsela[0, 0, 0],
                prod=glxs.xspro[0, 0, 0],
                quasielastic=glxs.xsqep[0, 0, 0]
                + glxs.xsqet[0, 0, 0]
                + glxs.xsqe2[0, 0, 0]
                + glxs.xsela[0, 0, 0],
            )
        if (kin.p1.is_nucleus and kin.p1.A > 1) or (kin.p2.is_nucleus and kin.p2.A > 1):
            # The value cached in dtglxs.xspro during initialisation is
            # valid only at the initialization kinematics, so it must not
            # be returned for arbitrary queries (issue #242). Run the
            # production-only Glauber MC for the requested kinematics,
            # saving and restoring the RNG state (all Fortran draws go
            # through the numpy bit generator) so that event generation
            # streams stay untouched.
            rng_state = self.random_state
            saved_lprod = self._lib.dtglgp.lprod
            try:
                self._run_glauber(kin, photon_x, prod_only=True)
                prod = self._lib.dtglxs.xspro[0, 0, 0]
            finally:
                self._lib.dtglgp.lprod = saved_lprod
                self.random_state = rng_state
            return CrossSectionData(
                prod=prod,
            )
        if kin.p1 == 22 and kin.p2.A == 1:
            stot, sine, _ = self._lib.dt_siggp(photon_x, kin.virt_p1, kin.ecm, 0)
            return CrossSectionData(total=stot, inelastic=sine, elastic=stot - sine)
        stot, sela = self._lib.dt_xshn(
            self._lib.idt_icihad(kin.p1), self._lib.idt_icihad(kin.p2), 0.0, kin.ecm
        )
        return CrossSectionData(total=stot, elastic=sela, inelastic=stot - sela)

    @property
    def glauber_trials(self):
        """Number of trials for Glauber model integration

        Default is 1000 (set at model initialisation).
        Larger number of `ntrials` reduces the fluctuations in the cross section,
        thus, making it more smooth. Smaller number of `ntrials` makes calculations of
        cross section faster.
        """
        return self._lib.dtglgp.jstatb

    @glauber_trials.setter
    def glauber_trials(self, ntrials):
        self._lib.dtglgp.jstatb = ntrials

    def _set_kinematics(self, kin):
        # Save maximal mass that has been initialized
        # (DPMJET sometimes crashes if higher mass requested than initialized)
        if not self._max_A1:
            # only do this once
            if kin.frame == EventFrame.FIXED_TARGET:
                self._lib.dtflg1.iframe = 1
                self._frame = EventFrame.FIXED_TARGET
            else:
                self._lib.dtflg1.iframe = 2
                self._frame = EventFrame.CENTER_OF_MASS
            self._max_A1 = kin.p1.A or 1
            self._max_A2 = kin.p2.A or 1
            self._max_plab = max(kin.plab, 100.0)
            self._lib.dt_init(
                -1,
                self._max_plab,
                kin.p1.A or 1,
                kin.p1.Z or 0,
                kin.p2.A or 1,
                kin.p2.Z or 0,
                kin.p1,
                iglau=0,
            )

        if (kin.p1.A or 1) > self._max_A1 or (kin.p2.A or 1) > self._max_A2:
            msg = (
                "Maximal initialization mass exceeded "
                f"{kin.p1.A}/{self._max_A1}, {kin.p2.A}/{self._max_A2}"
            )
            raise ValueError(msg)

        # PHOJET h-N tables end at the init energy; PHO_CSINT
        # log-extrapolates above it (warning only at LPRi > 4)
        if kin.plab > self._max_plab * (1.0 + 1e-9):
            msg = (
                f"plab = {kin.plab:.6g} GeV/c exceeds the initialization "
                f"momentum {self._max_plab:.6g} GeV/c; cross sections "
                "are tabulated only up to the initialization energy. "
                "Initialize at the highest energy of the simulation and "
                "set lower-energy kinematics afterwards via "
                "generator.kinematics = ..."
            )
            raise ValueError(msg)

        # AF: No idea yet, but apparently this functionality was around?!
        # if hasattr(k, 'beam') and hasattr(self._lib, 'init'):
        #     self._lib.dt_setbm(k.A1, k.Z1, k.A2, k.Z2, k.beam[0], k.beam[1])
        #     print 'OK'

    def _set_stable(self, pdgid, stable):
        kc = self._lib.pycomp(pdgid)
        self._lib.pydat3.mdcy[kc - 1, 0] = not stable

    def _generate(self):
        k = self.kinematics
        reject = self._lib.dt_kkinc(
            k.p1.A or 1,
            k.p1.Z or 0,
            k.p2.A or 1,
            k.p2.Z or 0,
            (
                self._lib.idt_icihad(2212)
                if (k.p1.A and k.p1.A > 1)
                else self._lib.idt_icihad(k.p1)
            ),
            k.elab,
            kkmat=-1,
        )
        self._lib.dtevno.nevent += 1
        return not reject

    def print_native_event(self, mode=1):
        if hasattr(self._lib, "dtflka"):
            saved_lpri = self._lib.dtflka.lpri
            self._lib.dtflka.lpri = 5
        self._lib.dt_evtout(mode)
        self._lib.dtflka.lpri = saved_lpri


class DpmjetIII193(DpmjetIIIRun):
    _version = "19.3"
    _library_name = "_dpmjet193"


class DpmjetIII307(DpmjetIIIRun):
    """DPMJET 3.0-7 with PHOJET 1.12 as the hadron/photon-nucleon engine.

    In addition to the standard hadronic projectiles and nuclei, DPMJET
    supports photons (PDG 22) as projectiles on nuclear targets
    (photon-induced interactions, VDM + Glauber). The Fortran code of
    DPMJET 3.0-7 selects the projectile identity (``IJPROJ``) from the
    ``/DTPRTA/`` common block instead of the ``IDP`` argument of
    ``DT_INIT``, and the photon cross-section tables in the Glauber
    module are only initialized if ``IJPROJ = 7`` (the BAMJET index of
    the photon) is already set when the Glauber initialization runs.
    Therefore, the first call to ``DT_INIT`` (triggered by setting the
    kinematics) is transparently redone here with ``IJPROJ = 7`` when a
    photon projectile is requested.

    Photons on nucleon targets (gamma + p / n) are not enabled, since
    then no Glauber initialization is performed and the production
    cross section that the interface caches when the kinematics is set
    would be meaningless. Note that PHOJET 1.12 alone
    (``Phojet112`` via :mod:`chromo.models.phojet`) supports photons
    on nucleon targets.
    """

    _version = "3.0-7"
    _library_name = "_dpmjet307"
    _projectiles = standard_projectiles | {PDGID(22)} | Nuclei(a_max=280)
    _param_file_name = "fitpar.dat"
    _data_url = (
        "https://github.com/impy-project/chromo"
        "/releases/download/zipped_data_v1.0/dpm3_v001.zip"
    )

    @classmethod
    def _pair_allowed(cls, p1, p2):
        if p1 == 22 and pdg2AZ(p2)[0] == 1:
            return False
        return True

    def _check_kinematics(self, kin):
        super()._check_kinematics(kin)
        if not self._pair_allowed(abs(kin.p1), abs(kin.p2)):
            msg = (
                "DpmjetIII307 supports photon projectiles only on nuclear "
                "targets (A > 1); use Phojet112 or Pythia8 for gamma + "
                "nucleon interactions."
            )
            raise ValueError(msg)

    def _set_kinematics(self, kin):
        super()._set_kinematics(kin)
        if abs(kin.p1) != 22:
            return
        # DPMJET 3.0-7 reads the projectile identity for the Glauber
        # initialization from the /DTPRTA/ common block (IJPROJ),
        # ignoring the IDP argument of DT_INIT, and the very first DT_INIT
        # resets IJPROJ to 1 (via DT_DEFAUL). Photon-nucleus cross
        # sections are only tabulated when IJPROJ = 7 (the BAMJET index
        # of the photon) while DT_SHMAKI/DT_XSGLAU run, so repeat the
        # initialization with the photon index set. EPN is kept at the
        # maximal lab momentum of the run (self._max_plab) so that the
        # hadronic path is not degraded, and /DTPRTA/ is restored
        # afterwards. The repeated initialization is skipped as long as
        # the target nucleus does not change.
        target_key = (kin.p2.A or 1, kin.p2.Z or 0)
        if getattr(self, "_photon_init_target", None) == target_key:
            return
        self._lib.dtprta.ijproj = 7
        self._lib.dtprta.ibproj = 7
        try:
            self._lib.dt_init(
                -1,
                self._max_plab,
                1,
                0,
                *target_key,
                22,
                iglau=0,
            )
        finally:
            self._lib.dtprta.ijproj = 1
            self._lib.dtprta.ibproj = 1
        self._photon_init_target = target_key

    def _run_glauber(self, kin, photon_x, prod_only):
        if abs(kin.p1) == 22:
            # Use target slot 2 (NIDX = 2) for photon-induced cross sections
            # so that the hadronic results in slot 1, which are tabulated
            # once during DT_INIT and reused by the base class, remain intact.
            self._lib.dtglgp.lprod = prod_only
            self._lib.dt_xsglau(
                1,  # photon has no nucleons
                kin.p2.A or 1,
                7,  # BAMJET index of the photon projectile
                photon_x,
                kin.virt_p1,
                kin.ecm,
                1,
                1,
                2,
            )
            return
        super()._run_glauber(kin, photon_x, prod_only)

    def _cross_section(self, kin=None, photon_x=0, max_info=False):
        kin = self.kinematics if kin is None else kin
        if abs(kin.p1) == 22 and kin.p2.A and kin.p2.A > 1:
            # Photon-nucleus cross sections are computed with the Glauber
            # module (DT_XSGLAU with IJPROJ=7, VDM); the DTGLXS arrays are
            # populated by the call below into target slot 2.
            self._run_glauber(kin, photon_x, prod_only=not max_info)
            if max_info:
                # mirror the base class: the Glauber MC consumed Fortran
                # RNG draws, so subsequent events are not reproducible
                def _generate():
                    raise RuntimeError(
                        "Do not generate events with DPMJET after "
                        "calculations of nuclear cross sections."
                    )

                self._generate = _generate
            glxs = self._lib.dtglxs
            stot = glxs.xstot[0, 0, 1]
            sela = glxs.xsela[0, 0, 1]
            return CrossSectionData(
                total=stot,
                elastic=sela,
                inelastic=stot - sela,
                prod=glxs.xspro[0, 0, 1],
                quasielastic=glxs.xsqep[0, 0, 1]
                + glxs.xsqet[0, 0, 1]
                + glxs.xsqe2[0, 0, 1]
                + sela,
            )
        return super()._cross_section(kin, photon_x=photon_x, max_info=max_info)


class DpmjetIII193_DEV(DpmjetIIIRun):
    _version = "19.3-dev"
    _library_name = "_dev_dpmjet193"
