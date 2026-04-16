"""FLUKA event-generator interface.

Wraps FLUKA 2025.1 via f2py-exported Fortran helpers linked against the
official FLUKA archives (libflukahp.a, libdpmmvax.a, librqmdmvax.a,
librqmd.a, libdpmjet*.a). Supports hN, hA, AA, photohadronic,
photonuclear, and EMD interactions, with nuclear remnants in HEPEVT.
"""

import logging
import os
import pathlib
import tempfile
from enum import IntEnum

import numpy as np
from particle import Particle
from particle import literals as lp

from chromo.common import CrossSectionData, MCEvent, MCRun
from chromo.constants import GeV, TeV, nucleon_mass, standard_projectiles
from chromo.kinematics import EventFrame
from chromo.util import (
    CompositeTarget,
    Nuclei,
    info,
    is_real_nucleus,
    pdg2AZ,
    process_particle,
)

log = logging.getLogger(__name__)


class InteractionType(IntEnum):
    """FLUKA IFLXYZ flag decoding.

    Digit positions (decimal): units=inelastic, tens=elastic, hundreds=EMD.
    """

    INELASTIC = 1
    ELASTIC = 10
    INELA_ELA = 11
    EMD = 100
    INELA_EMD = 101
    ELA_EMD = 110
    INELA_ELA_EMD = 111


# Default targets registered at construction. Covers cosmic-ray air,
# common heavy-ion experiments, and typical laboratory materials.
# Entries are PDG ids (ints or strings resolved via process_particle).
#
# IMPORTANT: FLUKA's STPXYZ allocates geometry with GLBCRD(WHAT(1)=10)
# which limits MEDFLK to 10 regions. Total unique materials (defaults +
# user + kinematics target) must not exceed 10. Keep this list ≤ 9 to
# leave room for one extra target not in the list.
_DEFAULT_MATERIALS = (
    2212,  # free proton
    "He4",
    "C12",
    "N14",
    "O16",
    "Ar40",
    "Fe56",
    "Cu63",
    "Pb208",
)


class FlukaEvent(MCEvent):
    """Event data from FLUKA's HEPEVT common block."""

    def _get_charge(self, npart):
        pids = self._lib.hepevt.idhep[:npart]
        # charge_from_pdg_arr (Fortran) returns 0 for nuclei because MCIHAD
        # returns 0 for |pid|>=1e9. Recover the nucleus charge Z directly
        # from the PDG ion code 10LZZZAAAI via vectorised numpy ops.
        result = self._lib.charge_from_pdg_arr(pids).astype(float)
        absp = np.abs(pids.astype(np.int64))
        nucleus = absp >= 1_000_000_000
        z = ((absp // 10000) % 1000).astype(float)
        result = np.where(nucleus, np.sign(pids) * z, result)
        return result

    def _history_zero_indexing(self):
        pass

    def _prepend_initial_beam(self):
        pass

    def _repair_initial_beam(self):
        pass


class Fluka(MCRun):
    """FLUKA event generator (2025.1).

    Supports hN, hA, AA, photohadronic (gamma+p), photonuclear (gamma+A), and
    electromagnetic dissociation (EMD). Generator is single-instantiation
    per Python process.

    Generator transitions
    ---------------------
    FLUKA dispatches to different hadronic models depending on the
    lab-frame kinetic energy per nucleon:

    * **hA** (hadron-nucleus): Peanut below ``transition_peanut_dpmjet``,
      DPMJET-3 above. ``transition_peanut_dpmjet_smearing`` randomises
      the switch-over to avoid discontinuities (set to 0 for a sharp
      handoff).
    * **AA** (nucleus-nucleus): rQMD below ``transition_rqmd_dpmjet``,
      DPMJET-3 above. Below roughly 125 MeV/n (``QMDIOT``, not user-
      configurable) FLUKA's BME model takes over from rQMD.
    * **DPMJET → UHE** (EPOS/SIBYLL/QGSJET) has no effect in chromo's
      FLUKA build: FLUKA's own runtime UHE thresholds sit at 1e+30, and
      no UHE model is linked in any case. The top of chromo's supported
      range is the class-level sanity ceiling ``_max_sqrt_s_nn``
      (500 TeV by default); kinematics above this raise ``ValueError``
      at construction.

    The class-attribute defaults mirror FLUKA's post-STPXYZ runtime
    defaults and can be overridden instance-by-instance via the keyword
    arguments below.

    Parameters
    ----------
    evt_kin : EventKinematicsBase
        Initial kinematics. Target must be registered (see ``targets``).
    seed : int or None
        Random seed for FLUKA's Ranmar generator.
    targets : iterable of (str|int|PDGID), optional
        Extra nuclei to register beyond ``_DEFAULT_MATERIALS``. Required
        if you plan to shoot at a nucleus outside the default list.
    interaction_type : InteractionType
        Which channels to generate/compute. Default INELASTIC.
    transition_peanut_dpmjet : float or None
        hA Peanut↔DPMJET-3 transition energy (GeV, lab kin. energy).
        ``None`` → class default (``_transition_peanut_dpmjet``).
    transition_peanut_dpmjet_smearing : float or None
        Smearing (±, GeV) of the hA Peanut↔DPMJET transition. Set to 0
        for a sharp switch. ``None`` → class default.
    transition_rqmd_dpmjet : float or None
        AA rQMD↔DPMJET-3 transition energy (GeV/n, lab kin. energy per
        nucleon). ``None`` → class default.
    transition_rqmd_dpmjet_smearing : float or None
        Smearing (±, GeV/n) of the AA rQMD↔DPMJET transition.
        ``None`` → class default.
    enable_quasielastic : bool
        Enable quasi-elastic scattering. Default False.
    rng_state_file : pathlib.Path or str, optional
        File used to persist Ranmar state across runs.
    """

    _name = "FLUKA"
    _event_class = FlukaEvent
    _frame = EventFrame.FIXED_TARGET
    _version = "2025.1"
    _library_name = "_fluka"
    _projectiles = standard_projectiles | Nuclei() | {lp.photon.pdgid}
    _targets = Nuclei() | {2212, 2112}

    # ------------------------------------------------------------------
    # Hadronic-generator transition defaults (from flukapro/(GENTHR),
    # cross-checked against FLUKA's post-STPXYZ runtime state — the
    # static PARAMETERs in (GENTHR) are superseded by BDNOPT/DFLTS init).
    # Values are lab kinetic energy (GeV for hA, GeV/n for AA).
    # ------------------------------------------------------------------

    # hA: Peanut handles below this, DPMJET-3 above. FLUKA's runtime
    # default is 20 TeV (scalar DPJHDT); Peanut is accurate well beyond
    # 1 TeV for most hadrons. Smearing disabled by default (LFDSMR=0).
    _transition_peanut_dpmjet = 20 * TeV
    _transition_peanut_dpmjet_smearing = 0.0

    # AA: rQMD handles the intermediate range (from QMDIOT ≈ 125 MeV/n
    # upwards) and DPMJET-3 takes over above this. FLUKA runtime
    # default: 12.5 GeV/n with 2.0 GeV/n smearing (LAASMR=1).
    _transition_rqmd_dpmjet = 12.5 * GeV
    _transition_rqmd_dpmjet_smearing = 2.0 * GeV

    # ------------------------------------------------------------------
    # Kinematics ceiling.
    # ------------------------------------------------------------------
    # FLUKA's runtime UHE boundaries (UHEHDT/UHEIOT) sit at 1e+30 out of
    # the box, so DPMJET naturally handles the entire high-energy range;
    # no UHE model is linked in chromo's FLUKA build in any case. This
    # ceiling is chromo's own sanity limit: above it we raise before
    # DPMJET has a chance to misbehave.
    _max_sqrt_s_nn = 500 * TeV
    # Convert: for ultra-relativistic kinematics, s ≈ 2·m_N·E_kin_lab.
    _ekin_per_nucleon_max = _max_sqrt_s_nn**2 / (2 * nucleon_mass)
    _ecm_min = 0.0  # liberal; FLUKA itself has no practical low-E floor

    def __init__(
        self,
        evt_kin,
        *,
        seed=None,
        targets=None,
        interaction_type=InteractionType.INELASTIC,
        transition_peanut_dpmjet=None,
        transition_peanut_dpmjet_smearing=None,
        transition_rqmd_dpmjet=None,
        transition_rqmd_dpmjet_smearing=None,
        enable_quasielastic=False,
        rng_state_file=None,
    ):
        super().__init__(seed)

        if (
            "FLUPRO" not in os.environ
            or not pathlib.Path(os.environ["FLUPRO"]).exists()
        ):
            raise RuntimeError(
                "FLUPRO environment variable is not set or points to a "
                "non-existing directory — run scripts/install_fluka.sh"
            )

        # Resolve transition thresholds (ctor override → class default).
        if transition_peanut_dpmjet is not None:
            self._transition_peanut_dpmjet = float(transition_peanut_dpmjet)
        if transition_peanut_dpmjet_smearing is not None:
            self._transition_peanut_dpmjet_smearing = float(
                transition_peanut_dpmjet_smearing
            )
        if transition_rqmd_dpmjet is not None:
            self._transition_rqmd_dpmjet = float(transition_rqmd_dpmjet)
        if transition_rqmd_dpmjet_smearing is not None:
            self._transition_rqmd_dpmjet_smearing = float(
                transition_rqmd_dpmjet_smearing
            )

        self._interaction_type = int(interaction_type)
        self._init_rng(rng_state_file, seed)
        self._set_quasielastic(enable_quasielastic)
        self._init_fluka_materials(evt_kin, targets or ())
        # GENTHR defaults are installed by FLUKA during STPXYZ; override
        # them AFTER the call so our values are final.
        self._set_generator_thresholds()

        self.kinematics = evt_kin
        self._set_final_state_particles()
        self._activate_decay_handler(on=True)

    # ------------------------------------------------------------------
    # Internal helpers — RNG
    # ------------------------------------------------------------------

    def _init_rng(self, rng_state_file, seed):
        if rng_state_file is None:
            rng_state_file = pathlib.Path(tempfile.gettempdir()) / "fluka_rng_state.dat"
        self._rng_state_file = pathlib.Path(rng_state_file)
        self._logical_unit = 888

        if self._rng_state_file.exists() and self._rng_state_file.stat().st_size > 0:
            self._lib.load_rng_state(str(self._rng_state_file), self._logical_unit)
        else:
            seed_int = 0 if seed is None else int(seed)
            self._lib.init_rng_state(
                str(self._rng_state_file),
                self._logical_unit,
                seed_int,
                0,
                0,
            )
            self._lib.save_rng_state(str(self._rng_state_file), self._logical_unit)

    def _set_quasielastic(self, enable):
        # FLUKA common-block toggles; names as exposed by f2py.
        self._lib.qelcmm.lxsqel = 0
        self._lib.qelcmm.lpqels = 1 if enable else 0
        self._lib.nucflg.lqecmp = 1 if enable else 0

    # ------------------------------------------------------------------
    # Internal helpers — materials
    # ------------------------------------------------------------------

    def _build_material_list(self, evt_kin, user_targets):
        """Return ordered list of PDG ids (ints) covering defaults + user + kin."""
        seen = []
        for src in (_DEFAULT_MATERIALS, tuple(user_targets)):
            for item in src:
                pdg = int(process_particle(item))
                if pdg not in seen:
                    seen.append(pdg)
        target = process_particle(evt_kin.p2)
        if isinstance(target, CompositeTarget):
            for comp in target.components:
                pdg = int(comp)
                if pdg not in seen:
                    seen.append(pdg)
        else:
            pdg = int(target)
            if pdg not in seen:
                seen.append(pdg)
        return seen

    def _init_fluka_materials(self, evt_kin, user_targets):
        materials = self._build_material_list(evt_kin, user_targets)
        # FLUKA STPXYZ is initialised with GLBCRD(WHAT(1)=10), limiting the
        # allocatable MEDFLK geometry array to 10 regions.
        if len(materials) > 10:
            msg = (
                f"Too many materials ({len(materials)} > 10). "
                "FLUKA STPXYZ is limited to 10 geometry regions. "
                "Reduce the default list or pass fewer targets= entries."
            )
            raise ValueError(msg)
        # Each entry is a single-element material; nelmfl[i]=1 for all.
        nelmfl = np.ones(len(materials), dtype=np.int32)
        # pdg2AZ returns (A, Z); we want Z. Free proton (2212) → (1, 1).
        izelfl = np.array([max(pdg2AZ(pdg)[1], 1) for pdg in materials], dtype=np.int32)
        wfelfl = np.ones(len(materials), dtype=np.float64)
        lprint = 0  # suppress FLUKA material printout

        # STPXYZ's ef2dp3/df2dp3 only cover the hA Peanut↔DPMJET
        # transition; pass sentinels and override via GENTHR afterwards
        # (see _set_generator_thresholds) so all transitions go through
        # a single mechanism.
        #
        # pptmax is the maximum momentum used to build DPMJET's
        # initialisation tables. Keep it at 1 EeV/c — empirically the
        # largest value FLUKA's DPMJET init accepts without regressing
        # low-energy event generation. The UHE dispatch ceiling (set in
        # _set_generator_thresholds) is separate from this.
        mt = self._lib.chromo_stpxyz(
            nelmfl,
            izelfl,
            wfelfl,
            1e9,  # pptmax (GeV/c)
            -1.0,  # ef2dp3: FLUKA default, overridden below
            -1.0,  # df2dp3: FLUKA default, overridden below
            int(self._interaction_type),
            lprint,
        )
        self._materials_pdg = tuple(materials)
        self._materials_idx = np.asarray(mt, dtype=np.int32)
        self._materials_map = dict(zip(materials, self._materials_idx.tolist()))

    def _set_generator_thresholds(self):
        """Write hadronic-generator transition thresholds into GENTHR.

        Must be called AFTER ``chromo_stpxyz`` so FLUKA's own
        initialisation can't overwrite our values. See
        ``flukapro/(GENTHR)`` for variable semantics.

        Touches only the interaction thresholds and their smearing
        settings — leaves cross-section thresholds (``DPHDXT``/
        ``DPIOXT``) and the UHE/rQMD lower bounds (``UHEHDT``/
        ``UHEIOT``/``QMDIOT``) at their FLUKA runtime defaults.
        """
        g = self._lib.genthr

        # hA: Peanut ↔ DPMJET-3 (scalar DPJHDT; particle-specific Peanut
        # thresholds PEANCT/PEAPIT/... are not touched).
        g.dpjhdt = self._transition_peanut_dpmjet
        g.fldpsm = self._transition_peanut_dpmjet_smearing
        g.lfdsmr = 1 if self._transition_peanut_dpmjet_smearing > 0 else 0

        # AA: rQMD ↔ DPMJET-3. QMDIOT is rQMD's *lower* bound
        # (≈125 MeV/n), not the rQMD↔DPMJET switch — leave it alone.
        g.dpjiot = self._transition_rqmd_dpmjet
        g.dpqmsm = self._transition_rqmd_dpmjet_smearing
        g.laasmr = 1 if self._transition_rqmd_dpmjet_smearing > 0 else 0

    # ------------------------------------------------------------------
    # Internal helpers — projectile & target codes
    # ------------------------------------------------------------------

    def _fluka_projectile_code(self, pdg_id):
        return int(self._lib.pdg_to_proj_code(int(pdg_id)))

    def _get_material_index(self, particle):
        p = process_particle(particle)
        if isinstance(p, CompositeTarget):
            # CompositeTarget should be handled by base class iteration.
            raise KeyError(
                "CompositeTarget should be handled by _temporary_kinematics; "
                "_get_material_index received composite."
            )
        try:
            return int(self._materials_map[int(p)])
        except KeyError:
            name = Particle.from_pdgid(int(p)).name
            msg = (
                f"target {name} (pdg={int(p)}) not initialised; "
                f"pass targets=['{name}'] to the Fluka constructor"
            )
            raise KeyError(msg)

    # ------------------------------------------------------------------
    # Abstract overrides
    # ------------------------------------------------------------------

    def _check_kinematics(self, kin):
        super()._check_kinematics(kin)
        # Single ceiling: FLUKA has no practical low-energy floor (its
        # transport threshold is 1 µeV/n, well below anything chromo
        # users would set), and above _ekin_per_nucleon_max DPMJET would
        # hand off to an unlinked UHE model and abort.
        a = kin.p1.A or 1
        ekin_per_n = kin.ekin / a
        if ekin_per_n > self._ekin_per_nucleon_max:
            msg = (
                f"kinetic energy/nucleon {ekin_per_n / GeV:.3g} GeV > "
                f"FLUKA ceiling {self._ekin_per_nucleon_max / GeV:.3g} GeV "
                f"(sqrt(s_NN) = {self._max_sqrt_s_nn / TeV:.0f} TeV). "
                f"DPMJET->UHE transition not supported in chromo's FLUKA "
                f"build - no UHE model is linked."
            )
            raise ValueError(msg)

    def _set_kinematics(self, kin):
        # Resolve the target material index (required before _generate).
        # For CompositeTarget, use the first component as the default;
        # _composite_plan will call _set_kinematics again with each component
        # before actually generating events.
        p2 = kin.p2
        if isinstance(p2, CompositeTarget):
            first_component = p2.components[0]
            self._current_target_idx = self._get_material_index(first_component)
        else:
            self._current_target_idx = self._get_material_index(p2)

    def _cross_section(self, kin=None, max_info=False):
        kin = self.kinematics if kin is None else kin
        proj_code = self._fluka_projectile_code(int(kin.p1))
        mat_idx = (
            self._current_target_idx
            if not isinstance(process_particle(kin.p2), CompositeTarget)
            else self._get_material_index(kin.p2.components[0])
        )
        ekin = kin.ekin
        # Three channels at most; route by interaction_type flag.
        flag = int(self._interaction_type)
        inel = el = emd = np.nan
        if flag % 10 == 1:
            inel = float(self._lib.chromo_sgmxyz(proj_code, mat_idx, ekin, 0.0, 1))
        if (flag // 10) % 10 == 1 and int(kin.p1) != 22:
            el = float(self._lib.chromo_sgmxyz(proj_code, mat_idx, ekin, 0.0, 10))
        if (flag // 100) % 10 == 1:
            emd = float(self._lib.chromo_sgmxyz(proj_code, mat_idx, ekin, 0.0, 100))
        # For nuclear targets, FLUKA's inelastic (with quasielastic disabled,
        # the chromo default) is the production cross section. Populate both
        # `inelastic` and `prod` so downstream code that keys on either works.
        prod = inel if is_real_nucleus(kin.p2) else np.nan
        return CrossSectionData(inelastic=inel, elastic=el, emd=emd, prod=prod)

    def _set_stable(self, pdgid, stable):
        info(
            2,
            f"Fluka._set_stable is a no-op (pdgid={pdgid}, stable={stable}). "
            "FLUKA decay settings are global.",
        )

    def _generate(self):
        k = self.kinematics
        # EVTXYZ takes a different projectile-code encoding than SGMXYZ:
        # heavy ions (A>4) must be registered via PDGION first and then
        # passed as the -2 HEAVYION sentinel. pdg_to_evt_code handles
        # both steps; for hadrons and light nuclei it is equivalent to
        # the SGMXYZ encoding produced by pdg_to_proj_code.
        proj_code = int(self._lib.pdg_to_evt_code(int(k.p1)))
        if proj_code == -2:
            msg = (
                "Heavy-ion projectiles (A > 4) are not supported for event "
                "generation in FLUKA — EVTXYZ aborts for nuclear projectiles. "
                "Use hadronic projectiles (p, pi, K, n, gamma) for events; "
                "cross_section() works for AA kinematics."
            )
            raise ValueError(msg)
        mat_idx = self._current_target_idx
        ekin = k.ekin
        self._lib.chromo_evtxyz(
            proj_code,
            mat_idx,
            ekin,
            0.0,  # pproj0; ekin takes precedence
            0.0,
            0.0,
            1.0,  # on-axis +z
            int(self._interaction_type),
        )
        self._lib.chromo_fllhep()
        # If chromo_fill_remnants was added (see Task 9), call it too.
        if hasattr(self._lib, "chromo_fill_remnants"):
            self._lib.chromo_fill_remnants()
        return True

    # ------------------------------------------------------------------
    # Cleanup
    # ------------------------------------------------------------------

    def _cleanup_fort(self):
        for fort_file in pathlib.Path(".").glob("fort.*"):
            fort_file.unlink(missing_ok=True)
        for f in (
            pathlib.Path(tempfile.gettempdir()) / "fluka_rng_state.dat",
            pathlib.Path(".") / ".timer.out",
        ):
            f.unlink(missing_ok=True)

    def __del__(self):
        try:
            self._cleanup_fort()
        except Exception:
            pass

    # ------------------------------------------------------------------
    # Public accessors
    # ------------------------------------------------------------------

    def save_rng_state(self, file=None):
        target = pathlib.Path(file) if file else self._rng_state_file
        self._lib.save_rng_state(str(target), self._logical_unit)

    def load_rng_state(self, file=None):
        target = pathlib.Path(file) if file else self._rng_state_file
        self._lib.load_rng_state(str(target), self._logical_unit)

    def fluka_rand(self):
        return float(self._lib.fluka_rand())

    @property
    def registered_targets(self):
        """Tuple of PDG ids registered in FLUKA's material tables."""
        return self._materials_pdg

    def register_target(self, pdg):
        """Register an additional target nucleus at runtime.

        Raises
        ------
        NotImplementedError
            FLUKA 2025.1 does not expose a safe post-STPXYZ material-
            registration API. STPXYZ itself aborts (FLABRT) on a second
            call. While FLKMAT common block variables (nmat, ztar, amss)
            are readable via F2PY, the downstream initialisation routines
            SETITB (registers material with PEANUT nuclear tables) and
            DFATWG (computes atomic weight) are not wrapped and therefore
            not callable from Python. Direct manipulation of FLKMAT
            without refreshing these tables would produce silent physics
            errors or segfaults.

        Note
        ----
        Pass ``targets=[...]`` to the ``Fluka`` constructor instead to
        register all required nuclei before the STPXYZ call.
        """
        raise NotImplementedError(
            "Runtime target extension not supported by FLUKA 2025.1. "
            "STPXYZ aborts on a second call, and SETITB/DFATWG are not "
            "wrapped. Pass targets=[...] to the Fluka constructor."
        )
