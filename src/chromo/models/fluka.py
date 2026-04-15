"""FLUKA event-generator interface.

Wraps FLUKA 2025.1 via f2py-exported Fortran helpers linked against the
official FLUKA archives (libflukahp.a, libdpmmvax.a, librqmdmvax.a,
librqmd.a, libdpmjet*.a). Supports hN, hA, AA, photohadronic,
photonuclear, and EMD interactions, with nuclear remnants in HEPEVT.
"""

import logging
import os
import pathlib
from enum import IntEnum

import numpy as np
from particle import Particle
from particle import literals as lp

from chromo.common import CrossSectionData, MCEvent, MCRun
from chromo.constants import GeV, TeV, standard_projectiles
from chromo.kinematics import EventFrame
from chromo.util import CompositeTarget, Nuclei, info, process_particle

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
    2212,       # free proton
    "He4",
    "C12",
    "N14",
    "O16",
    "Ar40",
    "Fe56",
    "Cu63",
    "Pb208",
)


def _pdg_to_zA(pdg_id):
    """Return (Z, A) for a PDG nucleus id or standard hadron.

    For the free proton (pdg=2212) returns (1, 1).
    For a PDG nucleus id 10LZZZAAAI returns (Z, A).
    """
    p = Particle.from_pdgid(pdg_id)
    return p.pdgid.Z or 0, p.pdgid.A or 0


class FlukaEvent(MCEvent):
    """Event data from FLUKA's HEPEVT common block."""

    def _get_charge(self, npart):
        pids = self._lib.hepevt.idhep[:npart]
        hadron_charges = self._lib.charge_from_pdg_arr(pids)
        # charge_from_pdg_arr returns 0 for nuclei (MCIHAD returns 0 for
        # |pid|>=1e9). Recompute charge for nuclei directly from the PDG
        # id: Z = (|pid| // 10000) % 1000.
        result = hadron_charges.astype(float)
        for i, pid in enumerate(pids):
            if abs(int(pid)) >= 1000000000:
                result[i] = (abs(int(pid)) // 10000) % 1000
                if pid < 0:
                    result[i] = -result[i]
        return result

    def _history_zero_indexing(self):
        pass

    def _prepend_initial_beam(self):
        pass

    def _repair_initial_beam(self):
        pass


class Fluka(MCRun):
    """FLUKA event generator (2025.1).

    Supports hN, hA, AA, photohadronic (γ+p), photonuclear (γ+A), and
    electromagnetic dissociation (EMD). Generator is single-instantiation
    per Python process.

    Parameters
    ----------
    evt_kin : EventKinematicsBase
        Initial kinematics. Target must be registered (see `targets`).
    seed : int or None
        Random seed for FLUKA's Ranmar generator.
    targets : iterable of (str|int|PDGID), optional
        Extra nuclei to register beyond `_DEFAULT_MATERIALS`. Required
        if you plan to shoot at a nucleus outside the default list.
    interaction_type : InteractionType
        Which channels to generate/compute. Default INELASTIC.
    transition_energy : float or None
        FLUKA→DPMJET transition energy (GeV). None → FLUKA defaults.
    transition_smearing : float or None
        Smearing (+/-) of the transition energy. None → FLUKA default.
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
    _targets = Nuclei()
    _ecm_min = 0.01 * GeV  # covers GDR region for γA
    _ekin_per_nucleon_max_hadron = 20 * TeV
    _ekin_per_nucleon_max_photon = 1 * TeV  # refined by investigation

    def __init__(
        self,
        evt_kin,
        *,
        seed=None,
        targets=None,
        interaction_type=InteractionType.INELASTIC,
        transition_energy=None,
        transition_smearing=None,
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

        self._interaction_type = int(interaction_type)
        self._init_rng(rng_state_file, seed)
        self._set_quasielastic(enable_quasielastic)
        self._init_fluka_materials(
            evt_kin,
            targets or (),
            pptmax=1e9 if transition_energy is None else max(transition_energy, 1e9),
            ef2dp3=-1.0 if transition_energy is None else float(transition_energy),
            df2dp3=-1.0 if transition_smearing is None else float(transition_smearing),
        )

        self.kinematics = evt_kin
        self._set_final_state_particles()
        self._activate_decay_handler(on=True)

    # ------------------------------------------------------------------
    # Internal helpers — RNG
    # ------------------------------------------------------------------

    def _init_rng(self, rng_state_file, seed):
        if rng_state_file is None:
            rng_state_file = pathlib.Path(__file__).parent / "fluka_rng_state.dat"
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

    def _init_fluka_materials(
        self,
        evt_kin,
        user_targets,
        pptmax,
        ef2dp3,
        df2dp3,
    ):
        materials = self._build_material_list(evt_kin, user_targets)
        # FLUKA STPXYZ is initialised with GLBCRD(WHAT(1)=10), limiting the
        # allocatable MEDFLK geometry array to 10 regions.
        if len(materials) > 10:
            raise ValueError(
                f"Too many materials ({len(materials)} > 10). "
                "FLUKA STPXYZ is limited to 10 geometry regions. "
                "Reduce the default list or pass fewer targets= entries."
            )
        # Each entry is a single-element material; nelmfl[i]=1 for all.
        nelmfl = np.ones(len(materials), dtype=np.int32)
        izelfl = np.array(
            [max(_pdg_to_zA(pdg)[0], 1) for pdg in materials], dtype=np.int32
        )
        wfelfl = np.ones(len(materials), dtype=np.float64)
        lprint = 0  # suppress FLUKA material printout

        mt = self._lib.chromo_stpxyz(
            nelmfl,
            izelfl,
            wfelfl,
            pptmax,
            ef2dp3,
            df2dp3,
            int(self._interaction_type),
            lprint,
        )
        self._materials_pdg = tuple(materials)
        self._materials_idx = np.asarray(mt, dtype=np.int32)
        self._materials_map = dict(zip(materials, self._materials_idx.tolist()))

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
        a = kin.p1.A or 1
        ekin_per_n = kin.ekin / a
        if int(kin.p1) == 22:
            upper = self._ekin_per_nucleon_max_photon
        else:
            upper = self._ekin_per_nucleon_max_hadron
        if ekin_per_n > upper:
            raise ValueError(
                f"kinetic energy/nucleon {ekin_per_n / GeV:.3g} GeV > "
                f"max {upper / GeV:.3g} GeV for this projectile class"
            )

    def _set_kinematics(self, kin):
        # Resolve the target material index (required before _generate).
        self._current_target_idx = self._get_material_index(kin.p2)

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
            inel = float(
                self._lib.chromo_sgmxyz(proj_code, mat_idx, ekin, 0.0, 1)
            )
        if (flag // 10) % 10 == 1 and int(kin.p1) != 22:
            el = float(
                self._lib.chromo_sgmxyz(proj_code, mat_idx, ekin, 0.0, 10)
            )
        if (flag // 100) % 10 == 1:
            emd = float(
                self._lib.chromo_sgmxyz(proj_code, mat_idx, ekin, 0.0, 100)
            )
        return CrossSectionData(inelastic=inel, elastic=el, emd=emd)

    def _set_stable(self, pdgid, stable):
        info(
            2,
            f"Fluka._set_stable is a no-op (pdgid={pdgid}, stable={stable}). "
            "FLUKA decay settings are global.",
        )

    def _generate(self):
        k = self.kinematics
        proj_code = self._fluka_projectile_code(int(k.p1))
        mat_idx = self._current_target_idx
        ekin = k.ekin
        self._lib.chromo_evtxyz(
            proj_code,
            mat_idx,
            ekin,
            0.0,   # pproj0; ekin takes precedence
            0.0,
            0.0,
            1.0,   # on-axis +z
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
            try:
                fort_file.unlink()
            except OSError:
                pass
        for f in (
            pathlib.Path(__file__).parent / "fluka_rng_state.dat",
            pathlib.Path(".") / ".timer.out",
        ):
            try:
                f.unlink()
            except OSError:
                pass

    def __del__(self):
        try:
            self._cleanup_fort()
        except Exception:  # noqa: BLE001
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
