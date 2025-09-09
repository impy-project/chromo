"""
FLUKA event generator interface.

This module contains the implementation of the FLUKA event generator interface,
including FlukaEvent and FlukaRun classes. FLUKA uses PEANUT at low energies
for hh and hA interactions, and DPMJET-III at higher energies and for AA
interactions respectively.
"""

import os
import pathlib
from enum import Enum

import numpy as np
from particle import literals as lp

from chromo.common import CrossSectionData, MCEvent, MCRun
from chromo.constants import GeV, TeV, standard_projectiles
from chromo.kinematics import EventFrame
from chromo.util import CompositeTarget, Nuclei, info, process_particle

# =============================================================================
# Constants and Enums
# =============================================================================


class InteractionType(Enum):
    """FLUKA interaction types for stpxyz, smgxyz, evtxyx."""

    INELASTIC = 1
    ELASTIC = 10
    INELA_ELA = 11
    EMD = 100
    ENELA_EMD = 101
    ELA_EMD = 110
    INELA_ELA_EMD = 111


# Default material components for FLUKA
_DEFAULT_MATERIALS = [
    2212,  # proton
    "N14",
    "O16",
    "Ar40",
    "Fe56",  # common nuclei
]


# =============================================================================
# Event and Run Classes
# =============================================================================


class FlukaEvent(MCEvent):
    """Wrapper class around FLUKA HEPEVT-style particle stack."""

    def _get_charge(self, npart):
        return self._lib.charge_from_pdg_arr(self._lib.hepevt.idhep[:npart])

    def _history_zero_indexing(self):
        pass

    def _prepend_initial_beam(self):
        pass

    def _repair_initial_beam(self):
        pass


class Fluka(MCRun):
    """FLUKA event generator implementation."""

    _name = "FLUKA"
    _event_class = FlukaEvent
    _frame = EventFrame.FIXED_TARGET
    _projectiles = (
        standard_projectiles
        | Nuclei()
        | {lp.photon.pdgid, lp.e_plus.pdgid, lp.e_minus.pdgid}
    )
    _targets = Nuclei()
    _ecm_min = 0.1 * GeV
    _ekin_max = 20 * TeV
    _version = "2025.1"
    _library_name = "_fluka"

    def __init__(
        self,
        evt_kin,
        *,
        seed=None,
        interaction_type=InteractionType.INELASTIC,
        max_momentum_p=1e11,
        fluka_dpmjet_transition=-1.0,
        transition_smearing=-1.0,
        material_print=True,
        enable_quasielastic=False,
        rng_state_file=None,
    ):
        super().__init__(seed)

        # Validate FLUKA environment
        if (
            "FLUPRO" not in os.environ
            or not pathlib.Path(os.environ["FLUPRO"]).exists()
        ):
            raise RuntimeError(
                "FLUPRO environment variable is not set or points to a "
                "non-existing directory"
            )

        # Store configuration
        self._interaction_type = interaction_type.value
        self._max_momentum_p = max_momentum_p
        self._fluka_dpmjet_transition = fluka_dpmjet_transition
        self._transition_smearing = transition_smearing
        self._material_print = material_print

        # Setup RNG
        self._init_rng(rng_state_file, seed)

        # Configure quasielastic interactions
        self._set_quasielastic(enable_quasielastic)

        if evt_kin.ekin >= self._ekin_max and not (self._fluka_dpmjet_transition > 0.0):
            msg = (
                f"The maximal energy kinetic energy {evt_kin.ekin} GeV exceeds "
                "FLUKA's maximal range without DPMJET"
            )
            raise RuntimeError(msg)

        # Initialize materials and FLUKA
        self._init_fluka_materials(evt_kin)

        self.kinematics = evt_kin
        self._set_final_state_particles()
        self._activate_decay_handler(on=True)

    def _init_rng(self, rng_state_file, seed):
        """Initialize FLUKA random number generator."""
        if rng_state_file is None:
            rng_state_file = pathlib.Path(__file__).parent / "fluka_rng_state.dat"

        self._rng_state_file = rng_state_file
        self._logical_unit = 888

        # Load existing state or create new one
        pfile = pathlib.Path(self._rng_state_file)
        if pfile.exists() and pfile.stat().st_size > 0:
            self._lib.load_rng_state(self._rng_state_file, self._logical_unit)
        else:
            seed = 0 if seed is None else seed
            self._lib.init_rng_state(
                self._rng_state_file, self._logical_unit, seed, 0, 0
            )
            self._lib.save_rng_state(self._rng_state_file, self._logical_unit)

    def _set_quasielastic(self, enable):
        """Configure quasielastic interactions."""
        self._lib.qelcmm.lxsqel = 0  # default 0
        self._lib.qelcmm.lpqels = 1 if enable else 0  # default 1
        self._lib.nucflg.lqecmp = 1 if enable else 0  # default 1

    def _init_fluka_materials(self, evt_kin):
        """Initialize FLUKA materials and setup."""
        # Prepare material lists
        materials = []
        for mat in _DEFAULT_MATERIALS:
            p = process_particle(mat)
            if p not in materials:
                materials.append(p)

        # Add target from kinematics
        target = process_particle(evt_kin.p2)
        if isinstance(target, CompositeTarget):
            for comp in target.components:
                if comp not in materials:
                    materials.append(comp)
        elif target not in materials:
            materials.append(target)

        # Build arrays for FLUKA initialization
        nelements = np.array([1] * len(materials))
        charges = np.array([p.Z for p in materials])
        weights = np.array([1.0] * len(materials))

        # Setup arguments
        setup_args = [
            self._max_momentum_p,
            self._fluka_dpmjet_transition,
            self._transition_smearing,
            self._interaction_type,
            self._material_print,
        ]

        # Initialize FLUKA
        fluka_material_idcs = self._lib.chromo_stpxyz(
            nelements, charges, weights, *setup_args
        )

        # Store material mapping
        self._materials = materials
        self._material_idcs = fluka_material_idcs

    def _get_material_index(self, particle):
        """Get FLUKA material index for a particle."""
        p = process_particle(particle)
        try:
            idx = self._materials.index(p)
            return self._material_idcs[idx]
        except ValueError:
            msg = f"{p} is not among initialized target materials"
            raise KeyError(msg)

    def _cross_section(self, kin=None, max_info=False):
        """Calculate cross sections for given kinematics."""
        kin = self.kinematics if kin is None else kin
        projectile = self._fluka_pid(kin.p1)
        target = self._get_material_index(kin.p2)
        p_momentum = 0  # not used

        inel = self._lib.chromo_sgmxyz(
            projectile, target, kin.ekin, p_momentum, InteractionType.INELASTIC.value
        )
        el = self._lib.chromo_sgmxyz(
            projectile, target, kin.ekin, p_momentum, InteractionType.ELASTIC.value
        )

        return CrossSectionData(inelastic=float(inel), elastic=float(el))

    def _set_kinematics(self, kin):
        """Set kinematics for the simulation."""
        # Validate that target material is available
        self._get_material_index(kin.p2)

    def _set_stable(self, pdgid, stable):
        """Set particle stability (not implemented for FLUKA)."""
        info(2, f"Set_stable method no effect can't set {pdgid} stable to {stable}.")

    def _generate(self):
        """Generate a single event."""
        k = self.kinematics
        projectile = self._fluka_pid(k.p1)
        target = self._get_material_index(k.p2)

        self._lib.chromo_evtxyz(
            projectile,
            target,
            k.ekin,
            0,  # p_momentum = 0 means only kinetic energy will be used
            0,
            0,
            1,  # direction should be by default +z (not random)
            self._interaction_type,
        )

        return True

    def _cleanup_fort(self):
        import pathlib

        # Remove fort.* files created during fluka runs
        for fort_file in pathlib.Path(".").glob("fort.*"):
            fort_file.unlink()

        this_parent = pathlib.Path(__file__).parent
        lib_parent = pathlib.Path(self._lib.__file__).parent
        other_files = [
            this_parent / "fluka_rng_state.dat",
            lib_parent / "fluka_rng_state.dat",
            pathlib.Path(".") / ".timer.out",
        ]
        for f in other_files:
            f.unlink(missing_ok=True)

    def __del__(self):
        """Cleanup when the Fluka object is destroyed."""
        try:
            self._cleanup_fort()
        except Exception:
            # Silently ignore cleanup errors during destruction
            pass

    def _fluka_pid(self, pdg_id):
        """Convert PDG ID to FLUKA particle ID."""
        return self._lib.icode_from_pdg(pdg_id)

    def _index6(self, pdg_id):
        """Convert PDG ID to internal FLUKA array indices."""
        return self._lib.part.kptoip[self._lib.icode_from_pdg(pdg_id) + 6] + 6

    # Particle property methods
    def particle_name(self, pdg_id):
        """Get particle name from PDG ID."""
        return str(self._lib.chpprp.prname[self._index6(pdg_id)], encoding="utf-8")

    def particle_short_name(self, pdg_id):
        """Get particle short name from PDG ID."""
        return str(
            self._lib.chpart.aname[self._index6(pdg_id)], encoding="utf-8"
        ).strip()

    def particle_mass(self, pdg_id):
        """Get particle mass from PDG ID in GeV/c^2."""
        return self._lib.part.aam[self._index6(pdg_id)]

    def particle_tau(self, pdg_id):
        """Get particle lifetime from PDG ID in seconds."""
        return self._lib.part.tau[self._index6(pdg_id)]

    def particle_charge(self, pdg_id):
        """Get particle charge from PDG ID."""
        return self._lib.part.iich[self._index6(pdg_id)]

    def save_rng_state(self, file=None):
        """Save RNG state to file."""
        state_file = self._rng_state_file if file is None else file
        self._lib.save_rng_state(state_file, self._logical_unit)

    def load_rng_state(self, file=None):
        """Load RNG state from file."""
        state_file = self._rng_state_file if file is None else file
        self._lib.load_rng_state(state_file, self._logical_unit)

    def fluka_rand(self):
        """Generate a random number using the FLUKA RNG."""
        return self._lib.fluka_rand()
