"""
Pure Python implementation of PythiaCascade functionality.

This module provides a Python implementation of the core PythiaCascade
nextColl functionality for hadron-nucleus interactions, equivalent to
the C++ PythiaCascade class but implemented entirely in Python using
chromo's existing Pythia8 integration.
"""

import warnings
import numpy as np
from particle import Particle

from chromo.common import CrossSectionData, EventData, MCRun
from chromo.constants import standard_projectiles
from chromo.kinematics import EventFrame
from chromo.util import Nuclei
from chromo.models.pythia8 import Pythia8
from chromo.decay_handler import Pythia8DecayHandler


class PythiaCascadePythonEvent(EventData):
    """Event wrapper for Python PythiaCascade implementation."""

    def __init__(self, generator):
        """
        Parameters
        ----------
        generator : PythiaCascadePython
            The generator instance
        """
        # Get the current event from the generator
        event_data = generator._current_event_data

        super().__init__(
            (generator.name, generator.version),
            generator.kinematics,
            generator.nevents,
            None,  # impact parameter not available in cascade mode
            None,  # n_wounded not available in cascade mode
            generator._inel_or_prod_cross_section,
            event_data["pid"],
            event_data["status"],
            event_data["charge"],
            event_data["px"],
            event_data["py"],
            event_data["pz"],
            event_data["en"],
            event_data["m"],
            event_data["vx"],
            event_data["vy"],
            event_data["vz"],
            event_data["vt"],
            event_data["mothers"],
            event_data["daughters"],
        )


class PythiaCascadePython(MCRun):
    """
    Pure Python implementation of PythiaCascade for hadron-nucleus interactions.

    This class provides the core functionality of the C++ PythiaCascade class
    but implemented entirely in Python using chromo's existing Pythia8 integration.
    It focuses on the nextColl functionality for generating hadron-nucleus collisions.
    """

    _name = "PythiaCascadePython"
    _version = "1.0"
    _library_name = "_pythia8"  # Reuse Pythia8 library
    _event_class = PythiaCascadePythonEvent
    _frame = EventFrame.FIXED_TARGET

    # Support nucleus types for cascade interactions
    _projectiles = standard_projectiles
    _targets = standard_projectiles | Nuclei()

    # Nucleus data tables from PythiaCascade.h
    _tab_A = np.array([1, 2, 4, 9, 12, 14, 16, 27, 40, 56, 63, 84, 107, 129, 197, 208])
    _tab_offset = np.array([0., 0.03, 0.08, 0.15, 0.20, 0.20, 0.20, 0.26, 
                           0.30, 0.34, 0.40, 0.40, 0.40, 0.50, 0.50, 0.60])
    _tab_slope = np.array([0., 0.0016, 0.0033, 0.0075, 0.0092, 0.0105, 0.012, 
                          0.017, 0.022, 0.027, 0.028, 0.034, 0.040, 0.044, 0.055, 0.055])
    _tab_slope_lo = np.array([0., 0.0031, 0.0073, 0.015, 0.0192, 0.0205, 0.022, 
                             0.03, 0.037, 0.044, 0.048, 0.054, 0.06, 0.069, 0.08, 0.085])
    _tab_border = 20.0

    # Constants
    _prob_sd = 0.3  # Fraction of single diffractive events beyond first collision
    _e_kin_min = 0.2  # Minimum kinetic energy threshold (GeV)
    _mp = 0.938272  # Proton mass (GeV)

    def __init__(
        self,
        evt_kin,
        *,
        seed=None,
        config=None,
        banner=True,
        rapid_decays=True,
        small_tau0=1e-10,
        list_final=False,
    ):
        """
        Parameters
        ----------
        rapid_decays : bool, optional
            Whether to decay particles with tau0 < small_tau0. Default: True
        small_tau0 : float, optional
            Lifetime threshold for rapid decays (mm/c). Default: 1e-10
        list_final : bool, optional
            Whether to return only final particles. Default: False
        """
        if evt_kin is None:
            raise ValueError("evt_kin is required for PythiaCascadePython")

        if evt_kin.frame != EventFrame.FIXED_TARGET:
            raise ValueError("PythiaCascadePython requires fixed target frame kinematics")

        super().__init__(seed)

        # Store cascade-specific parameters
        self._rapid_decays = rapid_decays
        self._small_tau0 = small_tau0
        self._list_final = list_final

        # Set default config for cascade mode
        if config is None:
            config = ["SoftQCD:all = on", "LowEnergyQCD:all = on"]

        # Create collision Pythia instance for individual collisions
        coll_config = [
            *config,
            "Beams:allowVariableEnergy = on",
            "Beams:allowIDAswitch = on", 
            "Beams:frameType = 3",  # Fixed target
            "Beams:pzB = 0.",
            "HadronLevel:Decay = off",  # Decays handled separately
            "Print:quiet = on",
            "Check:epTolErr = 0.01", 
            "Check:epTolWarn = 0.0001",
            "Check:mTolErr = 0.01",
            "Stat:showProcessLevel = off",
            "Stat:showPartonLevel = off",
        ]

        self._pythia_coll = Pythia8(None, seed=seed, config=coll_config, banner=banner)

        # Initialize kinematics
        self.kinematics = evt_kin

        # Initialize cross section and event cache
        self._sigma_now = 0.0
        self._ecm_now = 0.0
        self._current_event_data = None

        # Set up decay handler if needed
        if self._rapid_decays:
            from chromo.constants import long_lived
            self._decay_handler = Pythia8DecayHandler(long_lived, seed=seed)
        else:
            self._decay_handler = None

    def _n_coll_avg(self, A):
        """
        Calculate the average number of collisions for nucleus A using numpy interpolation.
        
        This implements the interpolation logic from PythiaCascade::nCollAvg().
        """
        if A < 1 or A > 208:
            return np.nan

        # Use the appropriate slope arrays based on sigma threshold
        if self._sigma_now < self._tab_border:
            slopes = self._tab_slope_lo
            offsets = np.zeros_like(self._tab_offset)
        else:
            slopes = self._tab_slope
            offsets = self._tab_offset

        # Calculate n_coll values for tabulated masses
        n_coll_values = 1.0 + offsets + slopes * self._sigma_now

        # Interpolate or extrapolate to get value for mass A
        if A <= self._tab_A[-1]:
            # Interpolate within tabulated range
            n_coll = np.interp(A, self._tab_A, n_coll_values)
        else:
            # Extrapolate beyond tabulated range
            scaling = (A / self._tab_A[-1]) ** (2.0/3.0)
            n_coll = 1.0 + scaling * (n_coll_values[-1] - 1.0)

        return n_coll

    def sigma_setup_hN(self):
        """
        Calculate hadron-proton collision cross section using event kinematics.
        
        Equivalent to PythiaCascade::sigmaSetuphN().
        """
        kin = self.kinematics
        
        # Get beam momentum and particle info from kinematics
        p_beam = kin.beams[0]  # 4-momentum of beam particle
        projectile_id = int(kin.p1)
        projectile_mass = kin.m1
        
        # Calculate kinetic energy
        e_kinetic = p_beam[3] - projectile_mass
        
        # Cannot handle low-energy hadrons
        if e_kinetic < self._e_kin_min:
            return False

        # Calculate CM energy with proton at rest
        p_target = np.array([0, 0, 0, self._mp])
        s = np.sum((p_beam + p_target)**2) - np.sum((p_beam[:3] + p_target[:3])**2)
        self._ecm_now = np.sqrt(s)

        # Get cross section from Pythia
        try:
            # Set up collision kinematics
            if not self._pythia_coll._pythia.setBeamIDs(projectile_id, 2212):
                return False
                
            if not self._pythia_coll._pythia.setKinematics(p_beam, p_target):
                return False
            
            # Get total cross section
            self._sigma_now = self._pythia_coll._pythia.getSigmaTotal(
                projectile_id, 2212, self._ecm_now, projectile_mass, self._mp
            )
            
            if self._sigma_now <= 0:
                if self._ecm_now - projectile_mass - self._mp > self._e_kin_min:
                    warnings.warn("Vanishing cross section", RuntimeWarning)
                return False
                
        except Exception as e:
            warnings.warn(f"Failed to calculate cross section: {e}", RuntimeWarning)
            return False

        return True

    def sigma_hA(self, A):
        """
        Calculate hadron-nucleus cross section.
        
        Equivalent to PythiaCascade::sigmahA().
        """
        if A < 1 or A > 208:
            warnings.warn("A is outside of valid range (1 <= A <= 208)", RuntimeWarning)
            return 0.0

        # Correction factor for number of h-nucleon collisions per h-nucleus one
        n_coll_avg = self._n_coll_avg(A)
        return A * self._sigma_now / n_coll_avg

    def next_coll(self, Z, A, vertex=None):
        """
        Generate a hadron-nucleus collision.
        
        Equivalent to PythiaCascade::nextColl().
        
        Parameters
        ----------
        Z : int
            Atomic number of target nucleus
        A : int
            Mass number of target nucleus
        vertex : array-like, optional
            Collision vertex (x, y, z, t). Default: origin
            
        Returns
        -------
        dict
            Event data dictionary with particle information
        """
        if vertex is None:
            vertex = np.array([0.0, 0.0, 0.0, 0.0])
        elif len(vertex) != 4:
            raise ValueError("Vertex array must have 4 components (x, y, z, t)")

        if A < 1 or A > 208:
            warnings.warn("A is outside of valid range (1 <= A <= 208)", RuntimeWarning)
            return self._create_empty_event()

        kin = self.kinematics
        
        # Get beam particle info from kinematics
        p_beam = kin.beams[0]  # 4-momentum of beam particle
        projectile_id = int(kin.p1)
        
        # Initialize event record with system particle and incoming hadron
        event_data = {
            "pid": [90, projectile_id],
            "status": [-11, 12],
            "charge": [0.0, self._get_charge(projectile_id)],
            "px": [p_beam[0], p_beam[0]],
            "py": [p_beam[1], p_beam[1]], 
            "pz": [p_beam[2], p_beam[2]],
            "en": [p_beam[3], p_beam[3]],
            "m": [kin.m1, kin.m1],
            "vx": [vertex[0], vertex[0]],
            "vy": [vertex[1], vertex[1]],
            "vz": [vertex[2], vertex[2]],
            "vt": [vertex[3], vertex[3]],
            "mothers": [np.array([0, 0]), np.array([0, 0])],
            "daughters": [np.array([0, 0]), np.array([0, 0])],
        }

        # Set up for collisions on nucleus
        np_nucleons = Z  # Number of protons
        nn_nucleons = A - Z  # Number of neutrons
        
        # Direction vector for finding leading particle
        dir_now = p_beam[:3] / np.linalg.norm(p_beam[:3])

        # Drop rate of geometric series for multiple collisions
        prob_more = 1.0 - 1.0 / self._n_coll_avg(A)

        # Loop over varying number of hit nucleons in target nucleus
        i_had = 1  # Index of initial hadron
        
        for i_coll in range(1, A + 1):
            if i_coll > 1 and np.random.random() > prob_more:
                break

            # Pick incoming projectile
            if i_coll == 1:
                i_proj = i_had
                proc_type = 0
            else:
                # Find highest-pLongitudinal particle from latest subcollision
                i_proj = self._find_leading_particle(event_data, dir_now)
                if i_proj == 0:
                    break
                    
                # Choose process type for subsequent collisions
                proj_p = np.array([event_data["px"][i_proj], event_data["py"][i_proj], 
                                 event_data["pz"][i_proj], event_data["en"][i_proj]])
                target_p = np.array([0, 0, 0, self._mp])
                s = np.sum((proj_p + target_p)**2) - np.sum((proj_p[:3] + target_p[:3])**2)
                ecm_sub = np.sqrt(s)
                
                proc_type = 4 if ecm_sub > 10.0 and np.random.random() < self._prob_sd else 1

            # Pick target nucleon (proton or neutron)
            do_proton = np.random.random() < (np_nucleons / (np_nucleons + nn_nucleons))
            if do_proton:
                np_nucleons -= 1
                id_nuc = 2212
            else:
                nn_nucleons -= 1  
                id_nuc = 2112

            # Perform projectile-nucleon subcollision
            try:
                proj_momentum = np.array([event_data["px"][i_proj], event_data["py"][i_proj],
                                        event_data["pz"][i_proj], event_data["en"][i_proj]])
                collision_result = self._perform_collision(
                    event_data["pid"][i_proj], id_nuc, proj_momentum, proc_type
                )
                
                if collision_result is None:
                    continue

                # Add target nucleon to event record
                self._append_particle(event_data, id_nuc, -181 if i_coll == 1 else -182,
                                    [0, 0, 0, self._mp if id_nuc == 2212 else 0.940], 
                                    self._mp if id_nuc == 2212 else 0.940, vertex,
                                    mothers=[0, i_proj])

                # Update system energy
                nucleon_mass = self._mp if id_nuc == 2212 else 0.940
                event_data["en"][0] += nucleon_mass
                event_data["m"][0] = np.sqrt(event_data["en"][0]**2 - 
                                           (event_data["px"][0]**2 + event_data["py"][0]**2 + event_data["pz"][0]**2))

                # Add secondary particles from collision
                start_idx = len(event_data["pid"])
                for particle in collision_result:
                    self._append_particle(
                        event_data, particle["pid"], particle["status"],
                        [particle["px"], particle["py"], particle["pz"], particle["en"]],
                        particle["m"], vertex, mothers=[len(event_data["pid"]) - 1, i_proj]
                    )
                
                # Update daughters info
                if len(event_data["pid"]) > start_idx:
                    event_data["daughters"][i_proj] = np.array([start_idx, len(event_data["pid"]) - 1])
                    event_data["daughters"][-2] = np.array([start_idx, len(event_data["pid"]) - 1])  # nucleon daughters
                
                # Mark projectile as interacted
                if event_data["status"][i_proj] > 0:
                    event_data["status"][i_proj] = -event_data["status"][i_proj]

            except Exception as e:
                warnings.warn(f"Collision failed: {e}", RuntimeWarning)
                continue

        # Apply rapid decays using DecayHandler if configured
        if self._rapid_decays and self._decay_handler:
            # Convert to EventData format for decay handler
            temp_event = type('obj', (object,), {
                'pid': np.array(event_data["pid"]),
                'status': np.array(event_data["status"]),
                'charge': np.array(event_data["charge"]),
                'px': np.array(event_data["px"]),
                'py': np.array(event_data["py"]),
                'pz': np.array(event_data["pz"]),
                'en': np.array(event_data["en"]),
                'm': np.array(event_data["m"]),
                'vx': np.array(event_data["vx"]),
                'vy': np.array(event_data["vy"]),
                'vz': np.array(event_data["vz"]),
                'vt': np.array(event_data["vt"]),
                'mothers': np.array(event_data["mothers"]),
                'daughters': np.array(event_data["daughters"]),
            })()
            
            # Apply decays
            self._decay_handler(temp_event)
            
            # Convert back to dict format
            event_data = {
                "pid": temp_event.pid,
                "status": temp_event.status,
                "charge": temp_event.charge,
                "px": temp_event.px,
                "py": temp_event.py,
                "pz": temp_event.pz,
                "en": temp_event.en,
                "m": temp_event.m,
                "vx": temp_event.vx,
                "vy": temp_event.vy,
                "vz": temp_event.vz,
                "vt": temp_event.vt,
                "mothers": temp_event.mothers,
                "daughters": temp_event.daughters,
            }

        # Optionally compress event record to show only final particles
        if self._list_final:
            event_data = self._compress_event(event_data)

        return event_data

    def _perform_collision(self, id_proj, id_target, p_proj, proc_type=0):
        """Perform individual hadron-nucleon collision using Pythia."""
        try:
            # Set up collision in Pythia
            if not self._pythia_coll._pythia.setBeamIDs(int(id_proj), int(id_target)):
                return None
                
            p_target = np.array([0, 0, 0, self._mp if id_target == 2212 else 0.940])
            if not self._pythia_coll._pythia.setKinematics(p_proj, p_target):
                return None

            # Generate collision
            if not self._pythia_coll._pythia.next(proc_type):
                return None

            # Extract final state particles using list comprehension
            event = self._pythia_coll._pythia.event
            particles = [
                {
                    "pid": event[i].id(),
                    "status": 1,  # Final state
                    "px": event[i].px(),
                    "py": event[i].py(), 
                    "pz": event[i].pz(),
                    "en": event[i].e(),
                    "m": event[i].m(),
                }
                for i in range(3, event.size()) if event[i].isFinal()
            ]
            
        except Exception as e:
            warnings.warn(f"Collision generation failed: {e}", RuntimeWarning)
            return None
        else:
            return particles

    def _append_particle(self, event_data, pid, status, momentum, mass, vertex, mothers=None):
        """Append a particle to the event record."""
        event_data["pid"].append(pid)
        event_data["status"].append(status)
        event_data["charge"].append(self._get_charge(pid))
        event_data["px"].append(momentum[0])
        event_data["py"].append(momentum[1])
        event_data["pz"].append(momentum[2])
        event_data["en"].append(momentum[3])
        event_data["m"].append(mass)
        event_data["vx"].append(vertex[0])
        event_data["vy"].append(vertex[1])
        event_data["vz"].append(vertex[2])
        event_data["vt"].append(vertex[3])
        
        if mothers is None:
            mothers = [0, 0]
        event_data["mothers"].append(np.array(mothers))
        event_data["daughters"].append(np.array([0, 0]))

    def _perform_rapid_decays(self, event_data):
        """Perform rapid decays of short-lived particles."""
        # This is a simplified implementation
        # In practice, you would use the main Pythia instance to handle decays
        return event_data

    def _compress_event(self, event_data):
        """Compress event record to show only final particles."""
        if not self._list_final:
            return event_data
            
        compressed = {key: [] for key in event_data.keys()}
        compressed["mothers"] = []
        compressed["daughters"] = []
        
        # Keep line 0 (system particle) and final particles
        for i, status in enumerate(event_data["status"]):
            if i == 0 or status > 0:
                for key in event_data.keys():
                    if key in ["mothers", "daughters"]:
                        if i == 0:
                            compressed[key].append(event_data[key][i])
                        else:
                            compressed[key].append(np.array([0, 0]))  # Remove history
                    else:
                        compressed[key].append(event_data[key][i])
        
        return compressed

    def _create_empty_event(self):
        """Create an empty event record."""
        return {
            "pid": np.array([]),
            "status": np.array([]),
            "charge": np.array([]),
            "px": np.array([]),
            "py": np.array([]),
            "pz": np.array([]),
            "en": np.array([]),
            "m": np.array([]),
            "vx": np.array([]),
            "vy": np.array([]),
            "vz": np.array([]),
            "vt": np.array([]),
            "mothers": np.array([]),
            "daughters": np.array([]),
        }

    def _get_charge(self, pid):
        """Get electric charge of particle from PDG ID."""
        # Simplified charge calculation - in practice use particle data
        charge_map = {
            2212: 1.0,   # proton
            2112: 0.0,   # neutron  
            211: 1.0,    # pi+
            -211: -1.0,  # pi-
            111: 0.0,    # pi0
            321: 1.0,    # K+
            -321: -1.0,  # K-
            130: 0.0,    # K_L
            310: 0.0,    # K_S
            13: -1.0,    # mu-
            -13: 1.0,    # mu+
            11: -1.0,    # e-
            -11: 1.0,    # e+
            22: 0.0,     # photon
            12: 0.0,     # nu_e
            14: 0.0,     # nu_mu
        }
        return charge_map.get(pid, 0.0)

    def _is_hadron(self, pid):
        """Check if particle is a hadron."""
        abs_pid = abs(pid)
        # Simplified hadron identification
        return (abs_pid > 100 and abs_pid != 130 and abs_pid != 310 and 
                abs_pid not in [11, 12, 13, 14, 15, 16, 22])

    def _set_kinematics(self, kin):
        """Set up kinematics for cascade collision."""
        if kin.frame != EventFrame.FIXED_TARGET:
            raise ValueError("PythiaCascadePython requires fixed target frame kinematics")

    def _cross_section(self, kin=None, max_info=False):
        """Calculate cross sections for hadron-nucleus interaction."""
        if kin is None:
            kin = self.kinematics

        # Calculate hadron-nucleon cross section using kinematics
        if not self.sigma_setup_hN():
            warnings.warn("Could not calculate hadron-nucleon cross section", RuntimeWarning)
            return CrossSectionData(total=0, inelastic=0, elastic=0)

        # Get hadron-nucleus cross section
        target_A = getattr(kin.p2, 'A', 1)
        sigma_hA = self.sigma_hA(target_A)

        return CrossSectionData(
            total=sigma_hA,
            inelastic=sigma_hA,  # Assume all inelastic for cascade
            elastic=0,
        )

    def _generate(self):
        """Generate a hadron-nucleus collision event."""
        kin = self.kinematics
        
        # Set up cross section calculation using kinematics
        if not self.sigma_setup_hN():
            return False

        # Generate collision
        target_Z = getattr(kin.p2, 'Z', 1)
        target_A = getattr(kin.p2, 'A', 1)
        
        try:
            self._current_event_data = self.next_coll(target_Z, target_A)
            return len(self._current_event_data["pid"]) > 0
        except Exception as e:
            warnings.warn(f"Event generation failed: {e}", RuntimeWarning)
            return False

    def _set_stable(self, pdgid, stable):
        """Set particle stability."""
        self._pythia_coll._pythia.particleData.mayDecay(pdgid, not stable)

    def _get_stable(self):
        """Get list of stable particles."""
        r = set()
        for p in self._pythia_coll._pythia.particleData.all():
            if p.tau0 > 1e-5 and not p.mayDecay:
                r.add(p.id)
                if p.hasAnti:
                    r.add(-p.id)
        return r

    def stat(self):
        """Print statistics."""
        import logging
        logging.info("Statistics from PythiaCascadePython:")
        logging.info("Collision Pythia instance:")
        self._pythia_coll.stat()
