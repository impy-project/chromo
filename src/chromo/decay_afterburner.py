import chromo
import numpy as np


class DecayAfterburner:
    def __init__(self, seed=None, stable_pids=None, decaying_pids=None):
        """
        Initialize the DecayAfterburner with Pythia8 and set particle types for decay.

        Parameters:
            seed (int): Random seed for Pythia8 initialization.
            stable_pids (list[int]): A list of PDG IDs corresponding to particles
                                     that should be considered stable.
            decaying_pids (list[int]): A list of PDG IDs corresponding to particles
                                       that should be considered decaying.

        Notes:
            - If both `stable_pids` and `decaying_pids` are `None`,
              Pythia8 defaults are applied.
            - If the same PDG ID appears in both `stable_pids` and `decaying_pids`,
              `decaying_pids` takes precedence over `stable_pids`.
        """
        # Initialize Pythia8 with specific configurations
        config = ["ProcessLevel:all = off", "ParticleDecays:tau0Max = 1e100"]
        pythia8 = chromo.models.Pythia8(None, seed=seed, config=config, banner=False)
        self.pythia = pythia8._pythia

        # Set particles as stable or decaying based on the provided lists
        if stable_pids is not None:
            for pid in stable_pids:
                self.pythia.particleData.mayDecay(pid, False)

        if decaying_pids is not None:
            for pid in decaying_pids:
                self.pythia.particleData.mayDecay(pid, True)

        # Store the provided stable and decaying PDG IDs for future reference
        self.stable_pids = stable_pids
        self.decaying_pids = decaying_pids

    def __call__(self, event):
        """
        Decay particles in the provided `event` that are marked as `decaying`.

        A particle is considered `decaying` if its PDG ID is present in the
        `decaying_pids` list provided during the constructor's initialization
        or if it is set to decay by default in Pythia8.
        The `stable_pids` specified during the constructor's initialization are
        used to set particles as stable.

        After the call, `event` contains both the initial particles and
        the newly produced particles through decay.
        The decayed particles have their `event.status` set to 2.

        Parameters:
            event (Event): The event containing particles to be decayed.

        Raises:
            RuntimeError: If any produced particles have an unknown PDG ID.
        """
        # Return if the event is empty
        if len(event.pid) == 0:
            return

        # Save the original status of non-final particles before decay
        # as they will become 2
        not_final = event.status != 1
        status = np.copy(event.status[not_final])

        # Particles with negative status are not processed by Pythia
        event.status[not_final] = -1

        # Save refs to original event parameters
        # that are not passed to Pythia8 stack
        charge = event.charge
        vx = event.vx
        vy = event.vy
        vz = event.vz
        vt = event.vt

        mothers = event.mothers if event.mothers is not None else None
        daughters = event.daughters if event.daughters is not None else None
        init_len = len(event.pid)

        # Put the event on the Pythia8 stack for decay
        self.pythia.event.fill(
            event.pid, event.status, event.px, event.py, event.pz, event.en, event.m
        )

        # Decay the event using Pythia8
        self.pythia.forceHadronLevel()

        # Get the decayed event results from Pythia8
        pevent = self.pythia.event
        event.pid = pevent.pid()
        event.status = pevent.status()
        event.charge, not_valid = self.pythia.charge_nothrow()
        if np.any(not_valid[init_len:]):
            raise RuntimeError("Produced particles have unknown to Pythia PDG ID")

        event.px = pevent.px()
        event.py = pevent.py()
        event.pz = pevent.pz()
        event.en = pevent.en()
        event.m = pevent.m()
        event.vx = pevent.vx()
        event.vy = pevent.vy()
        event.vz = pevent.vz()
        event.vt = pevent.vt()
        event.mothers = np.maximum(pevent.mothers() - 1, -1)
        event.daughters = np.maximum(pevent.daughters() - 1, -1)

        # Restore the original parameters to the event
        event.status[np.where(not_final)[0]] = status
        event.charge[0:init_len] = charge
        event.vx[0:init_len] = vx
        event.vy[0:init_len] = vy
        event.vz[0:init_len] = vz
        event.vt[0:init_len] = vt

        if mothers is not None:
            event.mothers[0:init_len] = mothers
        if daughters is not None:
            event.daughters[0:init_len] = daughters
