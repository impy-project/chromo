import chromo
import numpy as np


class Pythia8DecayHandler:
    def __init__(
        self,
        seed=None,
        extra_stable_pids=None,
        extra_decaying_pids=None,
        stable_pids=None,
        decaying_pids=None,
    ):
        """
        Initialize the Pythia8DecayHandler with Pythia8 and set particle types for decay.

        Parameters:
            seed (int): Random seed for Pythia8 initialization.
            extra_stable_pids (list[int]): List of additional PDG IDs
            for stable particles.
            extra_decaying_pids (list[int]): List of additional PDG IDs
            for decaying particles.
            stable_pids (list[int]): List of PDG IDs
            for exclusively stable particles (all other are decaying if they can decay).
            decaying_pids (list[int]): List of PDG IDs
            for exclusively decaying particles (all other are stable).

        Notes:
            - If both `extra_stable_pids` and `extra_decaying_pids`
            are provided, check for duplicates.
            - If `stable_pids` or `decaying_pids`
            are provided, ensure exclusivity.
        """
        # Initialize Pythia8 with specific configurations
        config = ["ProcessLevel:all = off", "ParticleDecays:tau0Max = 1e100"]
        pythia8 = chromo.models.Pythia8(None, seed=seed, config=config, banner=False)
        self.pythia = pythia8._pythia

        self.set_stable_decaying(
            extra_stable_pids,
            extra_decaying_pids,
            stable_pids,
            decaying_pids,
        )

    def set_stable_decaying(
        self,
        extra_stable_pids=None,
        extra_decaying_pids=None,
        stable_pids=None,
        decaying_pids=None,
    ):
        self.extra_stable_pids = extra_stable_pids
        self.extra_decaying_pids = extra_decaying_pids
        self.stable_pids = stable_pids
        self.decaying_pids = decaying_pids
        # Check for duplicates in extra_stable_pids and extra_decaying_pids
        if (extra_stable_pids is not None) and (extra_decaying_pids is not None):
            dup_pids = set(extra_stable_pids) & set(extra_decaying_pids)
            if dup_pids:
                raise ValueError(
                    f"Duplicates found in {dup_pids} between "
                    "extra_stable_pids and extra_decaying_pids"
                )

        # Set particles as stable or decaying based on the provided lists
        if extra_stable_pids is not None:
            for pid in extra_stable_pids:
                self.pythia.particleData.mayDecay(pid, False)

        if extra_decaying_pids is not None:
            for pid in extra_decaying_pids:
                self.pythia.particleData.mayDecay(pid, True)

        # Handle exclusive stability options
        if stable_pids is not None:
            self._set_exclusive_stability(stable_pids, False)

        if decaying_pids is not None:
            self._set_exclusive_stability(decaying_pids, True)

        # Store stable and decaying PDG IDs for future reference
        all_stable_pids = []
        all_decaying_pids = []
        for pentry in self.pythia.particleData.all():
            if pentry.mayDecay and pentry.canDecay:
                if pentry.hasAnti:
                    all_decaying_pids.append(pentry.antiId)
                all_decaying_pids.append(pentry.id)
            else:
                if pentry.hasAnti:
                    all_stable_pids.append(pentry.antiId)
                all_stable_pids.append(pentry.id)

        self.all_decaying_pids = np.array(all_decaying_pids, dtype=np.int64)
        self.all_stable_pids = np.array(all_stable_pids, dtype=np.int64)

    def _set_exclusive_stability(self, pids, is_decaying):
        """
        Set particles exclusively as stable or decaying based on the provided list.

        Parameters:
            pids (list[int]): List of PDG IDs for particles to be exclusively set.
            is_decaying (bool): True if the particles should be set as decaying,
            False for stable.

        Raises:
            ValueError: If exclusivity conditions are not met.
        """
        if (
            sum(
                x is not None
                for x in [
                    self.extra_stable_pids,
                    self.extra_decaying_pids,
                    self.stable_pids,
                    self.decaying_pids,
                ]
            )
            != 1
        ):
            raise ValueError(
                "stable_pids or decaying_pids can be set "
                "if all other options are None"
            )

        # Set particles exclusively as stable or decaying
        for pentry in self.pythia.particleData.all():
            self.pythia.particleData.mayDecay(pentry.id, not is_decaying)
        for pid in pids:
            self.pythia.particleData.mayDecay(pid, is_decaying)

    def __call__(self, event):
        """
        Decay particles in the provided `event` that are marked as `decaying`.

        After the call, `event` contains both the initial particles and
        the newly produced particles through decay.
        The decayed particles have their `event.status` set to 2.

        Parameters:
            event (Event): The event containing particles to be decayed.
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
        event.charge = self.pythia.charge()
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
