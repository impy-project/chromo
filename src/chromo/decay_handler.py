import chromo
import numpy as np


class Pythia8DecayHandler:
    def __init__(
        self,
        stable_pids,
        seed=None,
    ):
        """
        Initialize the Pythia8DecayHandler with Pythia8 and set particle types for decay.

        Parameters:
            seed (int): Random seed for Pythia8 initialization.
            stable_pids (list[int]): List of PDG IDs
            for stable particles (all other will decay).

        """
        # Initialize Pythia8 with specific configurations
        config = ["ProcessLevel:all = off"]
        pythia8 = chromo.models.Pythia8(None, seed=seed, config=config, banner=False)
        self.pythia = pythia8._pythia

        self.set_stable(stable_pids)

    def set_stable(self, stable_pids):
        # Set all particles as unstable
        for pentry in self.pythia.particleData.all():
            self.pythia.particleData.mayDecay(pentry.id, True)

        # Set stable_pids as stable
        for pid in stable_pids:
            self.pythia.particleData.mayDecay(pid, False)

        # Store stable and unstable PDG IDs for future reference
        all_stable_pids = []
        all_unstable_pids = []
        for pentry in self.pythia.particleData.all():
            if pentry.mayDecay and pentry.canDecay:
                if pentry.hasAnti:
                    all_unstable_pids.append(pentry.antiId)
                all_unstable_pids.append(pentry.id)
            else:
                if pentry.hasAnti:
                    all_stable_pids.append(pentry.antiId)
                all_stable_pids.append(pentry.id)

        self.all_unstable_pids = np.array(all_unstable_pids, dtype=np.int64)
        self.all_stable_pids = np.array(all_stable_pids, dtype=np.int64)

    def __call__(self, event):
        """
        Decay particles in the provided `event` that are not set as `stable`.

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
