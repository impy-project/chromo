import chromo
import numpy as np


class DecayAfterburner:
    def __init__(self, seed=None, stable_pids=None, decaying_pids=None):
        """Initialize pythia8 and set `stable` and `decaying` particle types.

        `decaying_pids` overrides `stable_pids` if the same pid is given
        in both lists. If both are `None`, Pythia8 defaults are applied
        """
        config = ["ProcessLevel:all = off", "ParticleDecays:tau0Max = 1e100"]
        pythia8 = chromo.models.Pythia8(None, seed=seed, config=config, banner=False)
        self.pythia = pythia8._pythia

        if stable_pids is not None:
            for pid in stable_pids:
                self.pythia.particleData.mayDecay(pid, False)

        if decaying_pids is not None:
            for pid in decaying_pids:
                self.pythia.particleData.mayDecay(pid, True)

        self.stable_pids = stable_pids
        self.decaying_pids = decaying_pids

    def __call__(self, event):
        """Decay particles in `event` which are set as `decaying`.

        Particle is `decaying` if its pdg id is passed
        in `decaying_pids` list of constructor or it is decaing by default
        in Pythia8. `stable_pids` of constructor is used to set stable
        particles.

        After the call `event` contains initial particles + produced particles
        The decayed particles has event.status == 2.
        """
        if len(event.pid) == 0:
            return

        # Save not_final status
        not_final = event.status != 1
        status = np.copy(event.status[not_final])

        # Particles with negative status are not processed by Pythia
        event.status[not_final] = -1

        # Save original event parameters which are not passed
        # to Pythia8 stack
        charge = event.charge
        vx = event.vx
        vy = event.vy
        vz = event.vz
        vt = event.vt

        mothers = event.mothers if event.mothers is not None else None
        daughters = event.daughters if event.daughters is not None else None
        init_len = len(event.pid)

        # Put on stack
        self.pythia.event.fill(
            event.pid, event.status, event.px, event.py, event.pz, event.en, event.m
        )

        # Decay
        self.pythia.forceHadronLevel()

        # Get results
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

        # Record back original parameters
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
