import numpy as np
from impy.constants import quarks_and_diquarks_and_gluons

BUFFER_SIZE = 100000
INT_TYPE = np.int32
FLOAT_TYPE = np.float32


# Differences to CRMC
# - No header yet
# - Particle tree
#   - Branch n instead of nPart, name cannot be chosen in uproot
#   - ImpactParameter renamed to impact
#   - Branch E is redundant, we skip this to save space
#   - Extra branches vx, vy, vz, parent
class Root:
    def __init__(self, file):
        import uproot

        # TODO also write header once
        # Seed
        # ProjectileId
        # ProjectileMomentum
        # TargetId
        # TargetMomentum
        # HEModel
        # sigmaPairTot
        # sigmaPairInel
        # sigmaPairEl
        # sigmaTot
        # sigmaInel
        # sigmaEl

        self._file = uproot.recreate(file)
        self._tree = None
        self._event_buffers = {
            "impact": np.empty(BUFFER_SIZE, FLOAT_TYPE),
        }
        self._particle_buffers = {
            "px": np.empty(BUFFER_SIZE, FLOAT_TYPE),
            "py": np.empty(BUFFER_SIZE, FLOAT_TYPE),
            "pz": np.empty(BUFFER_SIZE, FLOAT_TYPE),
            "m": np.empty(BUFFER_SIZE, FLOAT_TYPE),
            "vx": np.empty(BUFFER_SIZE, FLOAT_TYPE),
            "vy": np.empty(BUFFER_SIZE, FLOAT_TYPE),
            "vz": np.empty(BUFFER_SIZE, FLOAT_TYPE),
            "pdgid": np.empty(BUFFER_SIZE, INT_TYPE),
            "status": np.empty(BUFFER_SIZE, INT_TYPE),
            "parent": np.empty(BUFFER_SIZE, INT_TYPE),
        }
        self._lengths = []
        self._iparticle = 0

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if self._iparticle > 0:
            self._write_buffers()
        self._file.__exit__(exc_type, exc_value, traceback)

    def _write_buffers(self):
        import awkward as ak

        lengths = self._lengths
        b = self._iparticle
        chunk = {key: val[: len(lengths)] for (key, val) in self._event_buffers.items()}
        chunk[""] = ak.zip(
            {
                key: ak.unflatten(val[:b], lengths)
                for (key, val) in self._particle_buffers.items()
            }
        )
        if self._tree is None:
            self._file["event"] = chunk
            self._tree = self._file["event"]
        else:
            self._tree.extend(chunk)
        self._iparticle = 0
        self._lengths = []

    def write(self, event):
        mask = True
        # skip parton shower
        apid = np.abs(event.pid)
        for pdg in quarks_and_diquarks_and_gluons:
            mask &= apid != pdg
        # skip beam particles
        mask[:2] = False
        event = event[mask]

        ievent = len(self._lengths)
        self._event_buffers["impact"][ievent] = getattr(event, "impact_parameter", 0.0)

        a = self._iparticle
        self._iparticle += len(event)
        b = self._iparticle
        self._lengths.append(b - a)
        for key, val in self._particle_buffers.items():
            if key == "parent":
                val[a:b] = event.parents[:, 0] - 1
            elif key == "pdgid":
                val[a:b] = event.pid
            else:
                val[a:b] = getattr(event, key)

        if b > BUFFER_SIZE // 2:
            self._write_particle_buffers()
            self._lengths = []
            self._iparticle = 0


def lhe(file, mode):
    raise SystemExit("LHE not yet supported")
