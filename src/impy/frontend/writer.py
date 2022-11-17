import numpy as np
import contextlib
import gzip

BUFFER_SIZE = 100000
INT_TYPE = np.int32
FLOAT_TYPE = np.float32


class RootFile:
    def __init__(self, file):
        import uproot

        self._file = uproot.recreate(file)
        self._tree = None
        self._buffers = {
            "px": np.empty(BUFFER_SIZE, FLOAT_TYPE),
            "py": np.empty(BUFFER_SIZE, FLOAT_TYPE),
            "pz": np.empty(BUFFER_SIZE, FLOAT_TYPE),
            "m": np.empty(BUFFER_SIZE, FLOAT_TYPE),
            "vx": np.empty(BUFFER_SIZE, FLOAT_TYPE),
            "vy": np.empty(BUFFER_SIZE, FLOAT_TYPE),
            "vz": np.empty(BUFFER_SIZE, FLOAT_TYPE),
            "pid": np.empty(BUFFER_SIZE, INT_TYPE),
            "status": np.empty(BUFFER_SIZE, INT_TYPE),
            "parent": np.empty(BUFFER_SIZE, INT_TYPE),
        }
        self._lengths = []
        self._ievent = 0

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if self._ievent > 0:
            self._write_buffers()
        self._file.close()

    def _write_buffers(self):
        import awkward as ak

        lengths = self._lengths
        b = self._ievent
        chunk = {
            "": ak.zip(
                {
                    key: ak.unflatten(val[:b], lengths)
                    for (key, val) in self._buffers.items()
                }
            ),
        }
        if self._tree is None:
            self._file["event"] = chunk
            self._tree = self._file["event"]
        else:
            self._tree.extend(chunk)

    def write(self, event):
        b = self._ievent
        # apid = np.abs(event.pid)
        # for pdg in quarks_and_diquarks_and_gluons:
        #     mask &= apid != pdg
        # skip beam particles
        mask = np.ones(len(event), bool)
        mask[:2] = False
        event = event[mask]
        a = b
        b += len(event)
        self._lengths.append(b - a)
        for key, val in self._buffers.items():
            if key == "parent":
                val[a:b] = event.parents[:, 0] - 1
            else:
                val[a:b] = getattr(event, key)

        if b > BUFFER_SIZE // 2:
            self._write_buffers()
            self._lengths = []
            self._ievent = 0


@contextlib.contextmanager
def hepmc(file):
    from pyhepmc.io import _WrappedWriter, WriterAscii, pyiostream

    op = open
    if file.suffix == ".gz":
        op = gzip.open
    with op(file, "wb") as f:
        with pyiostream(f) as io:
            with _WrappedWriter(io, None, WriterAscii) as w:
                yield w


def lhe(file, mode):
    raise SystemExit("LHE not yet supported")
