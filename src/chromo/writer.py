import numpy as np
from chromo.constants import quarks_and_diquarks_and_gluons, millibarn, GeV
from chromo.kinematics import CompositeTarget
import dataclasses
from pathlib import Path
from abc import ABC, abstractmethod

INT_TYPE = np.int32
FLOAT_TYPE = np.float32


def _raise_import_error(name, task):
    raise ModuleNotFoundError(
        f"{name} not found, please install {name} (`pip install {name}`) to {task}"
    )


class Writer(ABC):
    @abstractmethod
    def __init__(self, file, model, **kwargs): ...

    @abstractmethod
    def write(self, event): ...

    # override as needed
    def __enter__(self):
        return self

    # override as needed
    def __exit__(self, *args):
        return


class Null(Writer):
    """
    Null writer.

    It does not do anything. Useful for benchmarks.
    """

    def __init__(self, file, model):
        pass

    def write(self, event):
        pass


# Root writer: Differences to CRMC
#
# - Names of trees and branches are in snake_case instead of CamelCase
#
# - Header tree was removed
#   - Metadata are stored in YAML format in title of particle tree
#   - Renamed branches
#     - HEModel -> model
#     - sigmaTot -> sigma_total
#     - sigmaEl -> sigma_elastic
#     - sigmaInel -> sigma_inelastic
#   - Branches sigmaPair* not included
#
# - Particle tree
#   - Branch n instead of nPart, name cannot be chosen in uproot
#   - Branch ImpactParameter renamed to impact
#   - Branch E is redundant, we skip this to save space
#   - Extra branches: parent
#
# For chromo in default configuration, the vertex locations are not interesting,
# so we don't write them. Long-lived particles are final state, and there is no
# interesting information in the vertices of very short-lived particles.
class Root(Writer):
    def __init__(self, file, model, write_vertices=False, buffer_size=100000):
        try:
            import uproot
        except ModuleNotFoundError:
            _raise_import_error("uproot", "write ROOT files")

        assert GeV == 1
        assert millibarn == 1

        kin = model.kinematics
        header = {
            "model": model.label,
            "seed": model.seed,
            "projectile_id": int(kin.p1),
            "projectile_momentum": kin.beams[0][2],
            "target_id": (
                repr(kin.p2) if isinstance(kin.p2, CompositeTarget) else int(kin.p2)
            ),
            "target_momentum": kin.beams[1][2],
        }
        header.update(
            {
                f"sigma_{k}": v
                for (k, v) in dataclasses.asdict(model.cross_section()).items()
                if not np.isnan(v)
            }
        )
        header.update(
            {
                "energy_unit": "GeV",
                "sigma_unit": "mb",
            }
        )

        self._file = uproot.recreate(file)

        self._event_buffers = {
            "impact": np.empty(buffer_size, FLOAT_TYPE),
        }
        self._particle_buffers = {
            "px": np.empty(buffer_size, FLOAT_TYPE),
            "py": np.empty(buffer_size, FLOAT_TYPE),
            "pz": np.empty(buffer_size, FLOAT_TYPE),
            "m": np.empty(buffer_size, FLOAT_TYPE),
            "pdgid": np.empty(buffer_size, INT_TYPE),
            "status": np.empty(buffer_size, INT_TYPE),
            "parent": np.empty(buffer_size, INT_TYPE),
        }
        if write_vertices:
            header["length_unit"] = "mm"
            self._particle_buffers.update(
                {
                    "vx": np.empty(buffer_size, FLOAT_TYPE),
                    "vy": np.empty(buffer_size, FLOAT_TYPE),
                    "vz": np.empty(buffer_size, FLOAT_TYPE),
                    "vt": np.empty(buffer_size, FLOAT_TYPE),
                }
            )

        self._header = "\n" + "\n".join(f"{k}: {v}" for (k, v) in header.items())
        self._tree = None
        self._lengths = []
        self._iparticle = 0

    def __exit__(self, *args):
        if self._iparticle > 0:
            self._write_buffers()
        return self._file.__exit__(*args)

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
            # uproot currently does not provide a high-level way
            # to write jagged arrays and set a title. This is really bad.
            # Until this is implemented, we force writing of a title
            # by monkeypatching WritableDirectory.mktree. Obviously, this
            # is a very bad and brittle solution, but you gotta do what
            # you gotta do.

            from uproot.writing import WritableDirectory

            mktree = WritableDirectory.mktree
            title = self._header

            def modded_mktree(self, name, metadata):
                return mktree(self, name, metadata, title=title)

            WritableDirectory.mktree = modded_mktree
            self._file["event"] = chunk
            WritableDirectory.mktree = mktree
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

        event_size = len(event)
        buffer_size = len(self._particle_buffers["px"])
        if self._iparticle + event_size > buffer_size:
            self._write_buffers()
            if event_size > buffer_size:
                # this should never happen
                raise RuntimeError("event is larger than buffer")

        ievent = len(self._lengths)
        self._event_buffers["impact"][ievent] = getattr(event, "impact_parameter", 0.0)

        a = self._iparticle
        self._iparticle += event_size
        b = self._iparticle
        self._lengths.append(b - a)
        for key, val in self._particle_buffers.items():
            if key == "parent":
                val[a:b] = event.mothers[:, 0]
            elif key == "pdgid":
                val[a:b] = event.pid
            else:
                val[a:b] = getattr(event, key)


class Svg(Writer):
    def __init__(self, file, model):
        self._idx = 0
        self._template = (file.parent, file.stem, file.suffix)

    def write(self, event):
        try:
            ge = event.to_hepmc3()
        except ModuleNotFoundError:
            _raise_import_error("pyhepmc", "write SVGs")
        if not hasattr(ge, "_repr_html_"):
            _raise_import_error("graphviz", "write SVGs")
        svg = ge._repr_html_()
        odir, name, ext = self._template
        fn = Path(odir) / f"{name}_{self._idx:03}{ext}"
        self._idx += 1
        # explicit utf-8 encoding is required on Windows
        with open(fn, "w", encoding="utf-8") as f:
            f.write(svg)


class Hepmc(Writer):
    def __init__(self, file, model):
        try:
            from pyhepmc._core import pyiostream
            from pyhepmc.io import _WrappedWriter, WriterAscii
        except ModuleNotFoundError:
            _raise_import_error("pyhepmc", "write HepMC files")

        import gzip

        # TODO add metadata to GenRunInfo, needs
        # fix in pyhepmc

        # TODO fix the following in pyhepmc, we should be able
        # to use the public API and not these secrets
        op = gzip.open if file.suffix == ".gz" else open
        self._file = op(file, "wb")
        self._ios = pyiostream(self._file)
        self._writer = _WrappedWriter(self._ios, None, WriterAscii)

    def __enter__(self):
        return self._writer

    def __exit__(self, *args):
        self._writer.__exit__(*args)
        self._ios.__exit__(*args)
        self._file.__exit__(*args)

    def write(self, event):
        self._writer.write(event)


def lhe(file):
    raise SystemExit("LHE not yet supported")
