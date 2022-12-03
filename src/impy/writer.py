import numpy as np
from impy.constants import quarks_and_diquarks_and_gluons, GeV, millibarn
import dataclasses
from pathlib import Path

BUFFER_SIZE = 100000
INT_TYPE = np.int32
FLOAT_TYPE = np.float32


def _raise_import_error(name, task):
    raise ModuleNotFoundError(
        f"{name} not found, please install {name} (`pip install {name}`) to {task}"
    )


# Differences to CRMC
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
# For impy in default configuration, the vertex locations are not interesting,
# so we don't write them. Long-lived particles are final state, and there is no
# interesting information in the vertices of very short-lived particles.
class Root:
    def __init__(self, file, config, model, write_vertices=False):
        try:
            import uproot
        except ModuleNotFoundError:
            _raise_import_error("uproot", "write ROOT files")

        # FIXME projectile_momentum and target_momentum are 0 if -S option is used,
        # to be fixed in follow-up PR that fixes EventKinematics
        header = {
            "seed": model.seed,
            "projectile_id": config.projectile_id,
            "projectile_momentum": config.projectile_momentum / GeV,
            "target_id": config.target_id,
            "target_momentum": config.target_momentum / GeV,
            "model": model.label,
        }
        header.update(
            {
                f"sigma_{k}": v / millibarn
                for (k, v) in dataclasses.asdict(model.cross_section()).items()
                if not np.isnan(v)
            }
        )
        header.update(
            {
                "energy_unit": "GeV",
                "length_unit": "mm",
                "sigma_unit": "mb",
            }
        )

        self._header = "\n".join(f"{k}: {v}" for (k, v) in header.items())
        self._file = uproot.recreate(file)

        self._event_buffers = {
            "impact": np.empty(BUFFER_SIZE, FLOAT_TYPE),
        }
        self._particle_buffers = {
            "px": np.empty(BUFFER_SIZE, FLOAT_TYPE),
            "py": np.empty(BUFFER_SIZE, FLOAT_TYPE),
            "pz": np.empty(BUFFER_SIZE, FLOAT_TYPE),
            "m": np.empty(BUFFER_SIZE, FLOAT_TYPE),
            "pdgid": np.empty(BUFFER_SIZE, INT_TYPE),
            "status": np.empty(BUFFER_SIZE, INT_TYPE),
            "parent": np.empty(BUFFER_SIZE, INT_TYPE),
        }
        if write_vertices:
            self._particle_buffers.update(
                {
                    "vx": np.empty(BUFFER_SIZE, FLOAT_TYPE),
                    "vy": np.empty(BUFFER_SIZE, FLOAT_TYPE),
                    "vz": np.empty(BUFFER_SIZE, FLOAT_TYPE),
                    "vt": np.empty(BUFFER_SIZE, FLOAT_TYPE),
                }
            )

        self._tree = None
        self._lengths = []
        self._iparticle = 0

    def __enter__(self):
        return self

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
                if key in ("px", "py", "pz", "m"):
                    val[a:b] *= 1e3

        if b > BUFFER_SIZE // 2:
            self._write_buffers()


class Svg:
    def __init__(self, file, config, model):
        self._idx = 0
        self._template = (file.parent, file.stem, file.suffix)

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return

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
        with open(fn, "w") as f:
            f.write(svg)


class Hepmc:
    def __init__(self, file, config, model):
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
