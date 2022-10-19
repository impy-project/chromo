import numpy as np
from impy.common import MCRun, MCEvent, RMMARDState, impy_config
from impy.util import info
import dataclasses


class SibyllEvent(MCEvent):
    """Wrapper class around SIBYLL 2.1 & 2.3 particle stack."""

    _jdahep = None  # no child info

    def _charge_init(self, npart):
        return self._lib.schg.ichg[:npart]

    # Nuclear collision parameters
    @property
    def impact_parameter(self):
        """Impact parameter for nuclear collisions."""
        return self._lib.cnucms.b

    @property
    def n_wounded_A(self):
        """Number of wounded nucleons side A"""
        return self._lib.cnucms.na

    @property
    def n_wounded_B(self):
        """Number of wounded nucleons side B"""
        return self._lib.cnucms.nb

    @property
    def n_NN_interactions(self):
        """Number of inelastic nucleon-nucleon interactions"""
        return self._lib.cnucms.ni


@dataclasses.dataclass
class RMMARDSib(RMMARDState):
    _gasdev_iset: np.ndarray = None

    def _record_state(self, generator):
        super()._record_state(generator)
        self._gasdev_iset = generator._lib.rndmgas.iset
        return self

    def _restore_state(self, generator):
        super()._restore_state(generator)
        generator._lib.rndmgas.iset = self._gasdev_iset
        return self

    def __eq__(self, other: object) -> bool:
        return super().__eq__(other) and np.array_equal(
            self._gasdev_iset, other._gasdev_iset
        )


class SIBYLLRun(MCRun):
    """Implements all abstract attributes of MCRun for the
    SIBYLL 2.1, 2.3 and 2.3c event generators."""

    _name = "SIBYLL"
    _event_class = SibyllEvent
    _output_frame = "center-of-mass"

    def __init__(self, event_kinematics, seed=None, logfname=None):
        super().__init__(seed, logfname)

        self._lib.s_debug.ndebug = impy_config["sibyll"]["debug_level"]

        self._lib.sibini(self._seed)
        # Set the internal state of GASDEV function (rng) to 0
        self._lib.rndmgas.iset = 0
        self._lib.pdg_ini()

        # This calls _set_event_kinematics which uses self._lib.isib_pdg2pid
        # which works only after an initialization call to self._lib.pdg_ini()
        self.event_kinematics = event_kinematics

        self._set_final_state_particles()

    def _sigma_inel(self, evt_kin):
        sigproj = None
        if abs(evt_kin.p1pdg) in [2212, 2112, 3112]:
            sigproj = 1
        elif abs(evt_kin.p1pdg) == 211:
            sigproj = 2
        elif abs(evt_kin.p1pdg) == 321:
            sigproj = 3
        else:
            raise ValueError(
                f"No cross section available for projectile {evt_kin.p1pdg}"
            )

        if evt_kin.p1_is_nucleus:
            raise ValueError(f"Nuclear projectiles not supported by {self.name}")

        if evt_kin.p2_is_nucleus:
            # Return production cross section for nuclear target
            try:
                return self._lib.sib_sigma_hnuc(sigproj, evt_kin.A2, evt_kin.ecm)[0]
            except AttributeError:
                raise ValueError(f"Nuclear projectiles not supported by {self.label}")

        return self._lib.sib_sigma_hp(sigproj, evt_kin.ecm)[2]

    def sigma_inel_air(self):
        """Inelastic cross section according to current
        event setup (energy, projectile, target)"""
        k = self.event_kinematics
        sigproj = None
        if abs(k.p1pdg) in [2212, 2112, 3112]:
            sigproj = 1
        elif abs(k.p1pdg) == 211:
            sigproj = 2
        elif abs(k.p1pdg) == 321:
            sigproj = 3
        else:
            info(0, "No cross section available for projectile", k.p1pdg)
            raise Exception("Input error")
        sigma = self._lib.sib_sigma_hair(sigproj, k.ecm)
        if isinstance(sigma, tuple):
            return sigma[0]
        return sigma

    def _set_event_kinematics(self, k):
        info(5, "Setting event kinematics")
        if k.p1_is_nucleus:
            raise ValueError(
                f"Projectile nuclei not natively supported in SIBYLL, A = {k.A1}"
            )
        if k.p2_is_nucleus and k.A2 > 20:
            raise ValueError(f"Target with A>20 not supported, A = {k.A2}")
        self._sibproj = self._lib.isib_pdg2pid(k.p1pdg)
        self._iatarg = k.A2

    def _attach_log(self, fname=None):
        """Routes the output to a file or the stdout."""
        fname = impy_config["output_log"] if fname is None else fname
        if fname == "stdout":
            self._lib.s_debug.lun = 6
            info(5, "Output is routed to stdout.")
        else:
            lun = self._attach_fortran_logfile(fname)
            self._lib.s_debug.lun = lun
            info(5, "Output is routed to", fname, "via LUN", lun)

    def _set_stable(self, pdgid, stable):
        sid = abs(self._lib.isib_pdg2pid(pdgid))
        if abs(pdgid) == 311:
            info(1, "Ignores K0. Using K0L/S instead")
            self.set_stable(130, stable)
            self.set_stable(310, stable)
            return
        idb = self._lib.s_csydec.idb
        if sid == 0 or sid > idb.size - 1:
            return
        if stable:
            idb[sid - 1] = -np.abs(idb[sid - 1])
        else:
            idb[sid - 1] = np.abs(idb[sid - 1])

    def _generate_event(self):
        k = self.event_kinematics
        self._lib.sibyll(self._sibproj, self._iatarg, k.ecm)
        self._lib.decsib()
        self._lib.sibhep3()
        return 0  # SIBYLL never rejects

    @property
    def random_state(self):
        return RMMARDSib()._record_state(self)

    @random_state.setter
    def random_state(self, rng_state):
        rng_state._restore_state(self)


class Sibyll21(SIBYLLRun):
    _version = "2.1"
    _library_name = "_sib21"

    def _generate_event(self):
        k = self.event_kinematics
        self._lib.sibyll(self._sibproj, self._iatarg, k.ecm)
        self._lib.decsib()
        self._lib.sibhep1()
        return 0  # SIBYLL never rejects


class Sibyll23(SIBYLLRun):
    _version = "2.3"
    _library_name = "_sib23"


class Sibyll23c(SIBYLLRun):
    _version = "2.3c"
    _library_name = "_sib23c01"


# undocumented patch version
class Sibyll23c00(SIBYLLRun):
    _version = "2.3c00"
    _library_name = "_sib23c00"


# identical to 2.3c
class Sibyll23c01(SIBYLLRun):
    _version = "2.3c01"
    _library_name = "_sib23c01"


# undocumented patch version
class Sibyll23c02(SIBYLLRun):
    _version = "2.3c02"
    _library_name = "_sib23c02"


# The c03 version was also in CORSIKA until 2020
class Sibyll23c03(SIBYLLRun):
    _version = "2.3c03"
    _library_name = "_sib23c03"


# The latest patch c04 was renamed to d, to generate less confusion
class Sibyll23d(SIBYLLRun):
    _version = "2.3d"
    _library_name = "_sib23d"
