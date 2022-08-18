"""
Created on 17.03.2014

@author: afedynitch
"""

import numpy as np
from impy.common import MCRun, MCEvent, impy_config
from impy.util import info


class SibyllEvent(MCEvent):
    """Wrapper class around SIBYLL 2.1 & 2.3 particle stack."""

    # Workaround for no data on vertext positions in SIBYLL
    _no_vertex_data = None

    def _charge_init(self, npart):
        return self._lib.schg.ichg[:npart]

    # Nuclear collision parameters
    @property
    def impact_parameter(self):
        """Impact parameter for nuclear collisions."""
        return self.lib.cnucms.b

    @property
    def n_wounded_A(self):
        """Number of wounded nucleons side A"""
        return self.lib.cnucms.na

    @property
    def n_wounded_B(self):
        """Number of wounded nucleons side B"""
        return self.lib.cnucms.nb

    @property
    def n_NN_interactions(self):
        """Number of inelastic nucleon-nucleon interactions"""
        return self.lib.cnucms.ni


class SIBYLLRun(MCRun):
    """Implements all abstract attributes of MCRun for the
    SIBYLL 2.1, 2.3 and 2.3c event generators."""

    def sigma_inel(self, *args, **kwargs):
        """Inelastic cross section according to current
        event setup (energy, projectile, target)"""
        k = self._curr_event_kin
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

        if k.p1_is_nucleus:
            raise Exception("Nuclear projectiles not supported by SIBYLL.")

        if k.p2_is_nucleus:
            # Return production cross section for nuclear target
            try:
                return self.lib.sib_sigma_hnuc(sigproj, k.A2, self._ecm)[0]
            except AttributeError:
                return "Nuclear cross section not supported for this SIBYLL version"

        return self.lib.sib_sigma_hp(sigproj, self._ecm)[2]

    def sigma_inel_air(self):
        """Inelastic cross section according to current
        event setup (energy, projectile, target)"""
        k = self._curr_event_kin
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
        sigma = self.lib.sib_sigma_hair(sigproj, self._ecm)
        if not isinstance(sigma, tuple):
            return sigma
        else:
            return sigma[0]

    def _set_event_kinematics(self, event_kinematics):
        """Set new combination of energy, momentum, projectile
        and target combination for next event."""

        info(5, "Setting event kinematics.")
        info(10, event_kinematics)
        k = event_kinematics
        if k.p1_is_nucleus:
            raise Exception("Projectile nuclei not natively supported in SIBYLL")
        elif k.p2_is_nucleus and k.A2 > 20:
            print(k.p2_is_nucleus, k.A2)
            raise Exception("Target nuclei with A>20 not supported in SIBYLL")
        self._sibproj = self.lib.isib_pdg2pid(k.p1pdg)
        self._iatarg = k.A2
        self._ecm = k.ecm
        self._curr_event_kin = event_kinematics

    def attach_log(self, fname=None):
        """Routes the output to a file or the stdout."""
        fname = impy_config["output_log"] if fname is None else fname
        if fname == "stdout":
            self.lib.s_debug.lun = 6
            info(5, "Output is routed to stdout.")
        else:
            lun = self._attach_fortran_logfile(fname)
            self.lib.s_debug.lun = lun
            info(5, "Output is routed to", fname, "via LUN", lun)

    def init_generator(self, event_kinematics, seed="random", logfname=None):
        from random import randint

        self._abort_if_already_initialized()

        if seed == "random":
            seed = randint(1000000, 10000000)
        else:
            seed = int(seed)
        info(5, "Using seed:", seed)

        self.lib.s_debug.ndebug = impy_config["sibyll"]["debug_level"]

        self.attach_log(fname=logfname)
        self.lib.sibini(int(seed))
        self.lib.pdg_ini()
        self.conv_hepevt = (
            self.lib.sibhep1 if "21" in self.lib.__name__ else self.lib.sibhep3
        )
        self._define_default_fs_particles()
        # _set_event_kinematics uses function self.lib.isib_pdg2pid
        # which works only after an initialization call to self.lib.pdg_ini()
        # therefore it should be called after this initalization
        self._set_event_kinematics(event_kinematics)

    def set_stable(self, pdgid, stable=True):
        sid = abs(self.lib.isib_pdg2pid(pdgid))
        if abs(pdgid) == 311:
            info(1, "Ignores K0. Use K0L/S 130/310 in final state definition.")
            return
        idb = self.lib.s_csydec.idb
        if sid == 0 or sid > idb.size - 1:
            return
        if stable:
            info(
                5, "defining as stable particle pdgid/sid = {0}/{1}".format(pdgid, sid)
            )
            idb[sid - 1] = -np.abs(idb[sid - 1])
        else:
            info(5, "pdgid/sid = {0}/{1} allowed to decay".format(pdgid, sid))
            idb[sid - 1] = np.abs(idb[sid - 1])

    def generate_event(self):
        self.lib.sibyll(self._sibproj, self._iatarg, self._ecm)
        self.lib.decsib()
        self.conv_hepevt()
        return 0  # SIBYLL never rejects


class Sibyll21(SIBYLLRun):
    def __init__(self, event_kinematics, seed="random", logfname=None):
        from impy.definitions import interaction_model_by_tag as models_dict

        interaction_model_def = models_dict["SIBYLL21"]
        super().__init__(interaction_model_def)
        self.init_generator(event_kinematics, seed, logfname)


class Sibyll23(SIBYLLRun):
    def __init__(self, event_kinematics, seed="random", logfname=None):
        from impy.definitions import interaction_model_by_tag as models_dict

        interaction_model_def = models_dict["SIBYLL23"]
        super().__init__(interaction_model_def)
        self.init_generator(event_kinematics, seed, logfname)


class Sibyll23c(SIBYLLRun):
    def __init__(self, event_kinematics, seed="random", logfname=None):
        from impy.definitions import interaction_model_by_tag as models_dict

        interaction_model_def = models_dict["SIBYLL23C"]
        super().__init__(interaction_model_def)
        self.init_generator(event_kinematics, seed, logfname)


class Sibyll23c00(SIBYLLRun):
    def __init__(self, event_kinematics, seed="random", logfname=None):
        from impy.definitions import interaction_model_by_tag as models_dict

        interaction_model_def = models_dict["SIBYLL23C00"]
        super().__init__(interaction_model_def)
        self.init_generator(event_kinematics, seed, logfname)


class Sibyll23c01(SIBYLLRun):
    def __init__(self, event_kinematics, seed="random", logfname=None):
        from impy.definitions import interaction_model_by_tag as models_dict

        interaction_model_def = models_dict["SIBYLL23C01"]
        super().__init__(interaction_model_def)
        self.init_generator(event_kinematics, seed, logfname)


class Sibyll23c02(SIBYLLRun):
    def __init__(self, event_kinematics, seed="random", logfname=None):
        from impy.definitions import interaction_model_by_tag as models_dict

        interaction_model_def = models_dict["SIBYLL23C02"]
        super().__init__(interaction_model_def)
        self.init_generator(event_kinematics, seed, logfname)


class Sibyll23c03(SIBYLLRun):
    def __init__(self, event_kinematics, seed="random", logfname=None):
        from impy.definitions import interaction_model_by_tag as models_dict

        interaction_model_def = models_dict["SIBYLL23C03"]
        super().__init__(interaction_model_def)
        self.init_generator(event_kinematics, seed, logfname)


class Sibyll23c04(SIBYLLRun):
    def __init__(self, event_kinematics, seed="random", logfname=None):
        from impy.definitions import interaction_model_by_tag as models_dict

        interaction_model_def = models_dict["SIBYLL23C04"]
        super().__init__(interaction_model_def)
        self.init_generator(event_kinematics, seed, logfname)


class Sibyll23d(SIBYLLRun):
    def __init__(self, event_kinematics, seed="random", logfname=None):
        from impy.definitions import interaction_model_by_tag as models_dict

        interaction_model_def = models_dict["SIBYLL23D"]
        super().__init__(interaction_model_def)
        self.init_generator(event_kinematics, seed, logfname)
