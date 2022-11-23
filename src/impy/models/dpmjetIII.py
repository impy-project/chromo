from impy.common import MCRun, MCEvent, CrossSectionData
from impy import impy_config
from impy.util import info, _cached_data_dir


class DpmjetIIIEvent(MCEvent):
    """Wrapper class around DPMJET-III HEPEVT-style particle stack."""

    _hepevt = "dtevt1"
    _phep = "phkk"
    _vhep = "vhkk"
    _nevhep = "nevhkk"
    _nhep = "nhkk"
    _idhep = "idhkk"
    _isthep = "isthkk"
    _jmohep = "jmohkk"
    _jdahep = "jdahkk"

    def _charge_init(self, npart):
        return self._lib.dtpart.iich[self._lib.dtevt2.idbam[:npart] - 1]

    def _get_impact_parameter(self):
        return self._lib.dtglcp.bimpac

    def _get_n_wounded(self):
        return self._lib.dtglcp.nwasam, self._lib.dtglcp.nwbsam

    # Unfortunately not that simple since this is bounced through
    # entire code as argument not in COMMON
    # @property
    # def n_inel_NN_interactions(self):
    #     """Number of inelastic nucleon-nucleon interactions"""
    #     return self._lib.dtglcp.nwtsum


# =========================================================================
# DpmjetIIIMCRun
# =========================================================================
class DpmjetIIIRun(MCRun):
    """Implements all abstract attributes of MCRun for the
    DPMJET-III series of event generators.

    It should work identically for the new 'dpmjet3' module and the legacy
    dpmjet306. No special constructor is necessary and everything is
    handled by the default constructor of the base class.
    """

    _name = "DPMJET-III"
    _event_class = DpmjetIIIEvent
    _output_frame = "center-of-mass"
    _param_file_name = "dpmjpar.dat"
    _evap_file_name = "dpmjet.dat"
    _data_url = (
        "https://github.com/impy-project/impy"
        + "/releases/download/zipped_data_v1.0/dpm3191_v001.zip"
    )

    def __init__(self, event_kinematics, seed=None, logfname=None):
        from impy.util import fortran_chars
        from impy.constants import sec2cm

        super().__init__(seed, logfname)

        self._lib.init_rmmard(self._seed)

        # Save maximal mass that has been inisialized
        # (DPMJET sometimes crashes if higher mass requested than initialized)
        self._max_A1 = event_kinematics.A1
        self._max_A2 = event_kinematics.A2

        self.event_kinematics = event_kinematics

        info(1, "Initializing DPMJET-III")

        data_dir = _cached_data_dir(self._data_url)

        # Set the dpmjpar.dat file
        if hasattr(self._lib, "pomdls") and hasattr(self._lib.pomdls, "parfn"):
            pfile = data_dir + self._param_file_name
            info(3, "DPMJET parameter file at", pfile)
            self._lib.pomdls.parfn = fortran_chars(self._lib.pomdls.parfn, pfile)

        # Set the data directory for the other files
        if hasattr(self._lib, "poinou") and hasattr(self._lib.poinou, "datdir"):
            pfile = data_dir
            info(3, "DPMJET data dir is at", pfile)
            self._lib.poinou.datdir = fortran_chars(self._lib.poinou.datdir, pfile)
            self._lib.poinou.lendir = len(pfile)

        if hasattr(self._lib, "dtimpy"):
            evap_file = data_dir + self._evap_file_name
            info(3, "DPMJET evap file at", evap_file)
            self._lib.dtimpy.fnevap = fortran_chars(self._lib.dtimpy.fnevap, evap_file)

        k = self.event_kinematics
        self._lib.dt_init(-1, k.plab, k.A1, k.Z1, k.A2, k.Z2, k.p1pdg, iglau=0)

        if impy_config["user_frame"] == "center-of-mass":
            self._lib.dtflg1.iframe = 2
            self._output_frame = "center-of-mass"
        elif impy_config["user_frame"] == "laboratory":
            self._lib.dtflg1.iframe = 1
            self._output_frame = "laboratory"

        # Relax momentum and energy conservation checks at very high energies
        if k.ecm > 5e4:
            # Relative allowed deviation
            self._lib.pomdls.parmdl[76] = 0.045
            # Absolute allowed deviation
            self._lib.pomdls.parmdl[77] = 0.395
            # Relax threshhold of rejected events for variable energy runs
            self._lib.pomdls.ipamdl[178] = 5000

        # Relax momentum and energy conservation checks at very high energies
        if k.ecm > 5e4:
            # Relative allowed deviation
            self._lib.pomdls.parmdl[76] = 0.05
            # Absolute allowed deviation
            self._lib.pomdls.parmdl[77] = 0.05

        self._set_final_state_particles()
        # Prevent DPMJET from overwriting decay settings
        # self._lib.dtfrpa.ovwtdc = False
        # Set PYTHIA decay flags to follow all changes to MDCY
        self._lib.pydat1.mstj[21 - 1] = 1
        self._lib.pydat1.mstj[22 - 1] = 2
        # # Set ctau threshold in PYTHIA for the default stable list
        self._lib.pydat1.parj[70] = impy_config["tau_stable"] * sec2cm * 10.0  # mm

    def _cross_section(self, evt_kin=None, precision=None):
        # we override to set precision
        if evt_kin is None:
            k = self.event_kinematics
        else:
            k = evt_kin
        info(10, "Cross section for", k.A1, k.A2, self._lib.idt_icihad(k.p1pdg))
        if precision is None:
            self._lib.dt_xsglau(
                k.A1, k.A2, self._lib.idt_icihad(k.p1pdg), 0, 0, k.ecm, 1, 1, 1
            )
        else:
            saved = self._lib.dtglgp.jstatb
            # Set number of trials for Glauber model integration
            self._lib.dtglgp.jstatb = precision
            self._lib.dt_xsglau(
                k.A1, k.A2, self._lib.idt_icihad(k.p1pdg), 0, 0, k.ecm, 1, 1, 1
            )
            self._lib.dtglgp.jstatb = saved

        return CrossSectionData._from_inel(self._lib.dtglxs.xspro[0, 0, 0])

    def _dpmjet_tup(self):
        """Constructs an tuple of arguments for calls to event generator
        from given event kinematics object."""
        k = self.event_kinematics
        info(
            20,
            "Request DPMJET ARGs tuple:\n",
            (k.A1, k.Z1, k.A2, k.Z2, self._lib.idt_icihad(k.p1pdg), k.elab),
        )
        return (k.A1, k.Z1, k.A2, k.Z2, self._lib.idt_icihad(k.p1pdg), k.elab)

    def _set_event_kinematics(self, k):
        # nothing to be done here except input validation, since
        # initialization and event generation is done in _generate_event
        info(5, "Setting event kinematics")
        if k.A1 > self._max_A1 or k.A2 > self._max_A2:
            raise ValueError(
                "Maximal initialization mass exceeded "
                f"{k.A1}/{self._max_A1}, {k.A2}/{self._max_A2}"
            )

        # AF: No idea yet, but apparently this functionality was around?!
        # if hasattr(k, 'beam') and hasattr(self._lib, 'init'):
        #     print(self.class_name + "::_set_event_kinematics():" +
        #           "setting beam params", k.beam)
        #     self._lib.dt_setbm(k.A1, k.Z1, k.A2, k.Z2, k.beam[0], k.beam[1])
        #     print 'OK'

    def _attach_log(self, fname=None):
        """Routes the output to a file or the stdout."""
        fname = impy_config["output_log"] if fname is None else fname
        if fname == "stdout":
            lun = 6
            info(5, "Output is routed to stdout.")
        else:
            lun = self._attach_fortran_logfile(fname)
            info(5, "Output is routed to", fname, "via LUN", lun)

        if hasattr(self._lib, "dtflka"):
            self._lib.dtflka.lout = lun
            self._lib.dtflka.lpri = 50
        elif hasattr(self._lib, "dtiont"):
            self._lib.dtiont.lout = lun
        else:
            raise Exception("Unknown DPMJET version, IO common block not detected.")

        self._lib.pydat1.mstu[10] = lun

    def _set_stable(self, pdgid, stable):
        if abs(pdgid) == 2212:
            return
        kc = self._lib.pycomp(pdgid)
        self._lib.pydat3.mdcy[kc - 1, 0] = 0 if stable else 1

    def _generate_event(self):
        reject = self._lib.dt_kkinc(*self._dpmjet_tup(), kkmat=-1)
        self._lib.dtevno.nevent += 1
        return not reject


class DpmjetIII191(DpmjetIIIRun):
    _version = "19.1"
    _library_name = "_dpmjetIII191"


class DpmjetIII193(DpmjetIIIRun):
    _version = "19.3"
    _library_name = "_dpmjetIII193"


class DpmjetIII306(DpmjetIIIRun):
    _version = "3.0-6"
    _library_name = "_dpmjet306"
    _param_file_name = "fitpar.dat"
    _data_url = (
        "https://github.com/impy-project/impy"
        + "/releases/download/zipped_data_v1.0/dpm3_v001.zip"
    )


class DpmjetIII193_DEV(DpmjetIIIRun):
    _version = "19.3-dev"
    _library_name = "_dev_dpmjetIII193"
