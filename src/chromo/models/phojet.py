from chromo.common import Model, MCEvent, CrossSectionData
from chromo.util import fortran_chars, _cached_data_dir
from chromo.kinematics import EventFrame
from chromo.constants import standard_projectiles
from particle import literals as lp


class PhojetEvent(MCEvent):
    """Wrapper class around (pure) PHOJET particle stack."""

    _hepevt = "poevt1"

    def _charge_init(self, npart):
        return self._lib.poevt2.icolor[0, :npart] / 3

    def _gen_cut_info(self):
        """Init variables tracking the number of soft and hard cuts"""

        self._mpi = self._lib.podebg.kspom + self._lib.podebg.khpom
        self._kspom = self._lib.podebg.kspom
        self._khpom = self._lib.podebg.khpom
        self._ksoft = self._lib.podebg.ksoft
        self._khard = self._lib.podebg.khard

    @property
    def mpi(self):
        """Total number of cuts"""
        return self._lib.podebg.kspom + self._lib.podebg.khpom

    @property
    def kspom(self):
        """Total number of soft cuts"""
        return self._lib.podebg.kspom

    @property
    def khpom(self):
        """Total number of hard cuts"""
        return self._lib.podebg.khpom

    @property
    def ksoft(self):
        """Total number of realized soft cuts"""
        return self._lib.podebg.ksoft

    @property
    def khard(self):
        """Total number of realized hard cuts"""
        return self._lib.podebg.khard

    # def elastic_t(self):
    #     """Squared momentum transfer t for elastic interaction.

    #     This only makes sense if the interaction is indeed elastic
    #     and only 4 particles are on the stack. The initial 2 and
    #     the final 2. Handle with care!!
    #     """
    #     return np.sum(
    #         np.square(self._lib.poevt1.phep[0:4, 0] -
    #                   self._lib.poevt1.phep[0:4, 5]))


class PHOJETRun(Model):
    """Implements all abstract attributes of MCRun for the
    PHOJET series of event generators.

    PHOJET is part of DPMJET and is run via the same fortran library.
    The results for pp and hadron-p `should` be the same.
    """

    _name = "PhoJet"
    _event_class = PhojetEvent
    _frame = None
    _projectiles = standard_projectiles
    _targets = {lp.proton.pdgid, lp.neutron.pdgid}
    _param_file_name = "dpmjpar.dat"
    _data_url = (
        "https://github.com/impy-project/chromo"
        + "/releases/download/zipped_data_v1.0/dpm3191_v001.zip"
    )

    def _once(self):
        import chromo

        data_dir = _cached_data_dir(self._data_url)
        # Set the dpmjpar.dat file
        if hasattr(self._lib, "pomdls") and hasattr(self._lib.pomdls, "parfn"):
            pfile = data_dir + self._param_file_name
            self._lib.pomdls.parfn = fortran_chars(self._lib.pomdls.parfn, pfile)

        # Set the data directory for the other files
        if hasattr(self._lib, "poinou") and hasattr(self._lib.poinou, "datdir"):
            pfile = data_dir
            self._lib.poinou.datdir = fortran_chars(self._lib.poinou.datdir, pfile)
            self._lib.poinou.lendir = len(pfile)

        # Set debug level of the generator
        for i in range(self._lib.podebg.ideb.size):
            self._lib.podebg.ideb[i] = chromo.debug_level

        # Setup logging
        lun = 6  # stdout
        verbosity = 5 if chromo.debug_level else 1
        # Standalone phojet IO block
        if hasattr(self._lib, "poinou"):
            print(1)
            self._lib.poinou.lo = lun
            self._lib.poinou.lpri = verbosity
        if hasattr(self._lib, "dtflka"):
            self._lib.dtflka.lout = lun
            self._lib.dtflka.lpri = verbosity
        elif hasattr(self._lib, "dtiont"):
            self._lib.dtiont.lout = lun
        else:
            assert (
                False
            ), "Unknown PHOJET (DPMJET) version, IO common block not detected"
        self._lib.pydat1.mstu[10] = lun

        # Initialize PHOJET's parameters
        # If PHOJET is run through the DPMJET library, this init needs -2 else -1
        if self._lib.pho_init(-1, lun, verbosity):
            raise RuntimeError("unable to initialize or set LUN")

        process_switch = self._lib.poprcs.ipron
        # non-diffractive elastic scattering (1 - on, 0 - off)
        process_switch[0, 0] = 1
        # elastic scattering
        process_switch[1, 0] = 0
        # quasi-elastic scattering (for incoming photons only)
        process_switch[2, 0] = 1
        # central diffration (double-pomeron scattering)
        process_switch[3, 0] = 1
        # particle 1 single diff. dissociation
        process_switch[4, 0] = 1
        # particle 2 single diff. dissociation
        process_switch[5, 0] = 1
        # double diff. dissociation
        process_switch[6, 0] = 1
        # direct photon interaction (for incoming photons only)
        process_switch[7, 0] = 1

        # Set PYTHIA decay flags to follow all changes to MDCY
        self._lib.pydat1.mstj[21 - 1] = 1
        self._lib.pydat1.mstj[22 - 1] = 2

    def _cross_section(self, kin):
        self._set_kinematics(kin)

        if hasattr(self._lib, "pho_setpcomb"):
            # The old 1.12 version of phojet doesn't have this function
            self._lib.pho_setpcomb()

        self._lib.pho_xsect(1, 0, kin.ecm)

        return CrossSectionData(
            total=float(self._lib.pocsec.sigtot),
            elastic=float(self._lib.pocsec.sigela),
            inelastic=float(self._lib.pocsec.sigine),
            diffractive_xb=float(
                self._lib.pocsec.siglsd[0] + self._lib.pocsec.sighsd[0]
            ),
            diffractive_ax=float(
                self._lib.pocsec.siglsd[1] + self._lib.pocsec.sighsd[1]
            ),
            diffractive_xx=float(self._lib.pocsec.sigldd + self._lib.pocsec.sighdd),
            diffractive_axb=float(self._lib.pocsec.sigcdf[0]),
        )

    def _set_stable(self, pdgid, stable):
        kc = self._lib.pycomp(pdgid)
        self._lib.pydat3.mdcy[kc - 1, 0] = not stable

    def _set_kinematics(self, k):
        if k.frame == EventFrame.FIXED_TARGET:
            self._lib.dtflg1.iframe = 1
            self._frame = EventFrame.FIXED_TARGET
        else:
            self._lib.dtflg1.iframe = 2
            self._frame = EventFrame.CENTER_OF_MASS
        self._lib.pho_setpar(1, k.p1, 0, 0.0)
        self._lib.pho_setpar(2, k.p2, 0, 0.0)
        self._beams = k.beams

        # Initialize kinematics and tables (only once needed)
        # LIMITATION: Once this has been called for a given set of particles,
        # it cannot be called again for different values of k.p1, k.p2. This
        # means we cannot change particles in Phojet.
        if self._lib.pho_event(-1, *self._beams)[1]:
            raise RuntimeError(
                "initialization failed with the current event kinematics"
            )

    def _generate(self):
        return not self._lib.pho_event(1, *self._beams)[1]


class Phojet112(PHOJETRun):
    _version = "1.12-35"
    _library_name = "_phojet112"
    _param_file_name = "fitpar.dat"
    _projectiles = {
        lp.photon.pdgid,
        lp.proton.pdgid,
    }  # FIXME: should allow photons and hadrons
    _data_url = (
        "https://github.com/impy-project/chromo"
        + "/releases/download/zipped_data_v1.0/dpm3_v001.zip"
    )


class Phojet191(PHOJETRun):
    _version = "19.1"
    _library_name = "_phojet191"


class Phojet193(PHOJETRun):
    _version = "19.3"
    _library_name = "_phojet193"
