'''
Created on 14.04.2014

@author: afedynitch
'''

import numpy as np
from impy.common import MCRun, MCEvent, impy_config, pdata
from impy.util import standard_particles, info


class DpmjetIIIEvent(MCEvent):
    """Wrapper class around DPMJET-III particle stack."""

    def __init__(self, lib, event_kinematics, frame):
        # HEPEVT (style) common block
        evt = lib.dtevt1
        # Number of entries on stack
        npart = evt.nhkk

        # Filter stack for charged particles if selected
        if impy_config["event_scope"] == 'charged':
            sel = np.where((evt.isthkk[:npart] == 1) & (
                lib.dtpart.iich[lib.dtevt2.idbam[:npart] - 1] != 0))
        elif impy_config["event_scope"] == 'stable':
            sel = np.where(evt.isthkk[:npart] == 1)
        else:
            raise Exception("not implemented, yet")

        # Save selector for implementation of on-demand properties
        self.sel = sel

        MCEvent.__init__(
            self,
            event_kinematics=event_kinematics,
            lib=lib,
            px=evt.phkk[0, sel][0],
            py=evt.phkk[1, sel][0],
            pz=evt.phkk[2, sel][0],
            en=evt.phkk[3, sel][0],
            p_ids=evt.idhkk[sel],
            npart=npart,
            frame=frame)

    @property
    def charge(self):
        return self.lib.dtpart.iich[self.lib.dtevt2.idbam[self.sel] - 1]

    @property
    def mass(self):
        return self.lib.dtevt1.phkk[4, self.sel][0]

    # Nuclear collision parameters

    @property
    def impact_parameter(self):
        """Returns impact parameter for nuclear collisions."""
        return self.lib.dtglcp.bimpac

    @property
    def n_wounded_A(self):
        """Number of wounded nucleons side A"""
        return self.lib.dtglcp.nwasam
    
    @property
    def n_wounded_B(self):
        """Number of wounded nucleons side B"""
        return self.lib.dtglcp.nwbsam

    # Unfortunately not that simple since this is bounced through
    # entire code as argument not in COMMON
    # @property
    # def n_inel_NN_interactions(self):
    #     """Number of inelastic nucleon-nucleon interactions"""
    #     return self.lib.dtglcp.nwtsum


#=========================================================================
# DpmjetIIIMCRun
#=========================================================================
class DpmjetIIIRun(MCRun):
    """Implements all abstract attributes of MCRun for the 
    DPMJET-III series of event generators.
    
    It should work identically for the new 'dpmjet3' module and the legacy
    dpmjet306.
    """

    def __init__(self, libref, event_class=None, **kwargs):

        if event_class is None:
            self._event_class = DpmjetIIIEvent
        else:
            self._event_class = event_class

        self._frame = 'center-of-mass'

        MCRun.__init__(self, libref, **kwargs)
    
    @property
    def frame(self):
        return self._frame

    @property
    def name(self):
        """Event generator name"""
        return "DPMJET"

    @property
    def version(self):
        """Event generator version"""
        # Needs some sort of smart handling here
        return self.lib.dtimpy.version

    def sigma_inel(self):
        """Inelastic cross section according to current
        event setup (energy, projectile, target)"""
        k = self._curr_event_kin
        info(10, 'Cross section for', k.A1, k.A2, self.lib.idt_icihad(k.p1pdg))
        self.lib.dt_xsglau(k.A1, k.A2, self.lib.idt_icihad(k.p1pdg), 0, 0,
                           k.ecm, 1, 1, 1)
        return self.lib.dtglxs.xspro[0, 0, 0]

    def _dpmjet_tup(self):
        """Constructs an tuple of arguments for calls to event generator
        from given event kinematics object."""
        k = self._curr_event_kin
        info(20, 'Request DPMJET ARGs tuple:\n', 
            (k.A1, k.Z1, k.A2, k.Z2, self.lib.idt_icihad(k.p1pdg), k.elab))
        return (k.A1, k.Z1, k.A2, k.Z2, self.lib.idt_icihad(k.p1pdg), k.elab)

    def set_event_kinematics(self, event_kinematics):
        """Set new combination of energy, momentum, projectile
        and target combination for next event."""
        info(5, 'Setting event kinematics')
        self._curr_event_kin = event_kinematics

        # AF: No idea yet, but apparently this functionality was around?!
        # if hasattr(k, 'beam') and hasattr(self.lib, 'init'):
        #     print(self.class_name + "::set_event_kinematics():" +
        #           "setting beam params", k.beam)
        #     self.lib.dt_setbm(k.A1, k.Z1, k.A2, k.Z2, k.beam[0], k.beam[1])
        #     print 'OK'

    def attach_log(self):
        """Routes the output to a file or the stdout."""
        fname = impy_config['output_log']
        if fname == 'stdout':
            lun = 6
            info(5, 'Output is routed to stdout.')
        else:
            lun = self._attach_fortran_logfile(fname)
            info(5, 'Output is routed to', fname, 'via LUN', lun)

        if hasattr(self.lib, 'dtflka'):
            self.lib.dtflka.lout = 6
            self.lib.dtflka.lpri = 50
        elif hasattr(self.lib, 'dtiont'):
            self.lib.dtiont.lout = lun
        else:
            raise Exception(
                'Unknown DPMJET version, IO common block not detected.')

        self.lib.pydat1.mstu[10] = lun

    def init_generator(self, event_kinematics):
        from impy.util import clear_and_set_fortran_chars

        self._abort_if_already_initialized()

        # Comprise DPMJET input from the event kinematics object
        self.set_event_kinematics(event_kinematics)
        k = self._curr_event_kin

        dpm_conf = impy_config['dpmjetIII']
        info(1, 'First initialization')
        # print self.class_name + "::init_generator(): calling init with:", \
        #     - 1, k.elab, PM, PCH, TM, TCH, k.p1pdg

        # Set the dpmjpar.dat file
        if hasattr(self.lib, 'pomdls.parfn'):
            pfile = dpm_conf['param_file']['III-2018.1']
            info(10, 'Using PHOJET parameter file at', pfile)
            clear_and_set_fortran_chars(self.lib.pomdls.parfn, pfile)
        
        if hasattr(self.lib, 'dtimpy'):
            evap_file = dpm_conf['evap_file']['3.0-6']
            info(10, 'Using DPMJET evap file at', evap_file)
            clear_and_set_fortran_chars(self.lib.dtimpy.fnevap, evap_file)

        self.lib.dt_init(
            -1, dpm_conf['e_max'], k.A1, k.Z1, k.A2, k.Z2, k.p1pdg, iglau=0)
        # Put protection to not run this stuff again
        self.lib.init = True

        if impy_config['user_frame'] == 'center-of-mass':
            self.lib.dtflg1.iframe = 2
            self._frame = 'center-of-mass'
        elif impy_config['user_frame'] == 'laboratory':
            self.lib.dtflg1.iframe = 1
            self._frame = 'laboratory'

        # if self.def_settings:
        #     print self.class_name + "::init_generator(): Using default settings:", \
        #         self.def_settings.__class__.__name__
        #     self.def_settings.enable()

        self._define_default_fs_particles()

    def set_stable(self, pdgid):
        kc = self.lib.pycomp(pdgid)
        self.lib.pydat3.mdcy[kc - 1, 0] = 0
        info(5, 'defining', pdgid, 'as stable particle')

    def generate_event(self):
        reject = self.lib.dt_kkinc(*self._dpmjet_tup(), kkmat=-1)
        self.lib.dtevno.nevent += 1
        return reject

