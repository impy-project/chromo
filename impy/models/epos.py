'''
Created on 03.05.2016

@author: afedynitch
'''

from impy.common import MCRun, MCEvent, EventGenerator, EventKinematics, standard_particles
import numpy as np

#=========================================================================
# EPOSMCEvent
#=========================================================================


class EPOSMCEvent(MCEvent):
    """Wrapper class around EPOS particle stack."""

    def __init__(self, lib, event_kinematics, frame):
        # HEPEVT (style) common block
        evt = lib.hepevt
        # Number of entries on stack
        npart = evt.nhep
        sel = None

        # Filter stack for charged particles if selected
        if impy_config["event_scope"] == 'charged':
            sel = np.where((evt.isthep[:nhep] == 1) & (
                np.abs(lib.charge_vect(evt.idhep[:nhep])) == 1))
        elif impy_config["event_scope"] == 'stable':
            sel = np.where(evt.isthep[:evt.nhep] == 1)
        else:
            raise Exception("not implemented, yet")

        if 'charge_info' in impy_config and impy_config['charge_info']:
            self.charge = lib.charge_vect(evt.idhep[sel])

        # Save selector for implementation of on-demand properties
        self.sel = sel

        MCEvent.__init__(
            self,
            event_kinematics=event_kinematics,
            lib=lib,
            px=evt.phep[0, sel][0],
            py=evt.phep[1, sel][0],
            pz=evt.phep[2, sel][0],
            en=evt.phep[3, sel][0],
            p_ids=evt.idhep[sel],
            npart=npart,
            frame=frame)

    @property
    def charge(self):
        return lib.charge_vect(evt.idhep[sel])

    @property
    def mass(self):
        return self.lib.hepevt.phep[4, self.sel][0]

    # Nuclear collision parameters

    @property
    def impact_parameter(self):
        """Returns impact parameter for nuclear collisions."""
        return self.lib.nuc3.bimp


#=========================================================================
# EPOSMCRun
#=========================================================================
class EPOSMCRun(MCRun):
    """Implements all abstract attributes of MCRun for the 
    EPOS-LHC series of event generators."""

    def __init__(self, libref, event_class=None, **kwargs):

        if event_class is None:
            self._event_class = EPOSMCEvent
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
        return "EPOS"

    @property
    def version(self):
        """Event generator version"""
        # Needs some sort of smart handling here
        return "LHC"

    def sigma_inel(self):
        """Inelastic cross section according to current
        event setup (energy, projectile, target)"""
        k = self._curr_event_kin
        return self.lib.xsection()[1]

    def _epos_tup(self):
        """Constructs an tuple of arguments for calls to event generator
        from given event kinematics object."""
        k = self._curr_event_kin
        info(20, 'Request EPOS ARGs tuple:\n', 
            (k.ecm, -1., k.p1pdg, k.p2pdg, k.A1, k.Z1, k.A2, k.Z2))
        return (k.ecm, -1., k.p1pdg, k.p2pdg, k.A1, k.Z1, k.A2, k.Z2)

    def set_event_kinematics(self, event_kinematics):
        """Set new combination of energy, momentum, projectile
        and target combination for next event."""
        k = event_kinematics
        self.lib.initeposevt(*self._epos_tup())
        info(5, 'Setting event kinematics')
        self._curr_event_kin = event_kinematics

    def attach_log(self):
        """Routes the output to a file or the stdout."""
        fname = impy_config['output_log']
        if fname == 'stdout':
            lun = 6
            info(5, 'Output is routed to stdout.')
        else:
            lun = self._attach_fortran_logfile(fname)
            info(5, 'Output is routed to', fname, 'via LUN', lun)

        self.lib.files.ifch = lun

    def init_generator(self, event_kinematics, datdir='./iamdata/'):
        self.abort_if_already_initialized()

        self.lib.aaset(0)
        self.lib.initializeepos(1., 1e6, datdir,
                                len(datdir), 1, 2212, 2212, 1, 1, 1, 1, 0,
                                6)
        self.lib.init = True

        # Set default stable
        for pid in [2112, 111, 211, -211, 321, -321, 310, 13, -13]:
            self.set_stable(pid)

        if 'stable' in self.impy_config:
            for pdgid in self.impy_config['stable']:
                self.set_stable(pdgid)

        # Comprise DPMJET input from the event kinematics object
        self.set_event_kinematics(event_kinematics)

    def set_stable(self, pdgid):
        if pdgid not in self.stable_list:
            self.lib.setstable(pdgid)
            self.stable_list.append(pdgid)
        info(5, 'defining', pdgid, 'as stable particle')

    def generate_event(self):
        self.lib.aepos(-1)
        self.lib.afinal()
        self.lib.hepmcstore()
        return False
