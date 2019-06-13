'''
Created on 03.05.2016

@author: afedynitch
'''

import numpy as np
from impy.common import MCRun, MCEvent, impy_config, root_dir
from impy.util import standard_particles, info

class EPOSEvent(MCEvent):
    """Wrapper class around EPOS particle stack."""

    def __init__(self, lib, event_kinematics, event_frame):
        # HEPEVT (style) common block
        evt = lib.hepevt

        # Save selector for implementation of on-demand properties
        px, py, pz, en, m = evt.phep
        vx, vy, vz, vt = evt.vhep

        MCEvent.__init__(
            self,
            lib=lib,
            event_kinematics=event_kinematics,
            event_frame=event_frame,
            nevent=evt.nevhep,
            npart=evt.nhep,
            p_ids=evt.idhep,
            status=evt.isthep,
            px=px,
            py=py,
            pz=pz,
            en=en,
            m=m,
            vx=vx,
            vy=vy,
            vz=vz,
            vt=vt,
            pem_arr=evt.phep,
            vt_arr=evt.vhep)

    def filter_final_state(self):
        self.selection = np.where(self.status == 1)
        self._apply_slicing()

    def filter_final_state_charged(self):
        self.selection = np.where((self.status == 1) & (self.charge != 0))
        self._apply_slicing()

    def filter_final_state_centrality(self, event_kinematics, centrality_range):
        '''Filters all events within the selected centrality bins, defined by the impact parameter and the radius of the nucleus (R = 1.2*A^(1/3)fm; works with error < 1%). This calculation works for all except for peripheral centrality bins (> 50%). 
        Source: https://www.hindawi.com/journals/ahep/2013/908046/
        This only works for MC models that support NN collisions (that have impact parameter calculations). 
        '''
        centrality_val = (self.impact_parameter)**2 / (2*1.2*(event_kinematics.A1**(1/3)))**2

        self.selection = np.where((self.status == 1) &
         (self.charge != 0) & (centrality_val > centrality_range[0]) &(centrality_val < centrality_range[1]))
        self._apply_slicing()

    @property
    def parents(self):
        if self._is_filtered:
            raise Exception(
                'Parent indices do not point to the' +
                ' proper particles if any slicing/filtering is applied.')
        return self.lib.hepevt.jmohep

    @property
    def children(self):
        if self._is_filtered:
            raise Exception(
                'Parent indices do not point to the' +
                ' proper particles if any slicing/filtering is applied.')
        return self.lib.hepevt.jdahep

    @property
    def charge(self):
        return self.lib.charge_vect(self.lib.hepevt.idhep[self.selection])

    # Nuclear collision parameters
    @property
    def impact_parameter(self):
        """Returns impact parameter for nuclear collisions."""
        return self.lib.nuc3.bimp


class EPOSRun(MCRun):
    """Implements all abstract attributes of MCRun for the 
    EPOS-LHC series of event generators."""

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
        self._curr_event_kin = k
        self.lib.initeposevt(*self._epos_tup())
        info(5, 'Setting event kinematics')

    def attach_log(self):
        """Routes the output to a file or the stdout."""
        fname = impy_config['output_log']
        if fname == 'stdout':
            lun = 6
            info(5, 'Output is routed to stdout.')
        else:
            lun = self._attach_fortran_logfile(fname)
            info(5, 'Output is routed to', fname, 'via LUN', lun)

        self._lun = lun


    def init_generator(self, event_kinematics, seed='random'):
        from random import randint
        from os import path

        self._abort_if_already_initialized()

        if seed == 'random':
            seed = randint(1000000, 10000000)
        else:
            seed = int(seed)
        info(5, 'Using seed:', seed)
        
        epos_conf = impy_config['epos']
        datdir = path.join(root_dir, epos_conf['datdir'])
        self.attach_log()
        info(1, 'First initialization')
        self.lib.aaset(0)
        self.lib.initializeepos(
            float(seed), 1e6, datdir, len(datdir), 1, 2212, 2212, 1, 1, 1, 1,
            0, self._lun)

        # Set default stable
        self._define_default_fs_particles()
        self.set_event_kinematics(event_kinematics)

        self.lib.charge_vect = np.vectorize(self.lib.getcharge)

    def set_stable(self, pdgid, stable=True):
        if stable:
            self.lib.setstable(pdgid)
            info(5, 'defining', pdgid, 'as stable particle')
        else:
            self.lib.setunstable(pdgid)
            info(5, pdgid, 'allowed to decay')

    def generate_event(self):
        self.lib.aepos(-1)
        self.lib.afinal()
        self.lib.hepmcstore()
        return False
