'''
Created on 17.03.2014

@author: afedynitch
'''

import numpy as np
from impy.common import MCRun, MCEvent, impy_config
from impy.util import standard_particles, info

class SibyllMCEvent(MCEvent):
    """Wrapper class around SIBYLL 2.1 & 2.3 particle stack."""

    def __init__(self, lib, event_config):
        # The stack common block
        s_plist = lib.s_plist
        # Number of entries on stack
        npart = lib.s_plist.np
        # Conversion to PDG ids from SIBYLL routine
        to_pdg = np.vectorize(lib.isib_pid2pdg)

        # Filter stack for charged particles if selected
        stable = np.nonzero(np.abs(s_plist.llist[:npart]) < 10000)[0]
        if impy_config["event_scope"] == 'charged':
            sel = stable[
                lib.s_chp.ichp[np.abs(s_plist.llist[stable]) - 1] != 0]
        elif impy_config["event_scope"] == 'stable':
            sel = stable
        else:
            raise Exception("not implemented, yet")

        # Save selector for implementation of on-demand properties
        self.sel = sel

        MCEvent.__init__(
            self,
            event_config=event_config,
            lib=lib,
            px=s_plist.p[sel, 0],
            py=s_plist.p[sel, 1],
            pz=s_plist.p[sel, 2],
            en=s_plist.p[sel, 3],
            p_ids=to_pdg(s_plist.llist[sel]),
            npart=npart)

    @property
    def charge(self):
        return self.lib.s_chp.ichp[self.lib.s_plist.llist[self.sel] - 1]

    @property
    def mass(self):
        return self.lib.s_plist.p[self.sel, 4]


class SIBYLLRun(MCRun):
    """Implements all abstract attributes of MCRun for the 
    SIBYLL 2.1, 2.3 and 2.3c event generators."""
        
    def __init__(self, libref, event_class=None, **kwargs):
        if event_class is None:
            self.event_class = SibyllMCEvent
        else:
            self.event_class = event_class

        MCRun.__init__(self, libref, **kwargs)

    @property
    def name(self):
        """Event generator name"""
        return "SIBYLL"

    @property
    def version(self):
        """Event generator version"""
        return "2.3"

    def sigma_inel(self):
        """Inelastic cross section according to current
        event setup (energy, projectile, target)"""
        k = self.evkin
        sigproj = None
        if abs(k.p1pdg) in [2212, 3112]:
            sigproj = 1
        elif abs(k.p1pdg) == 211:
            sigproj = 2
        elif abs(k.p1pdg) == 321:
            sigproj = 3
        else:
            info(0, "No cross section available for projectile", k.p1pdg)
            raise Exception('Input error')

        return self.lib.sib_sigma_hp(sigproj, self.ecm)[2]

    def set_event_kinematics(self, event_kinematics):
        """Set new combination of energy, momentum, projectile
        and target combination for next event."""

        k = event_kinematics

        self.sibproj = self.lib.isib_pdg2pid(k.p1pdg)
        self.eatarg = k.A2
        self.ecm = k.ecm

        self.event_config['event_kinematics'] = k
        self.evkin = event_kinematics

    def attach_log(self):
        """Routes the output to a file or the stdout."""
        fname = impy_config['output_log']
        if fname == 'stdout':
            self.lib.s_debug.lun = 6
            info(5, 'Output is routed to stdout.')
        else:
            lun = self._attach_fortran_logfile(fname)
            self.lib.s_debug.lun = lun
            info(5, 'Output is routed to', fname, 'via LUN', lun)

    def init_generator(self):
        from random import randint
        self.set_event_kinematics(self.evkin)

        self.abort_if_already_initialized()

        self.attach_log(self.output_file)

        self.lib.sibini(randint(1000000, 10000000))
        set_stable(self.lib, 2)
        self.lib.pdg_ini()

        self.set_stable(111)

        if 'stable' in self.event_config:
            for pdgid in self.event_config['stable']:
                self.set_stable(pdgid)

    def set_stable(self, pdgid):
        sid = abs(self.lib.isib_pdg2pid(pdgid))
        idb = self.lib.s_csydec.idb
        if sid == 0 or sid > idb.size - 1:
            return
        info(1, 'defining as stable particle',
             'pdgid/sid = {0}/{1}'.format(pdgid, sid))
        idb[sid - 1] = -np.abs(idb[sid - 1])

    def generate_event(self):
        self.lib.sibyll(self.sibproj, self.eatarg, self.ecm)
        self.lib.decsib()
        return True  # SIBYLL never rejects


#=========================================================================
# set_stable
#=========================================================================
def set_stable(lib, decay_mode):
    from particletools.tables import SibyllParticleTable
    idb = lib.s_csydec.idb

    info(1,"Setting standard particles stable.")

    if decay_mode < 0:
        info(1, "use default stable def.")
        return

    # fast-mode particles
    if decay_mode == 0:
        # Set all instable
        for i in range(4, 13):
            idb[i - 1] = np.abs(idb[i - 1])
        for i in range(23, 100):
            idb[i - 1] = np.abs(idb[i - 1])

        stab = SibyllParticleTable()
        for pdg_id in standard_particles:
            idb[stab.pdg2modid[pdg_id] - 1] = \
                -np.abs(idb[stab.pdg2modid[pdg_id] - 1])
        return

    # keep muons pions, kaons
    for i in range(4, 5 + 1):
        idb[i - 1] = -np.abs(idb[i - 1])
    for i in range(7, 18 + 1):
        idb[i - 1] = -np.abs(idb[i - 1])
    # K0 and K0-bar have to remain unstable to form K0S/L

    if decay_mode <= 1:
        return

    # Decay mode 2 for generation of decay spectra (all conventional with
    # lifetime >= K0S
    info(1, "Setting conventional Sigma-, Xi0,", 
    "Xi- and Lambda0 stable (decay mode).")
    for i in range(36, 39 + 1):
        idb[i - 1] = -np.abs(idb[i - 1])

    if decay_mode <= 2:
        return

    # Conventional mesons and baryons
    # keep eta, eta', rho's, omega, phi, K*
    info(1, "Setting all conventional stable.")
    # pi0
    idb[6 - 1] = -np.abs(idb[6 - 1])
    for i in range(23, 33 + 1):
        idb[i - 1] = -np.abs(idb[i - 1])

    # keep SIGMA, XI, LAMBDA
    for i in range(34, 49 + 1):
        idb[i - 1] = -np.abs(idb[i - 1])

    if decay_mode <= 3:
        return

    # Charmed particles (only for version >= 2.2)
    # keep all charmed
    info(1, "Setting all conventional and charmed stable.")
    for i in range(59, 61) + range(71, 99 + 1):
        idb[i - 1] = -np.abs(idb[i - 1])
