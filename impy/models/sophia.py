'''
Created on 17.03.2014

@author: afedynitch
'''

import numpy as np
from impy.common import MCRun, MCEvent, impy_config
from impy.util import standard_particles, info
import particletools
from particletools.tables import SibyllParticleTable
import time


class SophiaCascadeRun():
    def __init__(self,
                 lib_str,
                 label,
                 decay_mode,
                 n_events,
                 fill_subset=False,
                 p_debug=False,
                 nucleon_Ekin=None):
        from ParticleDataTool import SibyllParticleTable
        exec("import " + lib_str + " as siblib")
        self.lib = siblib  # @UndefinedVariable
        self.label = label
        self.nEvents = n_events
        self.fill_subset = fill_subset
        self.spectrum_hists = []
        self.dbg = p_debug
        self.nucleon_Ekin = nucleon_Ekin
        set_stable(self.lib, decay_mode, self.dbg)
        self.ptab = SibyllParticleTable()
        self.init_generator()

    def init_generator(self):
        from random import randint
        seed = randint(1000000, 10000000)
        print(self.__class__.__name__ + '::init_generator(): seed=', seed)
        self.lib.init_rmmard(seed)

    def get_hadron_air_cs(self, E_lab, projectile_sibid):
        raise Exception('Not implemented, yet')

    def start(self, evkin):

        swap = False
        nucleon_id = None
        nucleon_mass = None
        if evkin.p1pdg == 22 and evkin.p2pdg in [2212, 2112]:
            swap = True
            nucleon_id = self.ptab.pdg2modid[evkin.p2pdg]
            nucleon_mass = evkin.pmass2
        elif evkin.p2pdg == 22 and evkin.p1pdg in [2212, 2112]:
            swap = False
            nucleon_id = self.ptab.pdg2modid[evkin.p1pdg]
            nucleon_mass = evkin.pmass1

        swap = False
        nucleon_e = nucleon_mass + self.nucleon_Ekin

        hist_d = {}
        for hist in self.spectrum_hists:
            hist_d[hist.particle_id] = hist
        ngenerated = self.nEvents

        for i in range(self.nEvents):
            self.lib.eventgen(nucleon_id, nucleon_e, evkin.elab, 180., 0)
            if not (i % 10000) and i and self.dbg:
                print(i, "events generated.")

            event = SophiaCascadeEvent(self.lib, swap)

            unique_pids = np.unique(event.p_ids)
            if 0 in unique_pids:
                self.lib.sib_list(6)
                ngenerated = ngenerated - 1
                continue
            if not self.fill_subset:
                [hist_d[pid].fill_event(event) for pid in unique_pids]
            else:
                for pid in unique_pids:
                    if pid in list(hist_d.keys()):
                        hist_d[pid].fill_event(event)
        # Correct for selective filling of histograms
        for hist in self.spectrum_hists:
            hist.n_events_filled = ngenerated


class SophiaCascadeEvent():
    def __init__(self, lib, swap=False):
        npart = lib.s_plist.np
        stable = np.nonzero(np.abs(lib.s_plist.llist[:npart]) < 10000)[0]

        self.p_ids = lib.s_plist.llist[:npart][stable]
        self.E = lib.s_plist.p[:npart, 3][stable]
        if swap:
            self.pz = -lib.s_plist.p[:npart, 2][stable]
        else:
            self.pz = lib.s_plist.p[:npart, 2][stable]


#=========================================================================
# set_stable
#=========================================================================
def set_stable(lib, decay_mode, dbg=True):
    idb = lib.s_csydec.idb

    if dbg:
        print("SophiaCascadeRun::set_stable(): Setting standard" +
              " particles stable.")

    #fast-mode particles
    if decay_mode == 0:
        stab = SibyllParticleTable()
        for pdg_id in standard_particles:
            print('stable,', pdg_id)
            idb[i - 1] = -np.abs(idb[i - 1])
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
    if dbg:
        print("SophiaCascadeRun::set_stable(): Setting conventional " +
              "Sigma-, Xi0, Xi- and Lambda0 stable (decay mode).")
    for i in range(36, 39 + 1):
        idb[i - 1] = -np.abs(idb[i - 1])

    if decay_mode <= 2:
        return

    # Conventional mesons and baryons
    # keep eta, eta', rho's, omega, phi, K*
    if dbg:
        print("SophiaCascadeRun::set_stable(): Setting all " +
              "conventional stable.")
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
    if dbg:
        print("SophiaCascadeRun::set_stable(): Setting all " +
              "conventional and charmed stable.")
    for i in list(range(59, 61)) + list(range(71, 99 + 1)):
        idb[i - 1] = -np.abs(idb[i - 1])



class SophiaEvent(MCEvent):
    """Wrapper class around SIBYLL 2.1 & 2.3 particle stack."""
    # Workaround for no data on vertext positions in SIBYLL
    _no_vertex_data = None

    def __init__(self, lib, event_kinematics, event_frame):   
        # self.lib = lib
        # #lib.toevt()
        start_time = time.time()
        lib.prepare_event_data()
        time1 = time.time() - start_time
        start_time = time.time()
        lib.toevt()
        time2 = time.time() - start_time
        print("time1/time2 = %s" % (time1/time2))
        time.sleep(1)
        #hepevt = lib.hepevt
        
        # event_number = hepevt.nevhep
        number_of_particles = lib.hepevt.nhep
        # pdg_id = hepevt.idhep[0:number_of_particles]
        # status_code = hepevt.isthep[0:number_of_particles]
        # particle_momenta = hepevt.phep[:, 0:number_of_particles]
        # vertex = hepevt.vhep[:, 0:number_of_particles]
        
        # common /event_data/ event_number,number_of_particles,
        # status_codes,pdg_ids,charges,parents particle_momenta,vertices,jmohep,jdahep
        # event_number = lib.event_data.event_number
        #number_of_particles = lib.event_data.number_of_particles
        
        # particle_momenta = lib.event_data.particle_momenta[:, 0:number_of_particles]
        # px = particle_momenta[0]
        # py = particle_momenta[1]
        # pz = particle_momenta[2]
        # en = particle_momenta[3]
        # m =  particle_momenta[4]
        
        arr = lib.s_plist.p[:number_of_particles,:].T
        px = arr[0]
        py = arr[1]
        pz = arr[2]
        en = arr[3]
        ma = arr[4]
        
        # pdg_id = lib.event_data.pdg_ids[0:number_of_particles]
        # status_code = lib.event_data.status_codes[0:number_of_particles]
        #print(particle_momenta)
        
        # buf = lib.s_plist.p.data
        # px = lib.s_plist.p[:number_of_particles,0:1].T[0]
        # py = lib.s_plist.p[:number_of_particles,1:2].T[0]
        # pz = lib.s_plist.p[:number_of_particles,2:3].T[0]
        # en = lib.s_plist.p[:number_of_particles,3:4].T[0]
        # ma = lib.s_plist.p[:number_of_particles,4:5].T[0]
        
        # print('px =', px)
        # print('py =', py)
        # print('pz =', pz)
        # print('en =', en)
        # print('ma =', ma)
        # vertex = lib.event_data.vertices[:, 0:number_of_particles]
        
        #schg = lib.schg
        #self.decayed_parent = schg.decpar[0:number_of_particles]
        #self.particle_charge = schg.ichg[0:number_of_particles]
        
        #self.hepevt = hepevt
        #self.schg = schg
        
        
        # start_time = time.time()
        # # arr = np.array([lib.event_data.px, lib.event_data.py, 
        # #           lib.event_data.pz, lib.event_data.en,
        # #           lib.event_data.ma])[:, 0:number_of_particles]
        # arr = lib.s_plist.p[:number_of_particles,:].T
        # px = arr[0]
        # py = arr[1]
        # pz = arr[2]
        # en = arr[3]
        # ma = arr[4]

        # time_to_unpack = time.time() - start_time
        # #print("Unpack--- %s seconds ---" % (time.time() - start_time))
        # # px, py, pz, en, m = particle_momenta
        # start_time = time.time()
        # particle_momenta = lib.event_data.particle_momenta[:, 0:number_of_particles]
        # px = particle_momenta[0]
        # py = particle_momenta[1]
        # pz = particle_momenta[2]
        # en = particle_momenta[3]
        # m =  particle_momenta[4]
        # time_to_array = time.time() - start_time
        # print("Array--- %s seconds ---" % (time_to_unpack/time_to_array))
        # time.sleep(1)
        # vx, vy, vz, vt = vertex
        
    
        # MCEvent.__init__(self,
        #                  lib = lib,
        #                  event_kinematics = event_kinematics,
        #                  event_frame = event_frame,
        #                  nevent = event_number,
        #                  npart = number_of_particles,
        #                  p_ids = pdg_id,
        #                  status = status_code,
        #                  px = px,
        #                  py = py,
        #                  pz = pz,
        #                  en = en,
        #                  m = m,
        #                  vx = vx,
        #                  vy = vy,
        #                  vz = vz,
        #                  vt = vt,
        #                  pem_arr = particle_momenta,
        #                  vt_arr = vertex)

    def filter_final_state(self):
        self.selection = np.where(self.status == 1)
        self._apply_slicing()

    def filter_final_state_charged(self):
        self.selection = np.where((self.status == 1) & (self.charge != 0))
        self._apply_slicing()
    
    @property
    def charge(self):
        return self.particle_charge[self.selection]
    
    @property
    def decayed_parent(self):
        MCEvent.parents(self)
        return self.decayed_parent

    @property
    def parents(self):
        """In SOPHIA parents are difficult to obtain. This function returns 0."""
        MCEvent.parents(self)
        return self.hepevt.jmohep

    @property
    def children(self):
        """In SOPHIA daughters are difficult to obtain. This function returns 0."""
        MCEvent.children(self)
        return self.hepevt.jdahep

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

class SophiaRun(MCRun):
    """Implements all abstract attributes of MCRun for the 
    Sophia event generator.
    
    """

    def sigma_inel(self, *args, **kwargs):
        # """Inelastic cross section according to current
        # event setup (energy, projectile, target)"""
        # k = self._curr_event_kin
        # sigproj = None
        # if abs(k.p1pdg) in [2212, 2112, 3112]:
        #     sigproj = 1
        # elif abs(k.p1pdg) == 211:
        #     sigproj = 2
        # elif abs(k.p1pdg) == 321:
        #     sigproj = 3
        # else:
        #     info(0, "No cross section available for projectile", k.p1pdg)
        #     raise Exception('Input error')
        
        # if k.p1_is_nucleus:
        #     raise Exception('Nuclear projectiles not supported by SIBYLL.')

        # if k.p2_is_nucleus:
        #     # Return production cross section for nuclear target
        #     try:
        #         return self.lib.sib_sigma_hnuc(sigproj, k.A2, self._ecm)[0]
        #     except AttributeError:
        #         return 'Nuclear cross section not supported for this SIBYLL version'
        
        # return self.lib.sib_sigma_hp(sigproj, self._ecm)[2]
        # self.energy_of_photon is the energy in lab frame
        # where nucleon is at rest and photon is moving
        total_crossection_id = 3 # 3 is for total crossection
        # cross section in micro barn
        return self.lib.crossection(self.energy_of_photon, total_crossection_id, 
                                    self.nucleon_code_number)
        

    def sigma_inel_air(self):
        """Inelastic cross section according to current
        event setup (energy, projectile, target)"""
        raise Exception('SophiaRun.sigma_inel_air has no implementation')

    def set_event_kinematics(self, event_kinematics):
        """Set new combination of energy, momentum, projectile
        and target combination for next event."""

        info(5, 'Setting event kinematics.')
        info(10, event_kinematics)
        k = event_kinematics
        if k.p1_is_nucleus:
            raise Exception('Projectile nuclei not natively supported in SIBYLL')
        elif k.p2_is_nucleus and k.A2 > 20:
            print(k.p2_is_nucleus, k.A2)
            raise Exception('Target nuclei with A>20 not supported in SIBYLL')
        
        
        self.nucleon_code_number = self.lib.icon_pdg_sib(k.p2pdg)
        self.energy_of_nucleon = k.pmass2
        self.energy_of_photon = k.elab 
        # Here we consider laboratory frame where photon moves along z axis
        # and nucleon is at rest. The angle is counted from z axis.
        # However, because of the definitions in "eventgen" subroutine of
        # SOPHIA code (line "P_gam(3) = -EPS*COS(theta*pi/180.D0)")
        # this angle should be 180 for photon moving along z
        # (and 0 for photon moving in direction opposite to z)
        self.angle_between_nucleon_and_photon = 180
        # This is the output parameter which we initialize with 0
        self.code_of_interaction_mode = 0
        self._curr_event_kin = event_kinematics
        # Keep decayed particles in the history:
        # self.lib.eg_io.remdec == 0 removes all decayed particles
        # self.lib.eg_io.remdec != 0 leaves decayed particles
#        self.lib.eg_io.remdec = 1

    def attach_log(self, fname=None):
        """Routes the output to a file or the stdout."""
        fname = impy_config['output_log'] if fname is None else fname
        if fname == 'stdout':
            #self.lib.s_debug.lun = 6
            info(5, 'Output is routed to stdout.')
        else:
            lun = self._attach_fortran_logfile(fname)
            #self.lib.s_debug.lun = lun
            info(5, 'Output is routed to', fname, 'via LUN', lun)

    def init_generator(self, event_kinematics, seed='random', logfname=None):
        from random import randint

        self._abort_if_already_initialized()

        if seed == 'random':
            seed = randint(1000000, 10000000)
        else:
            seed = int(seed)
        info(5, 'Using seed:', seed)

        self.lib.s_plist.ideb = impy_config['sibyll']['debug_level']

        self.set_event_kinematics(event_kinematics)
        self.attach_log(fname=logfname)
        
        self.lib.init_rmmard(int(seed)) # setting random number generator seed
        self.lib.initial(self.nucleon_code_number) # setting parameters for cross-section
        #self.lib.pdg_ini()
        #self.conv_hepevt = (self.lib.sibhep1
        #                    if '21' in self.lib.__name__ else self.lib.sibhep3)
        self._define_default_fs_particles()

    def set_stable(self, pdgid, stable=True):
        sid = abs(self.lib.icon_pdg_sib(pdgid))
        if abs(pdgid) == 311:
            info(1, 'Ignores K0. Use K0L/S 130/310 in final state definition.')
            return
        idb = self.lib.s_csydec.idb
        if sid == 0 or sid > idb.size - 1:
            return
        if stable:
            info(
                5, 'defining as stable particle pdgid/sid = {0}/{1}'.format(
                    pdgid, sid))
            idb[sid - 1] = -np.abs(idb[sid - 1])
        else:
            info(5, 'pdgid/sid = {0}/{1} allowed to decay'.format(pdgid, sid))
            idb[sid - 1] = np.abs(idb[sid - 1])

    def generate_event(self):
        # Generate event (the final particles and their parameters) 
        # by underlying Fortran library 
        self.lib.eventgen(self.nucleon_code_number, 
                self.energy_of_nucleon, 
                self.energy_of_photon, 
                self.angle_between_nucleon_and_photon, 
                self.code_of_interaction_mode)
        return 0  # No rejection is implemented so far
    