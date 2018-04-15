class SIBYLLRun(MCRun):

    def __init__(self, libref,  event_class=None, **kwargs):
        if event_class is None:
            self.event_class = SibyllMCEvent
        else:
            self.event_class = event_class

        MCRun.__init__(self, libref, **kwargs)

    @property
    def generator(self):
        """Event generator name"""
        return "SIBYLL"

    @property
    def version(self):
        """Event generator version"""
        return "2.3"

    def sigma_inel(self):
        k = self.event_config['event_kinematics']
        sigproj = None
        if abs(k.p1pdg) in [2212, 3112]:
            sigproj = 1
        elif abs(k.p1pdg) == 211:
            sigproj = 2
        elif abs(k.p1pdg) == 321:
            sigproj = 3
        else:
            raise Exception(
                self.class_name + "::init_generator(): No " +
                "cross section available for projectile " + k.p1pdg)
        return self.lib.sib_sigma_hp(sigproj, self.ecm)[2]

    def set_event_kinematics(self, event_kinematics):
        k = event_kinematics

        self.sibproj = self.lib.isib_pdg2pid(k.p1pdg)
        self.eatarg = k.A2
        self.ecm = k.ecm

        self.event_config['event_kinematics'] = k
        self.evkin = event_kinematics

    def attach_log(self, fname):
        if fname == 'stdout':
            self.lib.s_debug.lun = 66
        else:
            

    def init_generator(self, config):
        self.abort_if_already_initialized()

        self.set_event_kinematics(self.evkin)

        from random import randint

        try:
            from MCVD.management import LogManager
            # initialize log manager
            self.log_man = LogManager(config, self)
            self.log_man.create_log(self.sibproj, self.eatarg, self.ecm)
            self.lib.s_debug.lun = 66
        except ImportError:
            print self.class_name + "::init_generator(): Running outside of MCVD,", \
                "the log will be printed to STDOUT."
        except AttributeError:
            print self.class_name + "::init_generator(): Logging not supported."

        self.lib.sibini(randint(1000000, 10000000))
        set_stable(self.lib, 2)
        self.lib.pdg_ini()
        self.set_event_kinematics(self.evkin)

        if self.def_settings:
            print self.class_name + "::init_generator(): Using default settings:", \
                self.def_settings.__class__.__name__
            self.def_settings.enable()

        # Set pi0 stable
        self.set_stable(111)

        if 'stable' in self.event_config:
            for pdgid in self.event_config['stable']:
                self.set_stable(pdgid)

    def set_stable(self, pdgid):
        sid = abs(self.lib.isib_pdg2pid(pdgid))
        idb = self.lib.s_csydec.idb
        if sid == 0 or sid > idb.size - 1:
            return
        print self.class_name + "::set_stable(): defining ", \
            pdgid, "as stable particle, sid =", sid
        idb[sid - 1] = -np.abs(idb[sid - 1])

    def generate_event(self):
        self.lib.sibyll(self.sibproj, self.eatarg, self.ecm)
        self.lib.decsib()
        return 0  # SIBYLL never rejects


# class SibyllCascadeEvent():
#     def __init__(self, lib):
#         npart = lib.s_plist.np
#         stable = np.where(np.abs(lib.s_plist.llist[:npart]) < 10000)[0]

#         self.p_ids = lib.s_plist.llist[:npart][stable]
#         self.E = lib.s_plist.p[:npart, 3][stable]
#         self.pz = lib.s_plist.p[:npart, 2][stable]


# class SibyllCascadeEventPQCDc():
#     def __init__(self, lib):
#         npart = lib.s_plist.np
#         spl = lib.s_plist
#         stable = np.where((np.abs(spl.llist[:npart]) < 10000) & np.logical_not(
#             (np.abs(spl.llist[:npart]) >= 59) & (lib.s_chist.jdif[0] == 0) & (
#                 (np.abs(lib.s_parto.nporig[:npart]) < 10) |
#                 (np.abs(lib.s_parto.nporig[:npart]) > 1000))))[0]

#         self.p_ids = spl.llist[:npart][stable]
#         # if np.any(self.p_ids >= 59):

#         #     for pid, co in zip(self.p_ids, lib.s_parto.nporig[:npart][stable]):
#         #         if abs(pid) < 59:
#         #             continue

#         #         print pid, co
#         #     print self.p_ids
#         self.E = spl.p[:npart, 3][stable]
#         self.pz = spl.p[:npart, 2][stable]


# class SibyllCascadeRun():
#     def __init__(self,
#                  lib_str,
#                  label,
#                  decay_mode,
#                  n_events,
#                  fill_subset=False,
#                  evt_class=SibyllCascadeEvent):
#         from ParticleDataTool import SibyllParticleTable
#         exec "import " + lib_str + " as siblib"
#         self.lib = siblib  # @UndefinedVariable
#         self.label = label
#         self.nEvents = n_events
#         self.fill_subset = fill_subset
#         self.evt_class = evt_class
#         self.spectrum_hists = []
#         if decay_mode > 2 and lib_str.find('sib21') != -1:
#             print "Limiting decay mode to 3 for SIBYLL 2.1"
#             decay_mode = 3

#             self.lib.s_debug.lun = 6
#         set_stable(self.lib, decay_mode)
#         self.ptab = SibyllParticleTable()
#         self.init_generator()

#     def init_generator(self):
#         from random import randint
#         seed = randint(1000000, 10000000)
#         print self.__class__.__name__ + '::init_generator(): seed=', seed
#         self.lib.sibini(seed)

#     def get_hadron_air_cs(self, E_lab, proj_name):
#         sqs = np.sqrt(2 + 2 * E_lab)

#         if proj_name.find('pi') != -1:
#             proj_name = 2
#         elif proj_name.find('K') != -1:
#             proj_name = 3
#         else:
#             proj_name = 1
#         try:
#             return self.lib.sib_sigma_hair(proj_name, sqs)[0]
#         except:
#             return self.lib.sib_sigma_hair(proj_name, sqs)

#     def get_hadron_p_cs(self, E_lab, proj_name):
#         sqs = np.sqrt(2 + 2 * E_lab)

#         if proj_name.find('pi') != -1:
#             proj_name = 2
#         elif proj_name.find('K') != -1:
#             proj_name = 3
#         else:
#             proj_name = 1

#         return self.lib.sib_sigma_hp(proj_name, sqs)[2]

#     def start(self, projectile, E_lab, sqs, Atarget=0):
#         templ = '''
#         Model     : {0}
#         N_events  : {1}
#         Projectile: {2}
#         E_lab     : {3:5.3e} GeV
#         E_cm      : {4:5.3e} GeV
#         Target    : {5}
#         '''
#         print templ.format(self.label, self.nEvents, projectile, E_lab, sqs,
#                            Atarget)

#         if projectile not in self.projectiles:
#             raise Exception("SibyllCascadeRun(): Projectile " + projectile +
#                             " not allowed.")
#         projectile = self.ptab.modname2modid[projectile]

#         hist_d = {}
#         for hist in self.spectrum_hists:
#             hist_d[hist.particle_id] = hist
#         ngenerated = self.nEvents
#         for i in xrange(self.nEvents):
#             self.lib.sibyll(projectile, Atarget, sqs)
#             self.lib.decsib()

#             if not (i % 10000) and i:
#                 print i, "events generated."
#             event = self.evt_class(self.lib)

#             unique_pids = np.unique(event.p_ids)
#             if 0 in unique_pids:
#                 self.lib.sib_list(6)
#                 ngenerated = ngenerated - 1
#                 continue
#             if not self.fill_subset:
#                 [hist_d[pid].fill_event(event) for pid in unique_pids]
#             else:
#                 for pid in unique_pids:
#                     if pid in hist_d.keys():
#                         hist_d[pid].fill_event(event)
#         # Correct for selective filling of histograms
#         for hist in self.spectrum_hists:
#             hist.n_events_filled = ngenerated

#     def generate_decay_spectrum(self, sibyll_particle_id, lab_energy):
#         sibid = sibyll_particle_id
#         asibid = abs(sibid)
#         llist = self.lib.s_plist.llist
#         p = self.lib.s_plist.p
#         lmass = self.lib.s_mass1.am
#         idb = self.lib.s_csydec.idb
#         # Save previous decay flag and set decaying particle unstable
#         save_decay_status = idb[asibid - 1]
#         idb[asibid - 1] = abs(idb[asibid - 1])
#         P0 = np.zeros(5)
#         P0[0] = 0.
#         P0[1] = 0.
#         P0[2] = lab_energy
#         P0[4] = lmass[asibid - 1]
#         P0[3] = np.sqrt(P0[2]**2 + P0[4]**2)
#         for i in range(self.nEvents):
#             llist[0] = sibid
#             p[0, :] = P0
#             self.lib.s_plist.np = 1
#             #             self.lib.sib_list(6)
#             self.lib.decsib()
#             if not (i % 10000) and i:
#                 print i, "decay events simulated."
#             event = SibyllCascadeEvent(self.lib)
#             for i, pid in enumerate(event.p_ids):
#                 try:
#                     self.spectrum_hists[pid].fill_event(event, i)
#                 except:
#                     pass

#         for hist in self.spectrum_hists.itervalues():
#             hist.n_events_filled = self.nEvents
#             hist.finalize_run()

#         idb[asibid - 1] = save_decay_status