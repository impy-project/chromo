
#=========================================================================
# PhoEventViewer
#=========================================================================
# class PhoEventViewer(PhojetRun):
#     def init_generator(self):

#         rejection = self.lib.pho_init(-1, 6)
#         if rejection:
#             raise Exception("PHOJET rejected the initialization!")

#         process_switch = self.lib.poprcs.ipron
#         # non-diffractive elastic scattering (1 - on, 0 - off)
#         process_switch[0, 0] = 1
#         # elastic scattering
#         process_switch[1, 0] = 0
#         # quasi-elastic scattering (for incoming photons only)
#         process_switch[2, 0] = 1
#         # central diffration (double-pomeron scattering)
#         process_switch[3, 0] = 1
#         # particle 1 single diff. dissociation
#         process_switch[4, 0] = 1
#         # particle 2 single diff. dissociation
#         process_switch[5, 0] = 1
#         # double diff. dissociation
#         process_switch[6, 0] = 1
#         # direct photon interaction (for incoming photons only)
#         process_switch[7, 0] = 1

#     def init_event(self):
#         reject = 0
#         self.p1, self.p2 = self.init_beam_setup(self.sqs, self.p1_type,
#                                                 self.p2_type)
#         # initial beam configuration
#         self.settings.enable()
#         print "Initializing PHOJET with", self.p1_type, self.p2_type, self.sqs
#         self.sigmax, self.reject = self.lib.pho_event(-1, self.p1, self.p2,
#                                                       self.sigmax, reject)
#         if reject:
#             raise Exception("Warning: PHOJET rejected the beam configuration.")
#         print self.label, ": Generation 1 event for", self.p1_type, self.p2_type, self.sqs

#     def next_event(self):
#         cursigma = 0.
#         reject = 0
#         cursigma, reject = self.lib.pho_event(1, self.p1, self.p2, cursigma,
#                                               reject)
#         if reject != 0:
#             print "Warning: PHOJET rejected the event configuration."
#             while reject != 0:
#                 cursigma, reject = self.lib.pho_event(1, self.p1, self.p2,
#                                                       cursigma, reject)
#         self.event = self.event_class(self.lib)


# No idea, yet, what to do with this one. It uses some link to
# the fastjet library to do more fancy computations

# class PhoEventJets(PhojetEvent):
#     def find_jets(self, dR, psrapcut):
#         # charged anti-kt jets
#         lib = self.lib
#         lib.pho_findjets(1, 1, dR, psrapcut)
#         njets = lib.pojets.njets

#         self.njets = njets
#         self.j_E = lib.pojets.jets[3, :njets]
#         self.j_p_t = np.sqrt(lib.pojets.jets[0, :njets]**2 +
#                              lib.pojets.jets[1, :njets]**2)
#         self.j_pz = lib.pojets.jets[1, :njets]
#         self.j_p_tot = np.sqrt(lib.pojets.jets[0, :njets]**2 +
#                                lib.pojets.jets[1, :njets]**2 +
#                                lib.pojets.jets[2, :njets]**2)

#         self.j_eta = np.log((self.j_p_tot + self.j_pz) / self.j_p_t)
#         self.j_y = 0.5 * \
#             np.log((self.j_E + self.j_pz) / (self.j_E - self.j_pz))
