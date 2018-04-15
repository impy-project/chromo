class MCRun():
    __metaclass__ = ABCMeta
    def __init__(self, libref, label, n_events, DEBUG, event_class,
                 event_config, event_kinematics, default_settings, settings):
        self.lib = libref
        self.label = label
        self.n_events = n_events
        self.triggers = {}
        self.generator = None
        self.version = None
        self.debug = DEBUG
        self.output_unit = -1
        self.event_class = event_class
        self.event_config = event_config
        self.evkin = event_kinematics
        self.class_name = self.__class__.__name__
        self._is_initialized = False

        def_settings_class, def_settings_args = default_settings
        if def_settings_class:
            self.def_settings = def_settings_class(self.lib, def_settings_args)
        else:
            self.def_settings = Settings(self.lib)

        settings_class, settings_args = settings
        if settings_class:
            self.settings = settings_class(self.lib, settings_args)
            print self.class_name + "::__init__(): Attaching custom settings", \
                  self.settings.get_label()
        else:
            self.settings = Settings(self.lib)

        # Make sure that only charged particles are kept
        if not self.event_config:
            self.event_config = {'charged_only': True}
        if 'charged_only' not in self.event_config:
            self.event_config['charged_only'] = True

        self.nondef_stable = None
        if 'stable' in self.event_config.keys():
            self.nondef_stable = np.sum(
                np.abs(np.array(self.event_config['stable'])))
        self.event_config['event_kinematics'] = self.evkin

    @abstractmethod
    def init_generator(self, config):
        pass

    @abstractmethod
    def generate_event(self):
        pass

    @abstractmethod
    def set_event_kinematics(self, evtkin):
        pass

    @abstractmethod
    def set_stable(self, pdgid):
        pass

    @abstractmethod
    def get_sigma_inel(self):
        pass

    def abort_if_already_initialized(self):
        assert not self._is_initialized
        self._is_initialized = True

    def get_model_label(self):
        return self.label

    def attach_log(self, log_fname):
        pass

    def deattach_log(self):
        pass

    def trigger_event(self, event):
        return [
            trigger for trigger in self.triggers.itervalues()
            if trigger.trigger_event(event)
        ]

    def _init_progress_bar(self, maxval=None):
        if maxval == None:
            maxval = self.n_events + 1
        try:
            from progressbar import ProgressBar, Percentage, Bar, ETA
        except ImportError:
            print "Failed to import 'progressbar' -- disabling progress indicator."
            print "Install the module with 'easy_install progressbar', or",
            print "get it from http://qubit.ic.unicamp.br/~nilton"
            raise Exception("It's easy do do this...")
        self.progressBar = ProgressBar(
            maxval=maxval, widgets=[Percentage(), ' ', Bar(), ' ', ETA()])

    def start(self):
        print "running start()"
        # Look for a definition of energy sweep the event kinematics object
        if self.evkin.e_range:
            self.start_variable()
        else:
            # Start fixed energy run
            self.start_fixed()

    def start_fixed(self):
        print "running start_fixed()"
        if not len(self.triggers):
            raise Exception('No trigger defined in ' + self.class_name + "!")

        print self.class_name + \
        "::start_fixed(): starting generation of", self.n_events, "events"
        print self.evkin
        self.settings.enable()

        self._init_progress_bar()
        self.progressBar.start()

        rej_counter = 0
        for i in xrange(self.n_events):  # @UnusedVariable
            reject = self.generate_event()
            self.progressBar.update(self.progressBar.currval + 1)
            if reject:
                rej_counter += 1
                continue

            event = self.event_class(self.lib, self.event_config)
            fired_triggers = self.trigger_event(event)
            if not fired_triggers:
                continue

            active_histograms = [
                hist
                for trigger in fired_triggers
                for hist in trigger.histogram_list.values()
            ]
            [hist.fill_event(event) for hist in active_histograms]

        # Finalize run
        sigmax = self.get_sigma_inel()
        sigma_gen = sigmax * (float(self.n_events - rej_counter)) / \
            float(self.n_events)

        for trigger in self.triggers.itervalues():
            for hist in trigger.histogram_list.itervalues():
                hist.finalize_run(sigma_gen)

        self.progressBar.finish()
        self.settings.reset()
        try:
            self.lib.pho_event(-2, self.p1, self.p2)[1]
        except:
            pass
        print(r"...completed {0:3.2f}\% of the events have been rejected."
              ).format(100. * float(rej_counter) / float(self.n_events))

        try:
            self.log_man.close_log()
        except AttributeError:
            print 'No logging facility defined.'

    def single_event(self):
        #        print "single_event()"
        if not len(self.triggers):
            raise Exception('No trigger defined in ' + self.class_name + "!")

#        print self.evkin

#        self.settings.enable()

        self.generate_event()
        event = self.event_class(self.lib, self.event_config)

        fired_triggers = self.trigger_event(event)
        if not fired_triggers:
            return

        active_histograms = [
            hist
            for trigger in fired_triggers
            for hist in trigger.histogram_list.values()
        ]
        [hist.fill_event(event) for hist in active_histograms]

        # Finalize run

    #        sigmax = self.get_sigma_inel()
    #        print sigmax

    def start_variable(self):
        if not len(self.triggers):
            raise Exception('No trigger defined in ' + self.class_name + "!")

        print self.class_name + \
        "::start_variable(): starting generation of", self.n_events, "events"
        print self.evkin

        self.settings.enable()

        maxval = self.n_events * len(self.evkin.e_range)
        self._init_progress_bar(maxval)
        self.progressBar.start()
        for evkin in self.evkin.e_range:
            self.set_event_kinematics(evkin)
            rej_counter = 0
            for i in xrange(self.n_events):  # @UnusedVariable
                reject = self.generate_event()
                self.progressBar.update(self.progressBar.currval + 1)
                if reject:
                    rej_counter += 1
                    continue

                event = self.event_class(self.lib, self.event_config)
                fired_triggers = self.trigger_event(event)

                if not fired_triggers:
                    continue

                active_histograms = [
                    hist
                    for trigger in fired_triggers
                    for hist in trigger.histogram_list.values()
                ]
                [
                    hist.fill_event(event, evkin.ecm)
                    for hist in active_histograms
                ]

            # Finalize run
            sigmax = self.get_sigma_inel()
            sigma_gen = sigmax * (float(self.n_events - rej_counter)) / \
                float(self.n_events)

            for trigger in self.triggers.itervalues():
                for hist in trigger.histogram_list.itervalues():
                    hist.finalize_run(sigma_gen, evkin.ecm)
            print "...completed. " + str(100 * float(rej_counter) / float(self.n_events)) + \
                  "% of the events have been rejected."

        self.progressBar.finish()
        self.settings.reset()
        self.log_man.close_log()
        del self.lib