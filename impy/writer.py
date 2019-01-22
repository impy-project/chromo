from abc import ABCMeta, abstractmethod, abstractproperty

class Writer(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def write(self, event): pass

    @abstractmethod
    def close(self): pass

    # support with statement
    def __enter__(self):
        return self

    # support with statement
    def __exit__(self, *args):
        self.close()


class HepMCWriter(Writer):
    def __init__(self, filename):
        self.hep = __import__('pyhepmc_ng') # delay import till instantiation
        self._writer = self.hep.WriterAscii(filename)
        self._evt = self.hep.GenEvent()

    def write(self, event):
        self._fill_hepmc_event(self.hep, self._evt, event)
        self._writer.write_event(self._evt)

    def close(self):
        self._writer.close()

    @staticmethod
    def _fill_hepmc_event(hep, hepmc_event, impy_event):
        hepmc_event.clear()
        # TODO:
        # - need to add full particle history, not only final state
        # - add beam particles
        # - add cross-section info
        # - add info about generator
        make_particle = hep.GenParticle
        for i in range(impy_event.px.shape[0]):
            p = make_particle((impy_event.px[i],
                               impy_event.py[i],
                               impy_event.pz[i],
                               impy_event.en[i]),
                              impy_event.p_ids[i],
                              3)
            hepmc_event.add_particle(p)