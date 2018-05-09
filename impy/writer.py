from abc import ABCMeta, abstractmethod, abstractproperty

class Writer(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def write(): pass

    @abstractmethod
    def close(): pass

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


class HepMCWriter(Writer):
    import pyhepmc
    def __init__(self, filename):
        self._writer = self.pyhepmc.WriterAscii(filename)
        self._evt = self.pyhepmc.GenEvent()

    def write(self, event):
        self._fill_hepmc_event(self._evt, event)
        self._writer.write_event(self._evt)

    def close(self):
        self._writer.close()

    @classmethod
    def _fill_hepmc_event(cls, hepmc_event, impy_event):
        hepmc_event.clear()

        # TODO:
        # - need to add full particle history, not only final state
        # - add beam particles
        # - add cross-section info
        # - add info about generator
        make_particle = cls.pyhepmc.GenParticle
        for i in range(impy_event.px.shape[0]):
            p = make_particle((impy_event.px[i],
                               impy_event.py[i],
                               impy_event.pz[i],
                               impy_event.en[i]),
                              impy_event.p_ids[i],
                              3)
            hepmc_event.add_particle(p)