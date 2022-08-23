from abc import ABC, abstractmethod


class Writer(ABC):
    @abstractmethod
    def write(self, event):
        pass

    @abstractmethod
    def close(self):
        pass

    # support with statement
    def __enter__(self):
        return self

    # support with statement
    def __exit__(self, *args):
        self.close()


class HepMCWriter(Writer):
    def __init__(self, filename):
        import pyhepmc  # delay import till instantiation

        self._writer = pyhepmc.WriterAscii(filename)
        self._genevent = pyhepmc.GenEvent()

    def write(self, event):
        # TODO:
        # - add cross-section info
        # - add info about generator

        event.to_hepmc3(self._genevent)
        self._writer.write_event(self._genevent)

    def close(self):
        self._writer.close()
