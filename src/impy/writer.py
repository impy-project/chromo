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
        import pyhepmc  # delay import until initialization

        self._writer = pyhepmc.WriterAscii(filename)
        self._genevent = pyhepmc.GenEvent()
        self._event_number = 0

    def write(self, event):
        self._fill_hepmc_event(event)
        self._writer.write_event(self._genevent)

    def close(self):
        self._writer.close()

    def _fill_hepmc_event(self, impy_event):
        # TODO:
        # - add cross-section info
        # - add info about generator

        pem = impy_event._pem_arr.T  # the need to transpose this is bad
        vt = impy_event._vt_arr.T  # the need to transpose this is bad
        n = pem.shape[0]

        breakpoint()

        # parents = impy_event.parents.T[:n]
        self._genevent.from_hepevent(
            event_number=self._event_number,
            momentum=pem[:, :4],
            mass=impy_event.m,
            position=vt,
            pid=impy_event.p_ids,
            parents=impy_event.parents.T[:n],
            children=impy_event.children.T[:n],
            status=impy_event.status,  # particle status
        )
        self._event_number += 1
