from abc import ABCMeta, abstractmethod, abstractproperty
from six import with_metaclass
import numpy as np


class Writer(object, with_metaclass(ABCMeta)):
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
        self._hep = __import__('pyhepmc_ng') # delay import till instantiation
        self._writer = self._hep.WriterAscii(filename)
        self._genevent = self._hep.GenEvent()
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

        self._genevent.clear()
        pem = impy_event.pem_arr.T # the need to transpose this is bad
        vt = impy_event.vt_arr.T # the need to transpose this is bad
        n = pem.shape[0]
        if impy_event.parents is None:
            parents = 0
        else:
            parents = impy_event.parents.T[:n]
        self._hep.fill_genevent_from_hepevt(
            self._genevent,
            event_number=self._event_number,
            p=pem[:, :4],
            m=impy_event.m,
            v=vt,
            pid=impy_event.p_ids,
            parents=parents,
            status=impy_event.status, # particle status
        )
        self._event_number += 1
