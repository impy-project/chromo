# Nuclear fragments and remnants

Event generators model hadron-nucleus and nucleus-nucleus collisions by first
knocking nucleons out of the colliding nuclei (the intranuclear cascade) and
then breaking up and de-exciting what is left over. The leftover nuclear
fragments and beam remnants are part of the simulation, but they are *not*
returned by `event.final_state()`. This page explains why and shows how to get
them out of the raw event stack. It answers the questions behind
[issue #183](https://github.com/impy-project/chromo/issues/183) ("ions missing
from `final_state()`") and
[issue #218](https://github.com/impy-project/chromo/issues/218)
("how do I get the prefragment with PDG ID 99999").

## Status codes: raw stack vs. final state

Every event is a flat stack of particle records, and each record has a status
code. `EventData.final_state()` and `final_state_charged()` select the records
with `status == 1`, i.e. long-lived particles that are the result of the
interaction and everything that decayed into them. This is by design: nuclear
fragments and beam remnants are not final-state particles in that sense, they
carry other status codes, so `final_state()` removes them together with the
parton-level history. If you need them, select them from the *raw* event with
numpy masks over `event.status`, `event.pid`, `event.charge`, and `event.m`.

chromo itself normalizes one code across all generators: the incoming beam
particles are stored as records 0 and 1 with `status == 4`, and for nuclear
beams or targets these are the beam nuclei. In models that keep the cascade
remnants in the stack (EPOS-LHC, DPMJET), status 4 also marks the nucleons and
nucleus fragments that belong to the beam remnants. For EPOS-LHC this is all
there is; DPMJET additionally hands its de-excitation model explicit
prefragment records, see below. All other status codes are model-specific.

## DPMJET prefragments: PDG ID 99999 and status codes >= 100

DPMJET (DpmjetIII193, DpmjetIII307, Phojet) writes the excited nuclear
fragments that its evaporation/fission/fragmentization model has to de-excite
into the event stack as placeholder entries with PDG ID **99999** and status
codes **103, 104, 106** (for DPMJET-III 3.07, higher codes in the 5xx/6xx/7xx
families appear as well, e.g. 503, 506; the pattern so far is always
`(status >= 100) and (pid == 99999)`). These are the records a
$\sigma(p + C \to B + X)$ study is after. The PDG ID and status codes are
written by the Fortran, see `DT_EVTFRG` in the bundled sources; note that
`DpmjetIIIEvent._prepare_for_hepmc` in `src/chromo/models/dpmjetIII.py` uses
exactly the mask `status in (1, 2, 4)` or `pid == 99999` to keep this
information when converting to HepMC3.

Two caveats when you select prefragments:

- In DPMJET-III 3.07, entries with `pid == 99999` and `status == 2` are *not*
  nuclear fragments. They are color-connected hadronization chains which the
  generator has already converted into the final-state hadrons. The nuclear
  prefragments are the ones with `status >= 100`. (Both kinds of placeholders
  are kept when exporting to HepMC3, so check twice if you count them there.)
- The 99999 entry is the *prefragment before* de-excitation, so it has no
  charge and no mass/energy bookkeeping for nucleons. Its rest mass is large
  (GeV-scale) but the mass/charge/energy of the physical fragment you measure
  (a Be from C, a C from O, ...) must be reconstructed from the stable decay
  products in the final state, which DPMJET does not group under a nucleus
  PDG code in the chromo stack. Only the incoming beam nucleus is stored with
  a proper nucleus PDG code (10LZZZAAAI) and `status == 4`.

## QGSJet: no fragment information

The QGSJet converters (`chepevt` in the bundled Fortran) mark every generated
particle with status 1 and chromo prepends the beam records with status 4, so
the only ion in the raw stack of a QGSJet h+A run is the incoming nucleus
(`status == 4`); `final_state()` contains no nuclei at all and there is no
prefragment record to select. QGSJet-III simply does not write the fragments
that its `hardq/frag` stage produces into the HEPEVT stack.

## Runnable example

[examples/extract_nuclear_fragments.py](../examples/extract_nuclear_fragments.py)
runs DPMJET-III 3.07 in p + O at 1e5 GeV/c and prints the observed status
codes, the beam records, and the nuclear prefragments:

```console
$ python examples/extract_nuclear_fragments.py   # generator banner omitted
status codes in the raw event stack:
  status   -42: 2 entries
  status   -22: 2 entries
  status   -21: 4 entries
  status     1: 25 entries
  status     2: 15 entries
  status     4: 2 entries
  status     9: 1 entries
  status    10: 1 entries
  status    12: 1 entries
  status    14: 14 entries
  status   503: 1 entries
  status   506: 1 entries

beam records (indices 0, 1):
  index 0: pdgid=2212  status=4
  index 1: pdgid=1000080160  status=4

nuclear prefragments (pid 99999, status >= 100): 2
  status=506  m=11.8 GeV  E=175.8 GeV
  status=503  m=10.3 GeV  E=41.6 GeV

final_state(): 25 particles, nuclei among them: 0
```

The selection logic of the example, applied to any event `ev` of a DPMJET
h+A run:

```python
import numpy as np

is_prefragment = (ev.pid == 99999) & (ev.status >= 100)   # excited nuclear fragments
is_beam_remnant = (ev.status == 4) & (np.abs(ev.pid) > 1000000000)  # beam nuclei
```
