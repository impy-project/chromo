# Nuclear fragments and remnants

Event generators model hadron-nucleus and nucleus-nucleus collisions by first
knocking nucleons out of the colliding nuclei (the intranuclear cascade) and
then breaking up and de-exciting what is left over. The leftover beam
remnants are part of the simulation, but they are *not* returned by
`event.final_state()`. This page explains why and shows what is actually in
the raw event stack. It answers the questions behind
[issue #183](https://github.com/impy-project/chromo/issues/183) ("ions missing
from `final_state()`") and
[issue #218](https://github.com/impy-project/chromo/issues/218)
("what is the particle with PDG ID 99999?").

## Status codes: raw stack vs. final state

Every event is a flat stack of particle records, and each record has a status
code. `EventData.final_state()` and `final_state_charged()` select the records
with `status == 1`, i.e. long-lived particles that are the result of the
interaction and everything that decayed into them. This is by design: the
event stack also contains incoming-beam records and parton-level history, and
those carry other status codes, so `final_state()` removes them. If you need
them, select them from the *raw* event with numpy masks over `event.status`,
`event.pid`, `event.charge`, and `event.m`.

chromo normalizes one code across all generators: the incoming beam particles
are stored as records 0 and 1 with `status == 4`, so for a nuclear beam or
target the incoming nucleus is in the raw stack as a nucleus record with a
proper PDG code (10LZZZAAAI) but with status 4, never status 1 — this is why
`final_state()` contains no ions (issue #183). The generator-specific remnant
conventions are described below.

## DPMJET: PDG ID 99999 marks hadronization chains, not nuclear fragments

In DpmjetIII193 and DpmjetIII307, the records with PDG ID **99999** are the
color-neutral hadronization chains (strings) of the underlying model *after*
their content has been converted into final-state hadrons: the Fortran writes
the string records with a status code encoding the scattering process and the
string type (see `DT_GETPJE`, `jstrg = 100 * IPROCE + NCODE`), and
`DT_EVTFRG` resets the PDG ID to 99999 once the chain is fragmented (this is
also the code 103/104/106 you may know from native DPMJET output; 3.07
additionally shows the 5xx/6xx/7xx families, e.g. 503, 506). Entries with
`status == 2` are chains that were too small to fragment and collapsed into a
single hadron. So 99999 is parton-level bookkeeping — counting or histogramming
these records answers nothing about nuclear fragments, and the HepMC3 export
keeps them (see `DpmjetIIIEvent._prepare_for_hepmc`) only so the string history
is not orphaned.

What *is* in the raw stack is the nucleon-level bookkeeping of the
intranuclear cascade: wounded participant nucleons (statuses 9/10/11/12,
17/18 when re-scattered) and spectator nucleons (13/14), with 15/16 marking
nucleons bound in the nuclear potential and their excitations, see
`DT_COORDI`, `DT_RESNCL`, `DT_SCN4BA`. CORSIKA's DPMJET interface (`DPMJST`)
counts exactly these records as projectile and target spectators. There are,
however, no nucleus records with A > 1 besides the incoming beam record: the
residual-nucleus/fragment records of the steering-card scheme (PDG ID 80000
with A and Z in the extended history `IDRES`/`IDXRES`, statuses -1/1001) are
not produced in chromo's default configuration (not once in 200 p+Pb events
from 10 GeV to 100 TeV), so a yield like $\sigma(p + C \to \mathrm{Be} + X)$
(issue #218) must be reconstructed by combining final-state nucleons (a Be is
the momentum-sum of 4 p + 4 n within a narrow kinematic window); there is no
Be record to select.

## EPOS-LHC: nucleon-level remnants with status 4

EPOS-LHC keeps the remnants on the nucleon level. Next to the incoming beam
records, the raw stack contains cascade nucleons (protons and neutrons,
`status == 4`) attached as daughters of the nucleus records, see
`_repair_initial_beam` in `src/chromo/models/epos.py`. These `status == 4`
nucleons are the beam-remnant content you can analyze; like in DPMJET, there
are no fragment nuclei with A > 1 besides the incoming beam record.

## QGSJet and SIBYLL: no fragment records in the stack

The QGSJet converters (`chepevt` in the bundled Fortran) mark every generated
particle with status 1 and chromo prepends the beam records with status 4, so
the only ion in the raw stack of a QGSJet h+A run is the incoming nucleus and
`final_state()` contains no nuclei at all. QGSJet-III does not write the
spectator fragments produced by its `qgfrgm` routine into the HEPEVT stack.
SIBYLL behaves the same in chromo (`sibnuc`/`sibhep` expose no A > 1 records);
when run under CORSIKA, its spectator fragments are kept separately in
CORSIKA's own `/S_PLNUC/` bookkeeping with code `1000 + A` (see `SIBNUC` in
CORSIKA's `sibyll2.3e.f`), not in the particle stack.

## Runnable example

[examples/extract_nuclear_fragments.py](../examples/extract_nuclear_fragments.py)
runs DPMJET-III 3.07 in p + O at 1e5 GeV/c and prints the observed status
codes, the beam records, and the chain placeholders:

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

chain placeholders (pid 99999): 3 total,
  1 collapsed to single hadrons (status 2),
  2 fragmented strings (statuses [503, 506])

final_state(): 25 particles, nuclei among them: 0
```
