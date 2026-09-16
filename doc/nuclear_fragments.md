# Nuclear fragments and remnants

In hadron-nucleus and nucleus-nucleus collisions the generators run an
intranuclear cascade, and the leftover beam remnants stay in the event stack
next to the final-state hadrons. They carry status codes assigned by chromo
and are selected with `event.final_state_with_nucl_frag()`:

```python
event.final_state()                 # status == 1
event.final_state_with_nucl_frag()  # status in (1, 4, 5)
```

The return value is an `EventData` with the usual arrays (`pid`, `status`,
`charge`, kinematics, ...). Nucleus records carry the standard PDG code
`10LZZZAAAI`. The default `final_state()`, `final_state_charged()`, and the
HepMC3 export keep their existing behavior.

## Status assignment

chromo normalizes the generator stacks to three codes:

* **1** — final state particles, as returned by `final_state()`.
* **4** — nucleus records: the incoming beam particles (always records 0 and
  1, nuclear beams and targets get their `10LZZZAAAI` PDG code), plus
  residual nuclei if the generator writes them to the stack.
* **5** — nucleon remnants of the cascade: spectator, wounded, and
  potential-bound nucleons, as protons and neutrons.

Since remnants and beam records have status != 1, `final_state()` contains
only hadrons (issue #183).

The mapping from native codes per generator:

| generator | native codes → chromo status |
| --- | --- |
| DPMJET-III 1.9.3 / 3.0.7 | wounded 9/10/11/12/17/18, spectators 13/14, potential-bound 15/16 (`DT_COORDI`, `DT_RESNCL`, `DT_SCN4BA`) → 5; residual nuclei (PDG 80000, A/Z in `IDRES`/`IDXRES`) → 4 + PDG code |
| EPOS-LHC / EPOS-LHC(R) | cascade nucleons attached to the beam nucleus records → 5 |
| Pythia8Angantyr | spectators/excited beams 13/15/16 (11 for beam remnants) with nucleon PDG → 5; partonic remnants (21–73) and the diffractive Pomeron (990) keep their Pythia code |
| Pythia8Cascade | final state only, no remnant records |
| QGSJet (all), SIBYLL (all), UrQMD | beam records status 4 only; spectator production exists in the model (`qgfrgm`, `sibnuc`), the HEPEVT interface omits those records |

This matches what the interfaces expose: CORSIKA's DPMJET interface counts
the same nucleon records as projectile/target spectators, and CORSIKA's
SIBYLL interface keeps that generator's spectators in its own `/S_PLNUC/`
bookkeeping. A yield like $\sigma(p + C \to \mathrm{Be} + X)$ (issue #218)
therefore requires coalescence over the status 5 nucleons, since none of the
currently built generators writes a Be record.

The records with PDG ID 99999 in DPMJET events (issue #218) are the
hadronization chains of the dual parton model after fragmentation:
`DT_EVTFRG` overwrites their ID once the chain content is in the stack, and
their status codes (`jstrg = 100 * IPROCE + NCODE` in `DT_GETPJE`, e.g. 103,
506) encode the scattering process and string type. They are parton-level
bookkeeping; `final_state_with_nucl_frag()` discards them, the HepMC3 export
keeps them so the string history is complete.

## Caveats on the remnant records

The status 5 records are the nucleon bookkeeping of the cascade. The
residual nucleus, if reported, is the cascade output: the de-excitation
stage is missing or truncated depending on the generator (the DPMJET
steering card in chromo prints "No evaporation performed since evaporation
modules not available"), so A, Z, and kinematics of the remnants can differ
from physical fragments. Spectator nucleons carry Fermi motion, and wounded
nucleons carry their cascade kinematics. Use the records as bookkeeping for
what the cascade consumed, and model evaporation/fragmentation on top if
your observable needs fragment nuclei.

## Runnable example

[examples/extract_nuclear_fragments.py](../examples/extract_nuclear_fragments.py)
runs DPMJET-III 3.07 in p + Pb at 100 GeV/c:

```console
$ python examples/extract_nuclear_fragments.py   # generator banner omitted
status codes in the raw event stack:
  status   -22: 4 entries
  status   -21: 4 entries
  status     1: 14 entries
  status     2: 21 entries
  status     4: 2 entries
  status     5: 223 entries
  status   104: 1 entries
  status   106: 1 entries

final_state(): 14 hadrons
final_state_with_nucl_frag(): 239 records = 14 hadrons + 2 nucleus records + 223 remnant nucleons
nucleus records (status 4):
  pdgid=2212  (beam nucleon)
  pdgid=1000822080  A=208  Z=82
remnant nucleons (status 5): Counter({2112: 135, 2212: 88})
```

The beam Pb keeps its nominal A = 208 while 223 nucleons are recorded as
status 5, a consequence of the caveats above.
