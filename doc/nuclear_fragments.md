# Nuclear fragments and remnants

In hadron-nucleus and nucleus-nucleus collisions the generators run an
intranuclear cascade, and the leftover beam remnants stay in the event stack
next to the final-state hadrons. They carry status codes assigned by chromo
and are selected with `event.final_state_with_nucl_frag()`:

```python
event.final_state()                 # status == 1
event.final_state_with_nucl_frag()  # physical final state + remnants
```

The returned event contains physical, terminal records only, and the
incoming beam records (always raw indices 0 and 1 with status 4) are
excluded. This means summing baryon number, charge, or energy over the
selection cannot double-count the initial state. The selection is a superset
of `final_state()`, and the default filters and the HepMC3 export keep their
existing behavior.

## Status assignment

chromo normalizes the generator stacks as follows:

* **1** — final state particles, as returned by `final_state()`. Some
  generators already write their nuclear fragments here (see table).
* **4** — residual/fragment nucleus records reported by the generator, with
  the standard PDG code `10LZZZAAAI`. Nucleon records and records with
  daughters are filtered out of the remnant selection.
* **5** — terminal nucleon remnants of the cascade: spectator and
  potential-bound nucleons that left the cascade unhit, at Fermi momentum.
  Never any record with daughters.

The mapping from native codes per generator, and what it implies for the
remnant selection:

| generator | native codes → chromo status | remnant content |
| --- | --- | --- |
| DpmjetIII 1.9.3 / 3.0.7 | terminal nucleons of 13/14/15/16 → 5; wounded 9–12/17/18 keep native code; residual nuclei (PDG 80000 with `IDRES`/`IDXRES`) → 4 + PDG code | status 5 nucleons; no residual nuclei (evaporation disabled, see below) |
| EposLHC / EposLHCR | nothing relabelled — nuclear fragments are already written at status 1 | status 1 fragments; selection = `final_state()` |
| Pythia8Angantyr | residual nuclei ("NucRem") are appended at status 1 with the isomer digit I=9; chromo rewrites them to a standard PDG code (I=0) | status 1 nuclei; selection = `final_state()` |
| Pythia8Cascade | final state only (`listFinalOnly`) | none |
| QGSJet (all), SIBYLL (all), UrQMD | beam records status 4 only; spectator production exists in the model (`qgfrgm`, `sibnuc`), the HEPEVT interface omits those records | none; selection = `final_state()` |

Conservation has been checked event-by-event over
`final_state_with_nucl_frag()` for all built models (see
`tests/test_common.py::test_final_state_with_frags` and
`tests/test_dpmjetIII.py::test_baryon_conservation`). EPOS-LHC, UrQMD, and
Pythia8Angantyr close exactly for baryon number and charge. DPMJET closes up
to the wounded nucleons absorbed into the residual nucleus: its de-excitation
stage is switched off in the chromo build (the banner prints "No evaporation
performed since evaporation modules not available"), so those nucleons are
not reported by any record; the DPMJET manual itself defines the event final
state as its status 1 only. The wounded nucleon records (native codes
9–12/17/18) are kept with their native status so they can be found in the
raw stack; they are intermediate records — they have daughters — and are
deliberately excluded from the remnant selection.

The records with PDG ID 99999 in DPMJET events (issue #218) are the
hadronization chains of the dual parton model after fragmentation:
`DT_EVTFRG` overwrites their ID once the chain content is in the stack, and
their status codes (`jstrg = 100 * IPROCE + NCODE` in `DT_GETPJE`, e.g. 103,
506) encode the scattering process and string type. They are parton-level
bookkeeping; `final_state_with_nucl_frag()` discards them, the HepMC3 export
keeps them so the string history is complete.

## Caveats

Status 5 kinematics are those of the cascade: nucleons are at rest in the
residual-nucleus frame up to Fermi motion and are not boosted, so
four-momentum closure over the remnant selection requires transforming the
remnants to the event frame. Where the generator's de-excitation stage is
skipped or truncated, A, Z, and kinematics of the remnants can differ from
physical fragments: DPMJET keeps the beam nucleus records nominal (A = 208
for Pb) while its nucleons are listed at status 5, and the fragment
distributions that evaporation or fission would produce are absent. For a
fragment yield like $\sigma(p + C \to \mathrm{Be} + X)$ (issue #218),
none of the currently built generators writes a Be record; it must be
reconstructed by coalescence over the status 5 nucleons (or the status 1
fragments of EPOS/Angantyr/UrQMD).

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
  status     5: 206 entries
  status     9: 1 entries
  status    10: 1 entries
  status    12: 1 entries
  status    18: 14 entries
  status   104: 1 entries
  status   106: 1 entries

final_state(): 14 hadrons
final_state_with_nucl_frag(): 220 records = 14 hadrons + 0 nucleus records + 206 remnant nucleons
nucleus records (status 4):
remnant nucleons (status 5): Counter({2112: 123, 2212: 83})
```

The wounded records with native codes and the status 4 beam records are
visible in the raw stack and absent from the remnant selection, which is the
point: the 123 n + 83 p are the physical breakup products of the target.
