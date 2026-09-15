# Nuclear fragments and remnants

Event generators model hadron-nucleus and nucleus-nucleus collisions by first
knocking nucleons out of the colliding nuclei (the intranuclear cascade) and
then breaking up and de-exciting what is left over. The leftover beam
remnants are part of the simulation, but they are *not* returned by
`event.final_state()`. This page explains what the event stack contains, how
chromo normalizes it, and — importantly — what is *missing* from it. It
answers the questions behind
[issue #183](https://github.com/impy-project/chromo/issues/183) ("ions missing
from `final_state()`") and
[issue #218](https://github.com/impy-project/chromo/issues/218)
("what is the particle with PDG ID 99999?").

## The chromo convention: status 1, 4, 5

Every event is a flat stack of particle records, and each record has a status
code. chromo normalizes the generator-specific stacks of all hadronic models
to a common convention:

* **status 1** — final state particles: long-lived hadrons and leptons that
  are the result of the interaction (and everything that decayed into them).
  `final_state()` and `final_state_charged()` select these.
* **status 4** — nucleus records: the incoming beam particles (always records
  0 and 1), plus residual/fragment nuclei when the generator reports them.
  Nuclear records carry the standard PDG code `10LZZZAAAI` (e.g. `1000822080`
  for a lead residual).
* **status 5** — nucleon-level remnants: spectator, wounded, and
  potential-bound nucleons left over from the intranuclear cascade, as
  protons and neutrons.

The remnant records never have status 1, so they do not appear in
`final_state()` — this is why it contains no ions (issue #183). Use the
companion filter to get them:

```python
event.final_state_with_nucl_frag()   # status in (1, 4, 5)
```

It returns a new `EventData` with the same API; the default `final_state()`
and the HepMC3 export are unaffected.

!!! warning "Remnant records may not be physical"
    Status 4/5 records are *bookkeeping*, not a complete fragmentation model.
    Generators typically record the cascade participants and leave the
    residual nucleus as it entered the cascade, and several skip the
    de-excitation stage entirely (DPMJET's default steering card in chromo
    prints "No evaporation performed since evaporation modules not
    available"). The remnants are therefore often missing evaporation
    fragments, fission, and nuclear breaks-up, and their A/Z/kinematics can
    be wrong or merely nominal. Treat them as "what the cascade used up",
    not as predicted fragment yields. The generator-specific caveats below
    are stricter in places.

## DPMJET: nucleon bookkeeping; PDG ID 99999 are *not* fragments

The DPMJET raw stack records the intranuclear cascade on the nucleon level:
wounded participant nucleons (native codes 9/10/11/12, 17/18 when
re-scattered), spectators (13/14), and nucleons bound in the nuclear potential
or their excitations (15/16); see `DT_COORDI`, `DT_RESNCL`, `DT_SCN4BA`.
chromo maps all of these to status 5. CORSIKA's DPMJET interface (`DPMJST`)
counts exactly the same records as projectile and target spectators, so this
bookkeeping is the community-recognized "remnant" content.

The residual-nucleus records of the DPMJET steering-card scheme (PDG ID 80000
with A and Z in the extended-history slots `IDRES`/`IDXRES`) are *not*
produced in chromo's default configuration — not once in 200 p+Pb events from
10 GeV to 100 TeV, in both 1.9.3 and 3.0.7. If they appear (non-default
options), chromo converts them to proper PDG nucleus codes with status 4.
Until then, a yield like $\sigma(p + C \to \mathrm{Be} + X)$ (issue #218) has
to be reconstructed by coalescing final-state nucleons (a Be is the
momentum-sum of 4 p + 4 n inside a narrow kinematic window); there is no Be
record to select.

The records with PDG ID **99999** that you may have noticed are *not* nuclear
fragments. They are the color-neutral hadronization chains (strings) of the
model *after* their content was converted into final-state hadrons: the
Fortran writes them with a status code encoding the scattering process and
the string type (`DT_GETPJE`, `jstrg = 100 * IPROCE + NCODE`), and
`DT_EVTFRG` resets the PDG ID to 99999 once the chain is fragmented (the
native 103/104/106 or 3.07's 5xx/6xx/7xx families). Status 2 chains were too
small to fragment and collapsed into a single hadron. 99999 is parton-level
bookkeeping; the HepMC3 export keeps those records (see
`DpmjetIIIEvent._prepare_for_hepmc`) only so the string history is not
orphaned. `final_state_with_nucl_frag()` discards them.

## EPOS-LHC: nucleon-level remnants

EPOS-LHC keeps the remnants on the nucleon level, too: next to the incoming
beam records, the raw stack contains cascade nucleons attached as daughters
of the beam nucleus records; chromo demotes them from the generator's status 4
to the remnant status 5 (see `_repair_initial_beam` in
`src/chromo/models/epos.py`). A p+O event at 100 TeV typically shows a few
dozen; there are no fragment nuclei with A > 1 besides the incoming beam
record.

## Pythia8 / Pythia8Angantyr: nucleon remnants, colored parton leftovers

The Angantyr heavy-ion mode reports spectators and excited beam nucleons with
Pythia status codes 13/15/16 (nucleons; 11 in some configurations). chromo
maps nucleon records among those codes to status 5; nucleus records get
status 4. The incoming beam nucleus is a status 4 record with a proper PDG
code. Be aware that the Angantyr stack also keeps *partonic* beam remnants at
codes 21–73 (colored records) and a diffractive Pomeron placeholder (PDG 990,
status 13); neither is a physical remnant, and both are excluded by
`final_state_with_nucl_frag()`. After the hadronization stage Pythia8 has no
residual nuclei — a Glauber picture of the collision, not a fragmentation
model. `Pythia8Cascade` exports final-state particles only and shows no
remnant records at all.

## QGSJet and SIBYLL: no fragment records in the stack

The QGSJet converters (`chepevt` in the bundled Fortran) mark every generated
particle with status 1 and chromo prepends the beam records with status 4, so
the only ion in the raw stack of a QGSJet h+A run is the incoming nucleus and
`final_state()` contains no nuclei at all. QGSJet-III does not write the
spectator fragments produced by its `qgfrgm` routine into the HEPEVT stack.
SIBYLL behaves the same in chromo (`sibnuc`/`sibhep` expose no A > 1 records);
when run under CORSIKA, its spectator fragments are kept separately in
CORSIKA's own `/S_PLNUC/` bookkeeping with code `1000 + A` (see `SIBNUC` in
CORSIKA's `sibyll2.3e.f`), not in the particle stack. UrQMD likewise exposes
only the beam records and final-state hadrons. For these generators
`final_state_with_nucl_frag()` just adds the incoming beam records (status 4)
to the final state; the nucleon-remnant status 5 does not appear.

## Summary of what each generator reports

| generator | status 4 (nucleus records) | status 5 (remnant nucleons) | A > 1 residual nuclei |
| --- | --- | --- | --- |
| DpmjetIII 1.9.3 / 3.0.7 | incoming beam only | yes (full cascade bookkeeping) | no (default config) |
| EposLHC / EposLHCR | incoming beam only | yes | no |
| Pythia8Angantyr | incoming beam only | yes (spectators/excited) | no |
| Pythia8Cascade | none (final state only) | no | no |
| QGSJet / SIBYLL / UrQMD | incoming beam only | no | no |

## Runnable example

[examples/extract_nuclear_fragments.py](../examples/extract_nuclear_fragments.py)
runs DPMJET-III 3.07 in p + Pb at 100 GeV/c and prints the normalized event
content:

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

chain placeholders (pid 99999): 2,
  in final_state_with_nucl_frag(): 0
```

Note the incoming beam nucleus still carries its full nominal A = 208: the
~223 status-5 nucleons were removed from it by the cascade (the Pb is left as
a broken-up pile of nucleons, per the warning above), and the two p+Pb chains
(99999) are correctly not part of the remnant selection.
