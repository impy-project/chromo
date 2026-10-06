# Nuclear fragments and remnants

`EventData.final_state()` returns status 1 records only. In hadron-nucleus and
nucleus-nucleus collisions, nuclear remnants are selected with

```python
event.final_state_with_nucl_frag()
```

which returns terminal records with status 1, 4, or 5. Incoming beam records
(indices 0 and 1) and records with daughters are excluded, so sums of baryon
number or charge over the selection do not double-count the initial state.
`final_state()`, `final_state_charged()`, and the HepMC3 export are unchanged.

## Status codes

| status | content |
| --- | --- |
| 1 | final state particles |
| 4 | residual nuclei, PDG code `10LZZZAAAI` |
| 5 | spectator nucleons of the intranuclear cascade |

## Generator support

| generator | remnant records |
| --- | --- |
| DpmjetIII | residual-nucleus records (native 13–16) → 5; residual nuclei (PDG 80000 with `IDRES`/`IDXRES`) → 4; wounded nucleons (9–12, 17, 18) keep native codes |
| EposLHC, EposLHCR | fragments at status 1 |
| Pythia8Angantyr | residual nuclei at status 1; isomer digit `I=9` reset to `0` |
| Pythia8Cascade | none (final state only) |
| QGSJet, SIBYLL, UrQMD | none |

## Caveats

- DPMJET is built without evaporation. Its residual nuclei are given as
  status 5 nucleons; baryon number and charge of the selection are exact.
- Status 5 nucleons carry cascade kinematics; light fragments (e.g. Be)
  require coalescence.

## Example

[examples/extract_nuclear_fragments.py](../examples/extract_nuclear_fragments.py)
selects the remnants of p + Pb events with DPMJET-III 3.0.7.
