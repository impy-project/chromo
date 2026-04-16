# FLUKA 2025.1 — questions on the STPXYZ / SGMXYZ / EVTXYZ API

I am calling FLUKA as a library through a thin Fortran wrapper around
`STPXYZ`, `SGMXYZ`, `EVTXYZ`, and `FLLHEP`. The link is equivalent to
what `flutil/ldpmqmd` produces (`libflukahp.a` + `libdpmmvax.a` +
`librqmdmvax.a` + `libdpmjet3.19.3.2.a` + `librqmd.a`). FLRNDM is used
for all random numbers (no external RNG override).

Call sequence: `STPXYZ` once at init (materials, `PPTMAX=1e9`,
`EF2DP3=-1`, `DF2DP3=-1`), then repeated `SGMXYZ`/`EVTXYZ` + `FLLHEP`.

---

## 1. EVTXYZ with heavy-ion projectile (A > 4): `KPPRCT=-1/-2`

Light nuclei work when passed with their dedicated PAPROP codes
(`-3` d, `-4` t, `-5` ³He, `-6` ⁴He) — event generation through rQMD
and DPMJET succeeds at all tested energies.

For heavier projectiles (e.g. C12, O16) I call `PDGION(1000060120)`
before `EVTXYZ` and pass `KPROJ = -2` (HEAVYION). `SGMXYZ` with the
same PDG-encoded `KPROJ` returns sensible cross sections (~3.6 b for
C+Pb). But `EVTXYZ` aborts:

```
STOP KPPRCT=-1/-2
STOP STOP: FLUKA ABORTED
```

This happens at every energy (1 GeV/n through 10 TeV/n). Without the
prior `PDGION` call, we get `STOP DUMMY EVEUHE CALLED` instead, so
`PDGION` does register something — but `EVTXYZ` still rejects the ion.

After `STPXYZ` the heavy-ion transport flags in `(THRSCM)` read
`LIONTR = 0`, `LHVTRN = 0`, `KHVTRN = 0`.

**Questions.**
- What is the correct call sequence to make `EVTXYZ` accept an A > 4
  projectile? Is `PDGION` + `KPROJ = -2` correct, or is `SETION`
  needed, or must `LIONTR`/`LHVTRN` be set before `STPXYZ`?
- Is AA event generation through `EVTXYZ` only available when an
  external UHE model (EPOS/SIBYLL/QGSJET) is linked? The `DUMMY
  EVEUHE CALLED` message suggests heavy-ion projectiles always
  dispatch to the UHE entry.

---

## 2. FPE in EVTXYZ via DPMJET above ~75 TeV CMS

`EVTXYZ` for hA (e.g. π⁺ + Pb) raises a floating-point exception when
the CMS energy exceeds roughly 75 TeV. `SGMXYZ` returns valid cross
sections at all energies. Standalone DPMJET-III 19.3.2 generates events
without issues above this energy — the FPE is specific to the FLUKA
adapter path.

Bisect with π⁺ + Pb208:

| E_lab (GeV)   | √s (TeV) | result |
|---------------|-----------|--------|
| 3.0 × 10⁹    | ~75       | OK     |
| 3.1 × 10⁹    | ~76       | FPE    |
| 5.3 × 10⁹    | 100       | FPE    |

Call chain at the FPE: `EVTXYZ → EVENTP → EVENTD → DPMRUN →
DT_KKINC → DT_KKEVNT → DT_GLAUBE → DT_DIAGR`.

`PPTMAX` was set to `1e9 GeV/c` in `STPXYZ`. `fort.11` reports
`Max. initialization energy for DPMJET 1.8e9 GeV`.

**Questions.**
- Is ~75 TeV CMS the expected upper limit of DPMJET-3.19 event
  generation inside FLUKA?
- Would raising `PPTMAX` in the `STPXYZ` call extend this range?

---

## 3. `PPCMSX` and `FLKIOM` zero after STPXYZ

After `STPXYZ`, the `(GENTHR)` variables `PPCMSX` ("maximum CMS
momentum for pp interactions, used for DPMJET initialization") and
`FLKIOM` ("FLUKA ion interaction max energy per nucleon") both read
`0.0`, despite `PPTMAX = 1e9` being passed.

**Question.** Are these set later (first `EVTXYZ` call), or only via
an input card? Could `PPCMSX = 0` be related to the FPE in entry #2?

---

## 4. `UHEHDT` / `UHEIOT` at 1e+30 after STPXYZ

All four UHE thresholds (`UHEHDT`, `UHEIOT`, `UHHDXT`, `UHIOXT`) read
`1e+30` after `STPXYZ`, meaning the DPMJET→UHE transition never fires.

**Question.** Is this the expected default when no UHE model is linked,
or should `STPXYZ` set these based on `LUHEEX`? What registers a UHE
model at runtime?

---

## 5. `GENTHR` state after STPXYZ

Dump of `(GENTHR)` after `STPXYZ` (PPTMAX=1e9, EF2DP3=-1, DF2DP3=-1):

| variable | value   | notes                                |
|----------|---------|--------------------------------------|
| DPJHDT   | 20000   | hA Peanut→DPMJET threshold           |
| FLDPSM   | 10000   | smearing (but LFDSMR=0 → disabled)   |
| DPJIOT   | 12.5    | AA rQMD→DPMJET threshold             |
| QMDIOT   | 0.125   | rQMD lower bound                     |
| DPQMSM   | 2.0     | AA smearing                          |
| LAASMR   | 1       |                                      |
| PEANCT   | 100000  | per-particle Peanut thresholds (all) |
| UHEHDT   | 1e+30   | see entry #4                         |
| PPCMSX   | 0       | see entry #3                         |
| FLKIOM   | 0       | see entry #3                         |

**Question.** `DPJHDT = 20 TeV` while the static `PARAMETER RQDP31 =
12.5` — is 20 TeV the intended Peanut→DPMJET threshold for all hadrons?
How does `DPJHDT` interact with `PEANCT` (100 TeV)?

---

## 6. STPXYZ re-entry

Calling `STPXYZ` a second time triggers `FLABRT`. I register all
materials in a single `STPXYZ` call, but would like to add materials
later (e.g. a target nucleus not in the initial list).

**Question.** Is there a way to register additional materials after the
first `STPXYZ` — via `SETITB`/`DFATWG` directly, or a re-entry flag?

---

## 7. `EVTXYZ` with `IFLXYZ = 100` (EMD-only)

`SGMXYZ(…, IFLXYZ=100)` returns valid EMD cross sections (e.g.
O16+Pb208 at 1.6 TeV). `EVTXYZ(…, IFLXYZ=100)` aborts. I work around
this with `IFLXYZ = 101` (inelastic + EMD).

**Question.** Is EMD-only event generation via `EVTXYZ` supported, or
must it always be combined with inelastic?

---

## 8. `EVTXYZ` with electron/positron projectile

Calling `EVTXYZ` with the electron PAPROP code aborts.

**Question.** Does `EVTXYZ` support leptonic projectiles, or should
users convert to equivalent-photon kinematics upstream?

---

## 9. Ranmar state completeness

I serialise Ranmar state via `RNINIT` / `RNWRIT` / `RNREAD` and
confirm event-level reproducibility across processes.

**Question.** Is Ranmar the only PRNG in FLUKA, or are there additional
stochastic states (Glauber tables, evaporation, …) that need separate
serialisation for full reproducibility?

---

## Build environment

- FLUKA 2025.1, `DPMVERS = 3.19.3.2`
- Archives linked: `ldpmqmd`-equivalent (all five: `libflukahp`,
  `libdpmmvax`, `librqmdmvax`, `libdpmjet3.19.3.2`, `librqmd`)
- gfortran (GCC toolset, Rocky Linux 9)
- FLRNDM used for all random numbers (no external RNG replacement)
