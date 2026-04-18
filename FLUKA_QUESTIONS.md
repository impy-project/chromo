# FLUKA 2025.1 — open questions for the maintainer

I call FLUKA as a library: `STPXYZ` once, then `SGMXYZ`/`EVTXYZ` +
`FLLHEP`. Link is `ldpmqmd`-equivalent. FLRNDM for all random numbers.

---

## 1. Heavy-ion projectiles (A > 4) abort in EVTXYZ

**Problem.** `EVTXYZ` with `KPROJ = -2` (HEAVYION) after `PDGION`
aborts with `KPPRCT=-1/-2`. Light nuclei (`-3`..`-6`) work.
`SGMXYZ` works for all ions with the PDG-extended encoding.

**Reproduce.** Call `PDGION(1000060120)`, then
`EVTXYZ(-2, mat, ekin, ...)` at any energy.

**Question.** What is the correct call sequence for A > 4 projectiles
in EVTXYZ? Is `LIONTR`/`LHVTRN` needed before STPXYZ?

---

## 2. `e+/e-` projectiles: EVTXYZ aborts on some targets

**Problem.** `SGMXYZ` returns valid cross sections for e+/e- on all
targets. `EVTXYZ` works on light/medium targets (p through Ar40,
Cr52, Mn55, Cu63) but aborts with `IBTAR/ICHTAR ???` on Ca40, Fe56,
Pb208. Not a simple A cutoff.

**Reproduce.** `EVTXYZ(3, mat_Fe56, 10.0, ...)` → abort.
`EVTXYZ(3, mat_N14, 10.0, ...)` → OK.

**Question.** Which targets does EVTXYZ support for e+/e- projectiles?
Is this a known limitation or a missing initialization step?

---

## 3. UHE thresholds at 1e+30 after STPXYZ

**Problem.** `UHEHDT`, `UHEIOT`, `UHHDXT`, `UHIOXT` all read `1e+30`
after STPXYZ. The DPMJET→UHE transition never fires.

**Question.** Is this expected when no UHE model is linked? What
registers a UHE model at runtime?

---

## 4. GENTHR state after STPXYZ

After STPXYZ (PPTMAX=plab, EF2DP3=-1, DF2DP3=-1):

| variable | value  | variable | value |
|----------|--------|----------|-------|
| DPJHDT   | 20000  | PEANCT   | 100000|
| DPJIOT   | 12.5   | QMDIOT   | 0.125 |
| DPQMSM   | 2.0    | LAASMR   | 1     |
| FLDPSM   | 10000  | LFDSMR   | 0     |

**Question.** Is DPJHDT=20 TeV the intended Peanut→DPMJET threshold?
How does it interact with PEANCT=100 TeV?

---

## 5. EMD-only event generation (IFLXYZ=100) aborts

**Problem.** `SGMXYZ(…, 100)` returns valid EMD cross sections.
`EVTXYZ(…, 100)` aborts.

**Reproduce.** `EVTXYZ(kproj, mat, ekin, 0, 0, 0, 1, 100)`.

**Question.** Is EMD-only event generation supported, or must it be
combined with inelastic (IFLXYZ=101)?

---

## 6. Ranmar state completeness

**Problem.** I serialise Ranmar via RNINIT/RNWRIT/RNREAD and get
event-level reproducibility.

**Question.** Is Ranmar the only PRNG in FLUKA, or are there other
stochastic states that need separate serialisation?

---

## Resolved

**DPMJET event gen FPE above ~75 TeV CMS.** Was caused by hardcoded
`PPTMAX=1e9`. Fix: set `PPTMAX = plab` from construction kinematics.
100 TeV CMS now works.

**DPMJET event gen segfault.** Was caused by chromo's numpy RNG shims
(`rangen.f`) shadowing DPMJET's `DT_RNDM.o` (which tail-calls
`FLRNDM`). Fix: exclude `rangen` from the FLUKA build entirely.

---

## Build environment

- FLUKA 2025.1, DPMVERS=3.19.3.2
- `ldpmqmd`-equivalent link (libflukahp + libdpmmvax + librqmdmvax +
  libdpmjet3.19.3.2 + librqmd)
- gfortran (GCC toolset, Rocky Linux 9)
- FLRNDM for all random numbers
