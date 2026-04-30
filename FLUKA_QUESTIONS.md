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

## 4. `QRDDCY` is in `dcytst.f` rather than `libflukahp.a`

**Problem.** `QRDDCY(IADCYP, IZDCYP, ISDCYP, IFLDCY, LNCMSS)` — Q-value
for a radioactive-decay channel — is shipped as user-supplied source in
the `dcytst.f` test harness, not as a library symbol. `nm libflukahp.a |
grep qrddcy` is empty across all archives at `$FLUPRO`.

**Workaround.** chromo's `chromo_fluka.f` ships a verbatim copy of
`QRDDCY` (vintage FLUKA 2025.1, last upstream change 25-Apr-26) so the
chromo `_fluka` extension can link `chromo_dcy_channels`. This
introduces a drift risk on FLUKA bumps; the embed has a
`TODO(FLUKA-bump)` marker to flag re-sync.

**Question.** Should `QRDDCY` be exported by a future FLUKA library
release, or is it intended to remain user-supplied? If user-supplied,
is there a recommended location (vendored, or sourced from a stable
auxiliary archive)?

---

## 3. EMD-only event generation (IFLXYZ=100) aborts

**Problem.** `SGMXYZ(…, 100)` returns valid EMD cross sections.
`EVTXYZ(…, 100)` aborts.

**Reproduce.** `EVTXYZ(kproj, mat, ekin, 0, 0, 0, 1, 100)`.

**Question.** Is EMD-only event generation supported, or must it be
combined with inelastic (IFLXYZ=101)?

---

## Build environment

- FLUKA 2025.1, DPMVERS=3.19.3.2
- `ldpmqmd`-equivalent link (libflukahp + libdpmmvax + librqmdmvax +
  libdpmjet3.19.3.2 + librqmd)
- gfortran (GCC toolset, Rocky Linux 9)
- FLRNDM for all random numbers
