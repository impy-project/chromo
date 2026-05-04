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

## 5. `SPDCEV` returns `LSUCCS=.FALSE.` for Ac-228 (and possibly others)

**Problem.** Calling `SPDCEV(228, 89, 0, ...)` returns `LSUCCS=.FALSE.`
even though Ac-228 has decay data in FLUKA's table (β⁻ to Th-228, T1/2
≈ 6.15 h). The Th-232 natural-decay series therefore terminates
prematurely at Ac-228 instead of running through to Pb-208.

**Reproduce.** After standard decay-table init (`NCDTRD`/`RDFLUO`/...),
call `SPDCEV(228, 89, 0, ZERZER, ZERZER, ZERZER, ZERZER, ONEONE,
ONEONE, -1.D9, .TRUE., .TRUE., LSUCCS)`. Returns `LSUCCS=.FALSE.`
without populating GENSTK.

**Question.** Is this a known limitation of `SPDCEV` (e.g., specific
isotopes excluded)? Is there a list of decay-data isotopes for which
`SPDCEV` is *not* a valid sampler? Is there an alternative sampler we
should call for the affected entries?

**Author response (2026-05).** SPDCEV is the *correlated* decay
sampler — it requires unambiguous knowledge of the daughter nuclear
level so the decay kinematics close exactly. Most α decays in heavy
nuclei (and many β chains, e.g. Ac-228) lack that level information,
so `LSUCCS=.FALSE.` is by design. FLUKA has a separate *uncorrelated*
routine that walks the chain to a stable isotope regardless of level
information; it's appropriate for inventory/dose-rate work but not for
in-flight kinematics, where uncorrelated decays produce buggy boost
behaviour. The author offered a simplified version of that routine
that would close the isotope chain (e.g. Th-232 → Pb-208) at the cost
of non-physical β/γ/α kinematics.

**Decision (chromo, 2026-05).** Do not invoke any uncorrelated
fallback yet. `chromo_dcy_sample` returns the `LSUCCS=.FALSE.` flag
unchanged; `FlukaDecay._sample_one` and `DecayChainHandler.expand`
emit a one-time per-(A,Z,m) `UserWarning` and yield an empty event /
terminate the chain branch. To revisit after further discussion with
the author about scoping (correlated-only vs. opt-in uncorrelated mode).

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

**Author response (2026-05).** The routine will be included in the
next FLUKA respin. Until then the author provided the source directly
(his email spelled it `QRRDCY` but the attached `dcytst.f` snapshot
spells it `QRDDCY` — a transcription typo, not an upstream rename).

**Status (chromo, 2026-05-04).** The embedded copy in `chromo_fluka.f`
already matches the author-supplied snapshot byte-for-byte (vintage
`25-Apr-26 by Alfredo Ferrari`); cross-checked against
`dcytst_xAnatoli.f`.  Embed stays with a `TODO(FLUKA-bump): drop in
favour of upstream library symbol` comment until the next FLUKA
respin.

---

## 3. EMD-only event generation (IFLXYZ=100) aborts

**Problem.** `SGMXYZ(…, 100)` returns valid EMD cross sections.
`EVTXYZ(…, 100)` aborts.

**Reproduce.** `EVTXYZ(kproj, mat, ekin, 0, 0, 0, 1, 100)`.

**Question.** Is EMD-only event generation supported, or must it be
combined with inelastic (IFLXYZ=101)?

---

## 6. `SPDCEV` for beta decays crashes after `EVTXYZ`

**Problem.** Calling `SPDCEV(A, Z, 0, ...)` for a beta-emitter
*after* one or more `EVTXYZ` hadronic events have run aborts with a
Fortran array-bound error in `lcendp.f` line 74:

```
At line 74 of file lcendp.f
Fortran runtime error: Index '0' of dimension 2 of array 'edpsco'
below lower bound of 1
```

Standalone (no prior `EVTXYZ`) `SPDCEV` works for the same isotopes.
Pure α-decay (e.g. U-238) sampling *after* `EVTXYZ` works. Only β−/β+
samplers go through `lcendp.f` and trip the bound. The same condition
appears for tritium, Cs-137, Co-60, and likely every β-emitter.

**Reproduce.**
```fortran
CALL STPXYZ(...)              ! standard FLUKA bootstrap
CALL chromo_dcy_init           ! decay-table init (NCDTRD/RDFLUO/...)
CALL EVTXYZ(...)               ! one hadronic event, e.g. p + O16 @ 100 GeV
CALL SPDCEV(137, 55, 0, 0d0, 0d0, 0d0, 0d0, 1d0, 1d0, -1d9,
            .TRUE., .TRUE., LSUCCS)  ! crashes
```

**Suspected cause.** `EVTXYZ` modifies a common-block array (`EDPSCO`?
shared between hadronic event scoring and β-endpoint tables) in a way
that leaves the index `LFNDP`/`KFNDP` (or similar) at an invalid
value when `lcendp.f:74` next reads it.

**Impact for chromo.** `Fluka(..., post_event=DecayChainHandler(...))`
crashes the moment a β-decayable residual (tritium, Be-7, C-14, P-32,
Cs-137, …) is encountered. Currently we work around this by
documenting that the pattern is unsafe; a clean fix needs the FLUKA
side to either reset the relevant common before SPDCEV or expose a
"reset" entry-point we can call between `EVTXYZ` and `SPDCEV`.

**Question.** Is there an init/reset routine we should call between
`EVTXYZ` and `SPDCEV`? Is interleaving the two within one process a
supported usage pattern?

**Author response (2026-05).** `IPRODC = 2` (in common `(TRACKR)`)
must hold *before each* `SPDCEV` call. It's not a one-time init flag;
`EVTXYZ` overwrites it during normal hadronic tracking.

**Fix (chromo).** `chromo_dcy_sample` (in `chromo_fluka.f`) now
includes `(TRACKR)` and resets `IPRODC = 2`, `MRTRCK = 1`,
`MMTRCK = MEDFLK(MRTRCK, 1)` immediately before `CALL SPDCEV`. The
two extra TRACKR fields cover a secondary abort observed in
`geoden.f:100` (`MEDFLK` indexed at 0 because `MRTRCK` was 0 after
EVTXYZ). Regression tests:
`tests/test_fluka_decay.py::test_dcy_sample_cs137_after_evtxyz`,
`...::test_dcy_sample_co60_after_evtxyz`. With this fix, interleaving
`EVTXYZ` and `SPDCEV` (the `Fluka(post_event=DecayChainHandler(...))`
pattern) is unblocked for correlated isotopes.

---

## Build environment

- FLUKA 2025.1, DPMVERS=3.19.3.2
- `ldpmqmd`-equivalent link (libflukahp + libdpmmvax + librqmdmvax +
  libdpmjet3.19.3.2 + librqmd)
- gfortran (GCC toolset, Rocky Linux 9)
- FLRNDM for all random numbers
