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

## 7. DPMJET-3 event generation silently aborts above the Peanut cut

**Problem.**  The moment `EVTXYZ`'s internal dispatcher hands off from
Peanut to DPMJET-3 (event lab momentum above
`_transition_peanut_dpmjet`, default 20 TeV), DPMJET-3 reaches a bare
`STOP` at `DT_KKINC.f:73` — a silent process exit (code 0, no
backtrace, no message on stderr; gfortran's `STOP` without a string
prints nothing).  This reproduces in **a minimal F77 driver linked
straight against `libflukahp + libdpmjet3.19.3.2 + librqmd + libdpmmvax
+ librqmdmvax`, with no chromo wrapping at all** — so it's a FLUKA
library-level issue, not a chromo issue.

**Reproduce.**  Self-contained F77 driver (≈120 lines, no chromo
includes):

```fortran
      PROGRAM TEST_DPMJET_DRV
      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IOUNIT)'
      INCLUDE '(BEAMCM)'
      INCLUDE '(FLKMAT)'
      INCLUDE '(GENSTK)'
      INCLUDE '(NUCDAT)'
      INCLUDE '(PAPROP)'
      INCLUDE '(PAREVT)'
      INCLUDE '(PART2)'

      INTEGER NELMFL(4), IZELFL(4), MTFLKA(4)
      DOUBLE PRECISION WFELFL(4)
      DOUBLE PRECISION CUMSGI(0:NELEMX), CUMSGE(0:NELEMX),
     &                 CUMSGM(0:NELEMX)
      LOGICAL LPRINT, LSEEDI
      CHARACTER CRVRCK*8

      OPEN(UNIT=LUNRDB, FILE='ranpeaxyz', STATUS='OLD')
      IJKL = -1
      CALL RNREAD(LUNRDB, IJKL, LSEEDI)
      CLOSE(UNIT=LUNRDB)
      ICR = NINT(AMPRMU * 1.D+12 - 1.0072D+12)
      WRITE(CRVRCK,'(I8)') ICR

      NMATFL    = 1
      NELMFL(1) = 1
      IZELFL(1) = 1
      WFELFL(1) = ONEONE
      PPTMAX    = 4.8D+04
      EF2DP3    = -ONEONE
      DF2DP3    = -ONEONE
      IFLXYZ    = 1
      LPRINT    = .TRUE.

      CALL STPXYZ(NMATFL, NELMFL, IZELFL, WFELFL, 4,
     &            PPTMAX, EF2DP3, DF2DP3, IFLXYZ, LPRINT,
     &            MTFLKA, CRVRCK)

      ELAB     = 4.8D+04            ! switch to 1.D+02 for low-E
      KP       = 1
      IPFLK_KP = KPTOIP(KP)
      EKIN     = ELAB - AMPRTN
      PPROJ    = SQRT(ELAB*ELAB - AMPRTN*AMPRTN)

      DO IEV = 1, 3
         CALL EVTXYZ(IPFLK_KP, MTFLKA(1), EKIN, PPROJ,
     &               ZERZER, ZERZER, ONEONE, IFLXYZ,
     &               CUMSGI, CUMSGE, CUMSGM)
         WRITE(6,*) ' event ', IEV, ' NP=', NP
      END DO
      STOP
      END
```

Build:

```bash
gfortran -ffixed-form -fno-automatic -finit-local-zero \
  -frecord-marker=4 -fd-lines-as-comments \
  -I$FLUPRO/flukapro -I$FLUPRO/aamodmvax \
  -c dpmjet_drv.f
gfortran -o dpmjet_drv dpmjet_drv.o fluka_all.a
```

Same driver, only `ELAB` differs (with `PPTMAX` fixed at 48 TeV
throughout):

| `ELAB` | Path inside FLUKA | Result |
|---|---|---|
| 100 GeV | Peanut | ✅ 3 events generated, `NP = 4, 17, 12` |
| 48 000 GeV (= 48 TeV) | DPMJET-3 | ❌ silent crash inside event 1 (only an `IEEE_UNDERFLOW_FLAG` note) |

**Backtrace at the crash** (lldb breakpoint on `_gfortran_stop_string`):

```
exit
_gfortran_stop_string + 20      ← STOP with NULL string, len=0 (silent)
dt_kkinc_     at DT_KKINC.f:73
dpmrun_       at DPMRUN.f:69
eventd_       at eventd.f:6
eventp_       at eventp.f:1354
evtxyz_       at evtxyz.f:528
```

**Diagnostic message** (with `LPri = 10` set in `(DTFLKA)` before the
event loop, so DPMJET writes to its own log unit `LOUDPM = 19`):

```
fort.19:
   Requested energy (0.480E+05 GeV) exceeds initialization energy (0.000E+00 GeV) !
```

The "initialization energy" reported as `0.000E+00 GeV` is the
suspicious value — DPMJET clearly thinks it was never initialised
for this energy.

**What's set after `STPXYZ` returns** (read from common `(DTFLKA)`):

| Variable | Value | Comment |
|---|---|---|
| `EHFLLO` | 15 GeV | low energy bound |
| `EHFLHI` | **30 GeV** | high energy bound — well below event plab 48 TeV |
| `AMXPFR` | 0.27 GeV | maximum Fermi momentum (looks correct) |
| `LPri` | 0 | DPMJET verbose flag (off by default) |

`fort.11` (FLUKA's main log) reports `Max. initialization energy for
DPMJET: 86331.4 GeV` — i.e. ≈ 1.8 × `PPTMAX`.  So FLUKA's outer init
*does* compute a sensible high-energy bound, but the value that
`DT_KKINC` checks against (and that `LPri≥4` prints as "0.000E+00 GeV"
in `fort.19`) is *not* `EHFLHI` and *not* the 86 TeV grid bound — it
is a DPMJET-internal variable set only by sub-inits that aren't
running.  Overriding `EHFLLO`/`EHFLHI`/`AMXPFR` to large values from
outside (via a small `INCLUDE '(DTFLKA)'` setter) does not unblock
the crash.

**`dt_init_` is called** from `dpmini_` line 117 with arguments:

```
NCASES = -1                ! variable-energy mode flag
EPN    = 86331.4 GeV       ! max projectile energy per nucleon (correct)
NPMASS = 12, NPCHAR = 6    ! C-12 placeholder (we asked for p projectile)
NTMASS = 0,  NTCHAR = 0    ! placeholder target
```

After `dt_init_` returns, **the only DPMJET sub-init invoked** (lldb
breakpoint regex `^dt_.*ini.*`) is `dt_ltini_` (at `DT_LTINI.f:33`).
Notably **none** of `dt_dtuini_`, `dt_phoini_`, `dt_glbini_`,
`dt_evtini_`, `dt_xhoini_` are called.  These are the routines that
populate the DPMJET-internal energy bounds, so the per-event check at
`DT_KKINC.f:73` always sees zero and rejects.

**Workaround that works in chromo today.**  `Fluka(...,
transition_peanut_dpmjet=100*TeV)` keeps Peanut on the path; DPMJET-3
is never dispatched; events run cleanly up to plab ≈ 192 TeV (Peanut's
own ceiling).  This is the documented workaround, but it caps the
useful range to Peanut's accuracy region.

**Things that did NOT unblock the crash:**

1. Bumping `PPTMAX` to 10 × `plab` (DPMJET Glauber grid → 863 TeV per
   `fort.11`) — it's a per-event init-energy state, not a
   grid-coverage problem.
2. Overriding `EHFLLO`/`EHFLHI`/`AMXPFR` in `(DTFLKA)` to large values
   from outside before the event loop.  The "init energy" `DT_KKINC`
   reads is a different DPMJET-private value (set only by sub-inits
   that aren't running).
3. Calling `SGMXYZ(KPROJ, MMAT, EKIN, PPROJ, IFLXYZ)` before
   `EVTXYZ` for each event — this is what the original
   `prinmvax/peaxyz.f` does (line 430).  `SGMXYZ` returns a sensible
   inelastic cross section (44 mb for p+p at 48 TeV plab — DPMJET's
   Glauber path is fine), but `EVTXYZ` immediately afterwards still
   trips the same bare `STOP` in `DT_KKINC.f:73`.

**Tested: linking `afedynitch/fluka_chromo`'s `prinmvax/stpxyz.f` in
place of `libflukahp.a`'s `STPXYZ` does NOT fix it.** That `STPXYZ`
sets `LFLUKA=.F.`, `LCRSKA=.T.`, `LHEPCM=.T.`, `LIONTR=.T.`,
coalescence flags, `IJBEAM=1`, `PBEAM=PPTMAX`, and calls `EVXINI`
(= current `INEVTI`) before the material loop — the long "Corsika-mode"
preamble we suspected was the missing piece.

Test setup: same F77 driver as above, with `prinmvax/stpxyz.f`
compiled and linked before `fluka_all.a` so its `_stpxyz_` wins
(confirmed by a runtime trace marker we added in the override). One
edit only: rename `EVXINI → INEVTI` to match the FLUKA 2025.1 symbol.
Driver run at `ELAB = 1.333E+09 GeV` (≈ 50 TeV CMS), `PPTMAX = 1.4E+09
GeV`, hydrogen target.

Result, identical to the `libflukahp.a` `STPXYZ` run:

```
fort.11:  Max. initialization energy for DPMJET   2.52E+09 GeV    (= 1.8 × PPTMAX, OK)
fort.19:  Requested energy (0.133E+10GeV) exceeds initialization energy (0.000E+00GeV) !
```

Silent exit at `DT_KKINC.f:73`, same diagnostic, same `0.000E+00 GeV`
internal value. `_stpxyz_` from the override is the one that ran;
the issue isn't whatever `STPXYZ` does. Same driver at `ELAB = 100
GeV` (Peanut path) generates the expected 3 events.

So the `STPXYZ`-side flags (`LCRSKA`, `IJBEAM`, `PBEAM`, …) don't
affect the per-event `0.000E+00 GeV` reject. The variable
`DT_KKINC` checks lives in `(DTVARE)` (binary scan: only
`DT_DEFAUL`, `DT_INIT`, `DT_SHMAKI` touch it), and is left at zero by
the `dt_init` path taken under `NCASES=-1, NPMASS=12, NPCHAR=6`.

**Question, narrowed.**  Under `NCASES=-1` (variable-energy mode) with
the C-12 placeholder that `DPMINI` passes to `dt_init`, what populates
the `(DTVARE)` energy bound that `DT_KKINC` then checks per event?
Is there a routine the library-mode caller is expected to invoke
between `INEVTI` and the first `EVTXYZ` (or per event) to bring this
bound up to `EPN`?

---

## Build environment

- FLUKA 2025.1, DPMVERS=3.19.3.2
- `ldpmqmd`-equivalent link (libflukahp + libdpmmvax + librqmdmvax +
  libdpmjet3.19.3.2 + librqmd)
- gfortran (GCC toolset, Rocky Linux 9)
- FLRNDM for all random numbers
