# FLUKA Radioactive-Decay Interface — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a `FlukaDecay` class to chromo that exposes FLUKA 2025.1's nuclear-decay tables and event-sampling primitives (catalog query, isotope inspection, inclusive decay sampling, recursive decay chains), plus an opt-in chain post-processor for `Fluka` hadronic events.

**Architecture:** Six new thin Fortran wrappers in `src/fortran/fluka/chromo_fluka.f` expose `ISMRCH`, table-walk, `BRDECY`/`IDCNUC`/`QRDDCY`, line-list pointers, and `SPDCEV`+`FLLHEP`. A new Python module `src/chromo/models/fluka_decay.py` builds `FlukaDecay`, `FlukaIsotope`, `DecayChannel`, `DecayLine`, `DecayChainHandler`, and `STABLE_DEFAULT` on top. `Fluka` gains an opt-in `post_event` callable for chain post-processing. Single-instantiation is enforced via a module-level guard shared with `Fluka`.

**Tech Stack:** Fortran (gfortran via meson + f2py), Python ≥ 3.9, numpy, particle, pytest, chromo's existing `MCRun`/`EventData` machinery.

**Spec:** [`docs/superpowers/specs/2026-04-30-fluka-decay-interface-design.md`](../specs/2026-04-30-fluka-decay-interface-design.md)

---

## File Structure

**Created:**
- `src/chromo/models/fluka_decay.py` — `FlukaDecay`, `FlukaIsotope`, `DecayChannel`, `DecayLine`, `DecayChainHandler`, `STABLE_DEFAULT`, the init guard, and lazy-fetch helpers.
- `tests/test_fluka_decay.py` — full test matrix (all run via `run_in_separate_process`).

**Modified:**
- `src/fortran/fluka/chromo_fluka.f` — append 6 wrapper subroutines (`CHROMO_DCY_INIT`, `CHROMO_DCY_LOOKUP`, `CHROMO_DCY_CATALOG`, `CHROMO_DCY_CHANNELS`, `CHROMO_DCY_LINES`, `CHROMO_DCY_SAMPLE`).
- `meson.build` — extend `fluka_syms` list (around line 432) with the 6 new symbols.
- `src/chromo/models/__init__.py` — re-export `FlukaDecay`.
- `src/chromo/models/fluka.py` — add `post_event` kwarg to `Fluka.__init__` and apply it inside `_generate()` before yielding.

**Build command after Fortran/meson edits** (run from repo root):
```bash
pip install --no-build-isolation -v -e .[test] 2>&1 | tail -20
```

---

## Task 1: Add `CHROMO_DCY_INIT` Fortran wrapper

**Files:**
- Modify: `src/fortran/fluka/chromo_fluka.f` (append at EOF)
- Modify: `meson.build` (extend `fluka_syms` list around line 432)
- Test: `tests/test_fluka_decay.py` (create)

- [ ] **Step 1: Write failing import test**

Create `tests/test_fluka_decay.py`:

```python
"""Tests for chromo's FLUKA radioactive-decay interface."""
from tests.util import run_in_separate_process


def _import_smoke():
    from chromo.models import _fluka

    assert hasattr(_fluka, "chromo_dcy_init")
    return True


def test_chromo_dcy_init_symbol_exists():
    assert run_in_separate_process(_import_smoke) is True


def _call_init():
    from chromo.models import _fluka

    _fluka.chromo_dcy_init()
    return True


def test_chromo_dcy_init_runs():
    assert run_in_separate_process(_call_init) is True
```

- [ ] **Step 2: Run test, verify failure**

```bash
python -m pytest tests/test_fluka_decay.py::test_chromo_dcy_init_symbol_exists -v
```
Expected: FAIL with `AttributeError: module ... has no attribute 'chromo_dcy_init'`.

- [ ] **Step 3: Append wrapper to chromo_fluka.f**

Append to `src/fortran/fluka/chromo_fluka.f` (after the existing `fluka_hepevt_summary` subroutine):

```fortran
!======================================================================!
!  Radioactive-decay tables: lightweight init (no STPXYZ).             !
!  Mirrors dcytst.f's main-program init sequence.  Idempotent: safe    !
!  to call after STPXYZ (which already loads the same tables).         !
!======================================================================!
      SUBROUTINE chromo_dcy_init()
      INCLUDE '(DBLPRW)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IOUNIT)'

      EXTERNAL BDNOPT, BDEVAP, BDPREE, BDINPT, BDTRNS, BDPART, BDPRDC

      INCLUDE '(EMFCMP)'
      INCLUDE '(EMFMAT)'
      INCLUDE '(EMFSWT)'
      INCLUDE '(EMFTHR)'
      INCLUDE '(EVAPRD)'
      INCLUDE '(FLUOXR)'
      INCLUDE '(FLKMAT)'
      INCLUDE '(FRBKCM)'
      INCLUDE '(NUCDAT)'
      INCLUDE '(PAPROP)'
      INCLUDE '(PAREVT)'
      INCLUDE '(TRACKR)'

      INTEGER MMAT, MREG
      LOGICAL, SAVE :: LDCYINI = .FALSE.

      IF (LDCYINI) RETURN

      CALL CMSPPR
      CALL ZEROIN
      LEVPRT = .TRUE.
      LDEEXG = .TRUE.
      LHEAVY = .TRUE.
      LGDHPR = .TRUE.
      LFRMBK = .FALSE.
      CALL NCDTRD
      CALL KPIXSR
      CALL INCINI
      AMUMEV = GEVMEV * AMUAMU

      LFLUKA = .FALSE.
      LEMFON = .TRUE.
      MMAT   = 3
      NMAT   = MMAT
      MREG   = 1
      IPRODC = 2
      MMTRCK = MMAT
      MRTRCK = MREG
      NMDEMF = 1
      METOFL(NMDEMF) = MMAT
      MFLTOE(MMAT)   = NMDEMF
      LXFLUO(NMDEMF) = .TRUE.
      MEDFLK(MREG,1) = MMAT
      MEDFLK(MREG,2) = MMAT
      EPHMIN(NMDEMF) = 1.D-07 * GV2EMF
      EEPMIN(NMDEMF) = (AMELCT + 1.D-07) * GV2EMF
      CALL RDFLUO

      LDCYINI = .TRUE.
      RETURN
      END SUBROUTINE chromo_dcy_init
```

- [ ] **Step 4: Add symbol to meson.build**

Modify `meson.build` around line 432 (the `fluka_syms = [...]` block):

```python
  fluka_syms = [
    'chromo_evtxyz', 'chromo_stpxyz', 'chromo_sgmxyz', 'chromo_fllhep',
    'pdg_to_proj_code', 'pdg_to_evt_code',
    'fluka_elem_properties', 'fluka_hepevt_summary',
    'random_direction',
    'icode_from_pdg', 'icode_from_pdg_arr', 'charge_from_pdg_arr',
    'pdg_from_icode',
    'init_rng_state', 'load_rng_state', 'save_rng_state', 'fluka_rand',
    'fluka_particle_scheme',
    'chromo_dcy_init',
  ]
```

- [ ] **Step 5: Rebuild**

```bash
pip install --no-build-isolation -v -e .[test] 2>&1 | tail -5
```
Expected: `Successfully installed chromo-...`.

- [ ] **Step 6: Run tests, verify pass**

```bash
python -m pytest tests/test_fluka_decay.py -v
```
Expected: both tests PASS.

- [ ] **Step 7: Commit**

```bash
git add src/fortran/fluka/chromo_fluka.f meson.build tests/test_fluka_decay.py
git commit -m "feat(fluka): add chromo_dcy_init wrapper"
```

---

## Task 2: Add `CHROMO_DCY_LOOKUP` Fortran wrapper

**Files:**
- Modify: `src/fortran/fluka/chromo_fluka.f`
- Modify: `meson.build`
- Modify: `tests/test_fluka_decay.py`

- [ ] **Step 1: Write failing test**

Append to `tests/test_fluka_decay.py`:

```python
def _lookup_u238():
    from chromo.models import _fluka

    _fluka.chromo_dcy_init()
    found, t12, exm, jsp, jpt = _fluka.chromo_dcy_lookup(238, 92, 0)
    return int(found), float(t12), float(exm)


def test_dcy_lookup_u238():
    found, t12, exm = run_in_separate_process(_lookup_u238)
    assert found == 1
    assert abs(t12 - 1.41e17) / 1.41e17 < 0.01      # ~4.5 Gyr
    assert abs(exm - 47.31) < 0.05                   # MeV


def _lookup_invalid():
    from chromo.models import _fluka

    _fluka.chromo_dcy_init()
    found, _t12, _exm, _jsp, _jpt = _fluka.chromo_dcy_lookup(2, 50, 0)
    return int(found)


def test_dcy_lookup_invalid_returns_zero():
    assert run_in_separate_process(_lookup_invalid) == 0
```

- [ ] **Step 2: Run, verify failure**

```bash
python -m pytest tests/test_fluka_decay.py::test_dcy_lookup_u238 -v
```
Expected: FAIL (`chromo_dcy_lookup` doesn't exist).

- [ ] **Step 3: Append wrapper**

Append to `src/fortran/fluka/chromo_fluka.f`:

```fortran
!======================================================================!
!  ISMRCH wrapper: returns IFOUND=1 if (A,Z,m) has decay data,         !
!  IFOUND=0 otherwise.  Caller is responsible for staying within the   !
!  safe (A,Z) probe band (see chromo_dcy_catalog).                     !
!======================================================================!
      SUBROUTINE chromo_dcy_lookup(ia, iz, im, ifound,
     &                             t12, exm, jsp, jpt)
      INCLUDE '(DBLPRW)'
      INCLUDE '(DIMPAR)'

      INTEGER ia, iz, im, ifound, jsp, jpt
      DOUBLE PRECISION t12, exm
Cf2py intent(in)  ia, iz, im
Cf2py intent(out) ifound, t12, exm, jsp, jpt

      INTEGER kisitp

      kisitp = 0
      ifound = 0
      t12 = 0.0D0
      exm = 0.0D0
      jsp = 0
      jpt = 0

!  Guard against ISMRCH OOB on pathological inputs:
      IF (ia .LT. 1) RETURN
      IF (iz .LT. 0) RETURN
      IF (ia .EQ. 1 .AND. iz .GT. 1) RETURN
      IF (ia .GT. 1 .AND. iz .LT. 1) RETURN
      IF (im .LT. 0 .OR. im .GT. 4) RETURN

      CALL ISMRCH(ia, iz, im, kisitp, t12, exm, jsp, jpt)

      IF (kisitp .GT. 0) ifound = 1

      RETURN
      END SUBROUTINE chromo_dcy_lookup
```

- [ ] **Step 4: Add to meson.build**

Add `'chromo_dcy_lookup'` to the `fluka_syms` list in `meson.build`.

- [ ] **Step 5: Rebuild**

```bash
pip install --no-build-isolation -v -e .[test] 2>&1 | tail -5
```

- [ ] **Step 6: Run tests, verify pass**

```bash
python -m pytest tests/test_fluka_decay.py -v
```
Expected: 4 PASS.

- [ ] **Step 7: Commit**

```bash
git add src/fortran/fluka/chromo_fluka.f meson.build tests/test_fluka_decay.py
git commit -m "feat(fluka): add chromo_dcy_lookup wrapper (ISMRCH)"
```

---

## Task 3: Add `CHROMO_DCY_CATALOG` Fortran wrapper

**Files:**
- Modify: `src/fortran/fluka/chromo_fluka.f`
- Modify: `meson.build`
- Modify: `tests/test_fluka_decay.py`

- [ ] **Step 1: Write failing test**

Append to `tests/test_fluka_decay.py`:

```python
def _catalog_size():
    import numpy as np
    from chromo.models import _fluka

    _fluka.chromo_dcy_init()

    max_n = 5500
    a = np.zeros(max_n, dtype=np.int32)
    z = np.zeros(max_n, dtype=np.int32)
    m = np.zeros(max_n, dtype=np.int32)
    t12 = np.zeros(max_n, dtype=np.float64)
    exm = np.zeros(max_n, dtype=np.float64)
    jsp = np.zeros(max_n, dtype=np.int32)
    jpt = np.zeros(max_n, dtype=np.int32)

    n = _fluka.chromo_dcy_catalog(max_n, a, z, m, t12, exm, jsp, jpt)
    return int(n), int(a[:n].max()), int(z[:n].max())


def test_dcy_catalog_size():
    n, a_max, z_max = run_in_separate_process(_catalog_size)
    assert 4000 <= n <= 5500
    assert 200 < a_max <= 295
    assert 80 < z_max <= 110
```

- [ ] **Step 2: Run, verify failure**

```bash
python -m pytest tests/test_fluka_decay.py::test_dcy_catalog_size -v
```
Expected: FAIL.

- [ ] **Step 3: Append wrapper**

Append to `src/fortran/fluka/chromo_fluka.f`:

```fortran
!======================================================================!
!  Bulk walk of the FLUKA decay-data catalogue.  Conservative valley-  !
!  of-stability Z-band per A; ground-state probe gates isomer probes.  !
!  Returns N total entries written to caller-allocated arrays.         !
!======================================================================!
      SUBROUTINE chromo_dcy_catalog(max_n, n_out,
     &                              a_out, z_out, m_out,
     &                              t12_out, exm_out,
     &                              jsp_out, jpt_out)
      INCLUDE '(DBLPRW)'
      INCLUDE '(DIMPAR)'

      INTEGER max_n, n_out
      INTEGER a_out(max_n), z_out(max_n), m_out(max_n)
      INTEGER jsp_out(max_n), jpt_out(max_n)
      DOUBLE PRECISION t12_out(max_n), exm_out(max_n)
Cf2py intent(in)  max_n
Cf2py intent(in,out) a_out, z_out, m_out, t12_out, exm_out, jsp_out, jpt_out
Cf2py intent(out) n_out

      INTEGER ia, iz, im, izmin, izmax, kisitp, jsp, jpt
      DOUBLE PRECISION t12, exm

      n_out = 0
      DO ia = 1, 295
         IF (ia .LE. 4) THEN
            izmin = 1
            izmax = MIN(ia, 110)
         ELSE IF (ia .LE. 20) THEN
            izmin = MAX(1, INT(0.30D0*ia) - 4)
            izmax = MIN(ia, 110)
         ELSE
            izmin = MAX(1, INT(0.30D0*ia) - 5)
            izmax = MIN(110, INT(0.55D0*ia) + 8)
         END IF
         DO iz = izmin, izmax
!  Probe ground state first; only try isomers if g.s. exists.
            kisitp = 0
            CALL ISMRCH(ia, iz, 0, kisitp, t12, exm, jsp, jpt)
            IF (kisitp .GT. 0) THEN
               IF (n_out .GE. max_n) RETURN
               n_out = n_out + 1
               a_out(n_out)   = ia
               z_out(n_out)   = iz
               m_out(n_out)   = 0
               t12_out(n_out) = t12
               exm_out(n_out) = exm
               jsp_out(n_out) = jsp
               jpt_out(n_out) = jpt
               IF (ia .GE. 5) THEN
                  DO im = 1, 4
                     kisitp = 0
                     CALL ISMRCH(ia, iz, im, kisitp,
     &                           t12, exm, jsp, jpt)
                     IF (kisitp .GT. 0) THEN
                        IF (n_out .GE. max_n) RETURN
                        n_out = n_out + 1
                        a_out(n_out)   = ia
                        z_out(n_out)   = iz
                        m_out(n_out)   = im
                        t12_out(n_out) = t12
                        exm_out(n_out) = exm
                        jsp_out(n_out) = jsp
                        jpt_out(n_out) = jpt
                     END IF
                  END DO
               END IF
            END IF
         END DO
      END DO

      RETURN
      END SUBROUTINE chromo_dcy_catalog
```

- [ ] **Step 4: Add to meson.build**

Add `'chromo_dcy_catalog'` to `fluka_syms`.

- [ ] **Step 5: Rebuild**

```bash
pip install --no-build-isolation -v -e .[test] 2>&1 | tail -5
```

- [ ] **Step 6: Run tests, verify pass**

```bash
python -m pytest tests/test_fluka_decay.py -v
```
Expected: 5 PASS.

- [ ] **Step 7: Commit**

```bash
git add src/fortran/fluka/chromo_fluka.f meson.build tests/test_fluka_decay.py
git commit -m "feat(fluka): add chromo_dcy_catalog wrapper"
```

---

## Task 4: Add `CHROMO_DCY_CHANNELS` Fortran wrapper

**Files:**
- Modify: `src/fortran/fluka/chromo_fluka.f`
- Modify: `meson.build`
- Modify: `tests/test_fluka_decay.py`

- [ ] **Step 1: Write failing test**

Append to `tests/test_fluka_decay.py`:

```python
def _channels_cs137():
    import numpy as np
    from chromo.models import _fluka

    _fluka.chromo_dcy_init()
    max_ch = 8
    kind = np.zeros(max_ch, dtype=np.int32)
    br = np.zeros(max_ch, dtype=np.float64)
    da = np.zeros(max_ch, dtype=np.int32)
    dz = np.zeros(max_ch, dtype=np.int32)
    dm = np.zeros(max_ch, dtype=np.int32)
    qv = np.zeros(max_ch, dtype=np.float64)

    n = _fluka.chromo_dcy_channels(137, 55, 0,
                                   max_ch, kind, br, da, dz, dm, qv)
    return int(n), kind[:n].tolist(), br[:n].tolist()


def test_dcy_channels_cs137():
    n, kinds, brs = run_in_separate_process(_channels_cs137)
    assert n >= 1
    # All recorded channels should be Beta- (kind=2) for Cs-137
    assert all(k == 2 for k in kinds), kinds
    assert abs(sum(brs) - 1.0) < 0.01
```

- [ ] **Step 2: Run, verify failure**

```bash
python -m pytest tests/test_fluka_decay.py::test_dcy_channels_cs137 -v
```
Expected: FAIL.

- [ ] **Step 3: Append wrapper**

Append to `src/fortran/fluka/chromo_fluka.f`:

```fortran
!======================================================================!
!  Per-isotope decay-channel data: branching ratios, daughter (A,Z,m), !
!  Q value (via QRDDCY).  KIND code: 1=alpha, 2=B-, 3=B+, 4=EC, 5=IT,  !
!  6=SF, 7=B-N, 8=B+P, 9=B-2N, 10=B-3N, 11=B-NA, 12=other.             !
!======================================================================!
      SUBROUTINE chromo_dcy_channels(ia, iz, im, max_ch, n_ch,
     &                               kind_ch, br_ch,
     &                               da_ch, dz_ch, dm_ch, q_ch)
      INCLUDE '(DBLPRW)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IDPPRM)'
      INCLUDE '(ISOTOP)'

      INTEGER ia, iz, im, max_ch, n_ch
      INTEGER kind_ch(max_ch), da_ch(max_ch),
     &        dz_ch(max_ch), dm_ch(max_ch)
      DOUBLE PRECISION br_ch(max_ch), q_ch(max_ch)
Cf2py intent(in)     ia, iz, im, max_ch
Cf2py intent(in,out) kind_ch, br_ch, da_ch, dz_ch, dm_ch, q_ch
Cf2py intent(out)    n_ch

      INTEGER kisitp, jsporu, jptoru, iia, jdcy, kdcy
      DOUBLE PRECISION t12, exm, brdctt
      DOUBLE PRECISION QRDDCY
      EXTERNAL QRDDCY

      n_ch = 0
      kisitp = 0
      CALL ISMRCH(ia, iz, im, kisitp, t12, exm, jsporu, jptoru)
      IF (kisitp .LE. 0) RETURN

      iia = KQMDIN(ia, KMSFVR)
      brdctt = 0.0D0

      DO jdcy = 1, MXDCIS
         IF (im .LE. 0) THEN
            kdcy = IDCNUC(jdcy, iia, kisitp)
            IF (jdcy .LT. MXDCIS) THEN
               br_ch(MIN(n_ch+1, max_ch)) = BRDECY(jdcy, iia, kisitp)
            ELSE
               br_ch(MIN(n_ch+1, max_ch)) = 1.0D0 - brdctt
            END IF
         ELSE
            kdcy = IDCISM(jdcy, kisitp)
            IF (jdcy .LT. MXDCIS) THEN
               br_ch(MIN(n_ch+1, max_ch)) = BRDISM(jdcy, kisitp)
            ELSE
               br_ch(MIN(n_ch+1, max_ch)) = 1.0D0 - brdctt
            END IF
         END IF

         IF (kdcy .LE. 0) CYCLE
         IF (n_ch .GE. max_ch) RETURN
         n_ch = n_ch + 1
         brdctt = brdctt + br_ch(n_ch)

!  Map FLUKA kdcy to a compact KIND code.  See (ISOTOP) for KDCY*
!  parameters: KDCYAL=2 (alpha), KDCYBM=22 (B-), KDCYBP=8 (B+),
!  KDCYEC=15 (EC), KDCYIT=1 (IT), KDCYSF=38 (SF), KDCBMN=23 (B-N),
!  KDCBPP=24 (B+P), KDCB2N=28 (B-2N), KDCB3N=29 (B-3N), KDCBNA=30 (B-NA).
         IF (kdcy .EQ. KDCYAL .OR. kdcy .EQ. KDCALM) THEN
            kind_ch(n_ch) = 1
         ELSE IF (kdcy .EQ. KDCYBM .OR. kdcy .EQ. KDCBMM) THEN
            kind_ch(n_ch) = 2
         ELSE IF (kdcy .EQ. KDCYBP .OR. kdcy .EQ. KDCBPM) THEN
            kind_ch(n_ch) = 3
         ELSE IF (kdcy .EQ. KDCYEC .OR. kdcy .EQ. KDCECM) THEN
            kind_ch(n_ch) = 4
         ELSE IF (kdcy .EQ. KDCYIT .OR. kdcy .EQ. KDCITM) THEN
            kind_ch(n_ch) = 5
         ELSE IF (kdcy .EQ. KDCYSF) THEN
            kind_ch(n_ch) = 6
         ELSE IF (kdcy .EQ. KDCBMN) THEN
            kind_ch(n_ch) = 7
         ELSE IF (kdcy .EQ. KDCBPP) THEN
            kind_ch(n_ch) = 8
         ELSE IF (kdcy .EQ. KDCB2N) THEN
            kind_ch(n_ch) = 9
         ELSE IF (kdcy .EQ. KDCB3N) THEN
            kind_ch(n_ch) = 10
         ELSE IF (kdcy .EQ. KDCBNA) THEN
            kind_ch(n_ch) = 11
         ELSE
            kind_ch(n_ch) = 12
         END IF

         IF (IDCYDA(kdcy) .GT. -100) THEN
            da_ch(n_ch) = ia + IDCYDA(kdcy)
            dz_ch(n_ch) = iz + IDCYDZ(kdcy)
            IF (kdcy .GT. NDCY1M) THEN
               dm_ch(n_ch) = 2
            ELSE IF (kdcy .GT. NDCYGS) THEN
               dm_ch(n_ch) = 1
            ELSE
               dm_ch(n_ch) = 0
            END IF
            q_ch(n_ch) = QRDDCY(ia, iz, im, kdcy, .FALSE.) * GEVMEV
         ELSE
            da_ch(n_ch) = -1
            dz_ch(n_ch) = -1
            dm_ch(n_ch) = -1
            q_ch(n_ch) = 0.0D0
         END IF
      END DO

      RETURN
      END SUBROUTINE chromo_dcy_channels
```

- [ ] **Step 4: Add to meson.build**

Add `'chromo_dcy_channels'` to `fluka_syms`.

- [ ] **Step 5: Rebuild**

```bash
pip install --no-build-isolation -v -e .[test] 2>&1 | tail -5
```

- [ ] **Step 6: Run tests**

```bash
python -m pytest tests/test_fluka_decay.py -v
```
Expected: 6 PASS.

- [ ] **Step 7: Commit**

```bash
git add src/fortran/fluka/chromo_fluka.f meson.build tests/test_fluka_decay.py
git commit -m "feat(fluka): add chromo_dcy_channels wrapper"
```

---

## Task 5: Add `CHROMO_DCY_LINES` Fortran wrapper

**Files:**
- Modify: `src/fortran/fluka/chromo_fluka.f`
- Modify: `meson.build`
- Modify: `tests/test_fluka_decay.py`

- [ ] **Step 1: Write failing test**

Append to `tests/test_fluka_decay.py`:

```python
def _gamma_lines_cs137():
    import numpy as np
    from chromo.models import _fluka

    _fluka.chromo_dcy_init()
    max_l = 64
    br = np.zeros(max_l, dtype=np.float64)
    e_mev = np.zeros(max_l, dtype=np.float64)
    nlev = np.zeros(max_l, dtype=np.int32)
    n = _fluka.chromo_dcy_lines(137, 55, 0, 1, max_l, br, e_mev, nlev)
    return int(n), e_mev[:n].tolist(), br[:n].tolist()


def test_dcy_lines_cs137_gamma():
    n, energies, brs = run_in_separate_process(_gamma_lines_cs137)
    assert n >= 1
    # Cs-137 famously has the 661.66 keV gamma
    assert any(abs(e - 0.66166) < 0.001 for e in energies), energies
```

- [ ] **Step 2: Run, verify failure**

```bash
python -m pytest tests/test_fluka_decay.py::test_dcy_lines_cs137_gamma -v
```
Expected: FAIL.

- [ ] **Step 3: Append wrapper**

Append to `src/fortran/fluka/chromo_fluka.f`:

```fortran
!======================================================================!
!  Per-isotope line lists.  KIND: 1=gamma (BR, E), 2=alpha (BR, E,     !
!  end-level), 3=CE/Auger (BR, E), 4=beta+/- (BR, <E>, end-point E,    !
!  end-level; positron flag encoded in sign of <E>).                   !
!  Energies returned in MeV.                                           !
!======================================================================!
      SUBROUTINE chromo_dcy_lines(ia, iz, im, kind, max_l, n_l,
     &                            br_l, e_l, nlev_l)
      INCLUDE '(DBLPRW)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(ISOTOP)'

      INTEGER ia, iz, im, kind, max_l, n_l
      INTEGER nlev_l(max_l)
      DOUBLE PRECISION br_l(max_l), e_l(max_l)
Cf2py intent(in)     ia, iz, im, kind, max_l
Cf2py intent(in,out) br_l, e_l, nlev_l
Cf2py intent(out)    n_l

      INTEGER kisitp, jsporu, jptoru, iia, kn_lines, kp_lines, ipt
      INTEGER ndata
      DOUBLE PRECISION t12, exm
      DOUBLE PRECISION SIGGTT
      EXTERNAL SIGGTT

      n_l = 0
      kisitp = 0
      CALL ISMRCH(ia, iz, im, kisitp, t12, exm, jsporu, jptoru)
      IF (kisitp .LE. 0) RETURN

      iia = KQMDIN(ia, KMSFVR)
      kn_lines = 0
      kp_lines = 0

      IF (im .LE. 0) THEN
         IF (kind .EQ. 1) THEN
            kn_lines = NGMLNS(iia, kisitp)
            kp_lines = KGMLNS(iia, kisitp)
            ndata = 2
         ELSE IF (kind .EQ. 2) THEN
            kn_lines = NALLNS(iia, kisitp)
            kp_lines = KALLNS(iia, kisitp)
            ndata = 3
         ELSE IF (kind .EQ. 3) THEN
            kn_lines = NCELNS(iia, kisitp)
            kp_lines = KCELNS(iia, kisitp)
            ndata = 2
         ELSE IF (kind .EQ. 4) THEN
            kn_lines = NBTSPC(iia, kisitp)
            kp_lines = KBTSPC(iia, kisitp)
            ndata = 4
         ELSE
            RETURN
         END IF
      ELSE
         IF (kind .EQ. 1) THEN
            kn_lines = NGMISM(kisitp)
            kp_lines = KGMISM(kisitp)
            ndata = 2
         ELSE IF (kind .EQ. 2) THEN
            kn_lines = NALISM(kisitp)
            kp_lines = KALISM(kisitp)
            ndata = 3
         ELSE IF (kind .EQ. 3) THEN
            kn_lines = NCEISM(kisitp)
            kp_lines = KCEISM(kisitp)
            ndata = 2
         ELSE IF (kind .EQ. 4) THEN
            kn_lines = NBTISM(kisitp)
            kp_lines = KBTISM(kisitp)
            ndata = 4
         ELSE
            RETURN
         END IF
      END IF

      DO ipt = 1, kn_lines
         IF (n_l .GE. max_l) RETURN
         n_l = n_l + 1
         br_l(n_l) = ABS(SIGGTT(kp_lines + ndata*(ipt-1) + 1))
         e_l(n_l)  = SIGGTT(kp_lines + ndata*(ipt-1) + 2) * GEVMEV
         IF (ndata .EQ. 3) THEN
            nlev_l(n_l) = NINT(SIGGTT(kp_lines + ndata*(ipt-1) + 3))
         ELSE IF (ndata .EQ. 4) THEN
!  Beta: index +2 holds <E> (signed), +3 end-point E, +4 end-level
            e_l(n_l) = SIGGTT(kp_lines + ndata*(ipt-1) + 3) * GEVMEV
            nlev_l(n_l) = NINT(SIGGTT(kp_lines + ndata*(ipt-1) + 4))
         ELSE
            nlev_l(n_l) = 0
         END IF
      END DO

      RETURN
      END SUBROUTINE chromo_dcy_lines
```

- [ ] **Step 4: Add to meson.build**

Add `'chromo_dcy_lines'` to `fluka_syms`.

- [ ] **Step 5: Rebuild**

```bash
pip install --no-build-isolation -v -e .[test] 2>&1 | tail -5
```

- [ ] **Step 6: Run tests**

```bash
python -m pytest tests/test_fluka_decay.py -v
```
Expected: 7 PASS.

- [ ] **Step 7: Commit**

```bash
git add src/fortran/fluka/chromo_fluka.f meson.build tests/test_fluka_decay.py
git commit -m "feat(fluka): add chromo_dcy_lines wrapper"
```

---

## Task 6: Add `CHROMO_DCY_SAMPLE` Fortran wrapper

**Files:**
- Modify: `src/fortran/fluka/chromo_fluka.f`
- Modify: `meson.build`
- Modify: `tests/test_fluka_decay.py`

- [ ] **Step 1: Write failing test**

Append to `tests/test_fluka_decay.py`:

```python
def _sample_cs137_one():
    from chromo.models import _fluka

    _fluka.chromo_dcy_init()
    success, kdcy, ilv = _fluka.chromo_dcy_sample(137, 55, 0)
    return bool(success), int(kdcy), int(ilv), int(_fluka.hepcmm.nhep)


def test_dcy_sample_cs137_succeeds():
    ok, kdcy, ilv, nhep = run_in_separate_process(_sample_cs137_one)
    assert ok
    assert kdcy == 2  # B-
    assert nhep >= 2   # at least e- + ν̄
```

- [ ] **Step 2: Run, verify failure**

```bash
python -m pytest tests/test_fluka_decay.py::test_dcy_sample_cs137_succeeds -v
```
Expected: FAIL.

- [ ] **Step 3: Append wrapper**

Append to `src/fortran/fluka/chromo_fluka.f`:

```fortran
!======================================================================!
!  Sample one correlated radioactive-decay event for (A,Z,m).  Clears  !
!  GENSTK / EMFSTK / FHEAVY counters, calls SPDCEV, then FLLHEP to put !
!  decay products into HEPEVT.  Returns LSUCCS, channel code, daughter !
!  level.                                                              !
!======================================================================!
      SUBROUTINE chromo_dcy_sample(ia, iz, im, lsuccs_out,
     &                             kdcy_out, ilv_out)
      INCLUDE '(DBLPRW)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IOUNIT)'

      INCLUDE '(DCYFLG)'
      INCLUDE '(EMFSTK)'
      INCLUDE '(FHEAVY)'
      INCLUDE '(GENSTK)'

      INTEGER ia, iz, im, lsuccs_out, kdcy_out, ilv_out
      DOUBLE PRECISION tdelay
      LOGICAL lsuccs
Cf2py intent(in)  ia, iz, im
Cf2py intent(out) lsuccs_out, kdcy_out, ilv_out

      lsuccs_out = 0
      kdcy_out = -1
      ilv_out = 0

      np     = 0
      np0    = 0
      npheav = 0
      npemf  = 0
      lsuccs = .TRUE.
      tdelay = -1.0D+09

      CALL SPDCEV(ia, iz, im, ZERZER, ZERZER, ZERZER, ZERZER,
     &            ONEONE, ONEONE, tdelay, .TRUE., .TRUE., lsuccs)

      IF (.NOT. lsuccs) RETURN
      lsuccs_out = 1
      ilv_out = ILVDCY

      IF (LALDCY) THEN
         kdcy_out = 1
      ELSE IF (LBMDCY) THEN
         kdcy_out = 2
      ELSE IF (LBPDCY) THEN
         kdcy_out = 3
      ELSE IF (LECDCY) THEN
         kdcy_out = 4
      ELSE IF (LITDCY) THEN
         kdcy_out = 5
      ELSE IF (LSFDCY) THEN
         kdcy_out = 6
      END IF

      CALL chromo_fllhep
      RETURN
      END SUBROUTINE chromo_dcy_sample
```

- [ ] **Step 4: Add to meson.build**

Add `'chromo_dcy_sample'` to `fluka_syms`.

- [ ] **Step 5: Rebuild**

```bash
pip install --no-build-isolation -v -e .[test] 2>&1 | tail -5
```

- [ ] **Step 6: Run tests**

```bash
python -m pytest tests/test_fluka_decay.py -v
```
Expected: 8 PASS.

- [ ] **Step 7: Commit**

```bash
git add src/fortran/fluka/chromo_fluka.f meson.build tests/test_fluka_decay.py
git commit -m "feat(fluka): add chromo_dcy_sample wrapper (SPDCEV+FLLHEP)"
```

---

## Task 7: Python module scaffold + dataclasses

**Files:**
- Create: `src/chromo/models/fluka_decay.py`
- Modify: `tests/test_fluka_decay.py`

- [ ] **Step 1: Write failing test**

Append to `tests/test_fluka_decay.py`:

```python
def _scaffold_imports():
    from chromo.models.fluka_decay import (
        DecayChannel,
        DecayLine,
        FlukaIsotope,
        FlukaDecay,
        STABLE_DEFAULT,
    )
    return True


def test_module_scaffold_imports():
    assert run_in_separate_process(_scaffold_imports) is True


def _dataclass_smoke():
    from chromo.models.fluka_decay import DecayChannel, DecayLine

    ch = DecayChannel(name="B-", branching=1.0,
                      daughter_A=137, daughter_Z=56, daughter_m=1,
                      q_value=0.514)
    ln = DecayLine(energy=0.66166, branching=0.85,
                   end_level=0, is_positron=False)
    return ch.name, ch.daughter_A, ln.energy


def test_dataclass_smoke():
    name, a, e = run_in_separate_process(_dataclass_smoke)
    assert name == "B-"
    assert a == 137
    assert abs(e - 0.66166) < 1e-6
```

- [ ] **Step 2: Run, verify failure**

```bash
python -m pytest tests/test_fluka_decay.py::test_module_scaffold_imports -v
```
Expected: FAIL (`ModuleNotFoundError`).

- [ ] **Step 3: Create `fluka_decay.py`**

Create `src/chromo/models/fluka_decay.py`:

```python
"""FLUKA radioactive-decay interface.

Exposes FLUKA 2025.1's nuclear-decay tables and `SPDCEV` event sampler
to chromo users via a separate, kinematics-free `FlukaDecay` class.
See ``docs/superpowers/specs/2026-04-30-fluka-decay-interface-design.md``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

# --------------------------------------------------------------------- #
# Data containers                                                       #
# --------------------------------------------------------------------- #

# FLUKA channel-kind codes returned by chromo_dcy_channels / chromo_dcy_sample.
# Keep in sync with the Fortran wrapper in chromo_fluka.f.
_CHANNEL_NAMES = {
    1: "alpha",
    2: "B-",
    3: "B+",
    4: "EC",
    5: "IT",
    6: "SF",
    7: "B-N",
    8: "B+P",
    9: "B-2N",
    10: "B-3N",
    11: "B-NA",
    12: "other",
}


@dataclass(frozen=True, slots=True)
class DecayChannel:
    """One radioactive-decay channel from FLUKA's table.

    Attributes
    ----------
    name : str
        Channel label (``"B-"``, ``"alpha"``, ``"EC"``, ...).
    branching : float
        Branching ratio in [0, 1].
    daughter_A, daughter_Z, daughter_m : int
        Daughter nuclide. ``-1`` if the channel has no single daughter
        (e.g. spontaneous fission).
    q_value : float
        Q value in MeV (atomic, computed via FLUKA's ``QRDDCY``).
    """

    name: str
    branching: float
    daughter_A: int
    daughter_Z: int
    daughter_m: int
    q_value: float


@dataclass(frozen=True, slots=True)
class DecayLine:
    """One discrete line from FLUKA's gamma/alpha/CE/beta tables.

    Attributes
    ----------
    energy : float
        Line energy in MeV. For β±, this is the end-point.
    branching : float
        Branching ratio in [0, 1] of the parent decay.
    end_level : int
        For α / β±, daughter level index reached. ``0`` for γ / CE.
    is_positron : bool
        For β± only — true if positron, false if electron.
    """

    energy: float
    branching: float
    end_level: int
    is_positron: bool


# Default "stop the chain here" set for `FlukaDecay.chain`. Populated on
# first FlukaDecay instantiation by `_compute_stable_default`.
STABLE_DEFAULT: set[int] = set()


# --------------------------------------------------------------------- #
# FlukaIsotope and FlukaDecay placeholders — filled in later tasks      #
# --------------------------------------------------------------------- #


class FlukaIsotope:
    """Single nuclide record from FLUKA's decay catalogue.

    Filled out across Tasks 8, 11, 12.
    """

    __slots__ = (
        "A",
        "Z",
        "m",
        "t_half",
        "mass_excess",
        "symbol",
        "j_spin",
        "j_parity",
        "_channels",
        "_lines",
        "_owner",
    )


class FlukaDecay:
    """FLUKA radioactive-decay generator and database.

    Filled out across Tasks 9, 10, 13, 14, 15.
    """
```

- [ ] **Step 4: Run tests**

```bash
python -m pytest tests/test_fluka_decay.py::test_module_scaffold_imports tests/test_fluka_decay.py::test_dataclass_smoke -v
```
Expected: 2 PASS.

- [ ] **Step 5: Commit**

```bash
git add src/chromo/models/fluka_decay.py tests/test_fluka_decay.py
git commit -m "feat(fluka): add fluka_decay module scaffold + dataclasses"
```

---

## Task 8: `FlukaIsotope.__init__` and `short()` repr

**Files:**
- Modify: `src/chromo/models/fluka_decay.py`
- Modify: `tests/test_fluka_decay.py`

- [ ] **Step 1: Write failing test**

Append to `tests/test_fluka_decay.py`:

```python
def _make_isotope():
    from chromo.models.fluka_decay import FlukaIsotope

    iso = FlukaIsotope(
        owner=None,
        A=137, Z=55, m=0,
        t_half=9.49e8, mass_excess=-86.546,
        symbol="Cs", j_spin=7, j_parity=1,
    )
    return (iso.A, iso.Z, iso.m, iso.symbol, iso.short())


def test_isotope_init_and_short():
    A, Z, m, sym, short = run_in_separate_process(_make_isotope)
    assert A == 137 and Z == 55 and m == 0 and sym == "Cs"
    assert "Cs137" in short
    assert "9.49" in short or "9.490" in short
```

- [ ] **Step 2: Run, verify failure**

```bash
python -m pytest tests/test_fluka_decay.py::test_isotope_init_and_short -v
```
Expected: FAIL (`FlukaIsotope` has no `__init__`).

- [ ] **Step 3: Implement `FlukaIsotope.__init__` + `short()`**

Replace the placeholder `class FlukaIsotope:` block in
`src/chromo/models/fluka_decay.py` with:

```python
class FlukaIsotope:
    """Single nuclide record from FLUKA's decay catalogue.

    Heavy data (decay channels, line lists) is fetched lazily on first
    access via the parent ``FlukaDecay`` instance.
    """

    __slots__ = (
        "A",
        "Z",
        "m",
        "t_half",
        "mass_excess",
        "symbol",
        "j_spin",
        "j_parity",
        "_channels",
        "_lines",
        "_owner",
    )

    def __init__(self, *, owner, A, Z, m, t_half, mass_excess,
                 symbol, j_spin, j_parity):
        self._owner = owner
        self.A = A
        self.Z = Z
        self.m = m
        self.t_half = t_half
        self.mass_excess = mass_excess
        self.symbol = symbol
        self.j_spin = j_spin
        self.j_parity = j_parity
        self._channels = None
        self._lines = {}

    def short(self) -> str:
        """One-line summary."""
        m_tag = "" if self.m == 0 else f"m{self.m}"
        return (
            f"{self.symbol}{self.A}{m_tag} (Z={self.Z}, m={self.m}): "
            f"T1/2={self.t_half:.3e} s, ExM={self.mass_excess:.3f} MeV"
        )

    def __repr__(self):
        return f"<FlukaIsotope {self.short()}>"
```

- [ ] **Step 4: Run tests**

```bash
python -m pytest tests/test_fluka_decay.py::test_isotope_init_and_short -v
```
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/chromo/models/fluka_decay.py tests/test_fluka_decay.py
git commit -m "feat(fluka): FlukaIsotope dataclass with short() repr"
```

---

## Task 9: `FlukaDecay.__init__` + init guard + `lookup`

**Files:**
- Modify: `src/chromo/models/fluka_decay.py`
- Modify: `tests/test_fluka_decay.py`

- [ ] **Step 1: Write failing tests**

Append to `tests/test_fluka_decay.py`:

```python
def _construct_and_lookup():
    from chromo.models.fluka_decay import FlukaDecay

    dcy = FlukaDecay(seed=42)
    iso = dcy.lookup("Cs137")
    return (iso.A, iso.Z, iso.m, iso.symbol,
            round(iso.t_half / 1e8, 1))


def test_fluka_decay_construct_and_lookup_cs137():
    A, Z, m, sym, t12_e8 = run_in_separate_process(_construct_and_lookup)
    assert (A, Z, m) == (137, 55, 0)
    assert sym == "Cs"
    assert 9.0 < t12_e8 < 10.0


def _lookup_isomer():
    from chromo.models.fluka_decay import FlukaDecay

    dcy = FlukaDecay()
    iso = dcy.lookup("Tc99m")
    return iso.A, iso.Z, iso.m, round(iso.t_half, 0)


def test_fluka_decay_lookup_isomer():
    A, Z, m, t12 = run_in_separate_process(_lookup_isomer)
    assert (A, Z, m) == (99, 43, 1)
    assert 20000 <= t12 <= 23000  # ~ 6.0 hr ≈ 21600 s


def _lookup_tuple_and_missing():
    from chromo.models.fluka_decay import FlukaDecay

    dcy = FlukaDecay()
    iso_t = dcy.lookup(238, 92, 0)
    iso_none = dcy.lookup(2, 50, 0)
    return iso_t.A, iso_t.Z, iso_none


def test_fluka_decay_lookup_tuple_and_missing():
    A, Z, none_ = run_in_separate_process(_lookup_tuple_and_missing)
    assert (A, Z) == (238, 92)
    assert none_ is None
```

- [ ] **Step 2: Run, verify failure**

```bash
python -m pytest tests/test_fluka_decay.py::test_fluka_decay_construct_and_lookup_cs137 -v
```
Expected: FAIL.

- [ ] **Step 3: Implement `FlukaDecay` init + lookup**

Replace the placeholder `class FlukaDecay:` block in
`src/chromo/models/fluka_decay.py` with:

```python
import re

from particle import Particle


# Module-level guard — shared with chromo.models.fluka.Fluka via
# `_mark_fluka_init_done()` (set in fluka.py after STPXYZ).
_INITIALIZED = False
_INSTANCE_COUNT = 0
_NAME_RE = re.compile(r"^([A-Za-z]{1,2})(\d{1,3})(?:m(\d*))?$")


def _mark_fluka_init_done():
    """Called by Fluka.__init__ after STPXYZ to share the init state."""
    global _INITIALIZED
    _INITIALIZED = True


def _ensure_init(lib):
    global _INITIALIZED
    if _INITIALIZED:
        return
    lib.chromo_dcy_init()
    _INITIALIZED = True


def _z_to_symbol(z: int) -> str:
    """Element symbol for atomic number Z."""
    try:
        return Particle.from_pdgid(1_000_000_000 + z * 10_000).name[:2].rstrip("0123456789")
    except Exception:
        return _Z_FALLBACK.get(z, "?")


# Compact fallback Z->symbol table for the cases where `particle` does
# not know the nuclide (unbound, theoretical, neutron):
_Z_FALLBACK = {
    0: "n", 1: "H", 2: "He", 3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N",
    8: "O", 9: "F", 10: "Ne", 11: "Na", 12: "Mg", 13: "Al", 14: "Si",
    15: "P", 16: "S", 17: "Cl", 18: "Ar", 19: "K", 20: "Ca", 21: "Sc",
    22: "Ti", 23: "V", 24: "Cr", 25: "Mn", 26: "Fe", 27: "Co", 28: "Ni",
    29: "Cu", 30: "Zn", 31: "Ga", 32: "Ge", 33: "As", 34: "Se", 35: "Br",
    36: "Kr", 37: "Rb", 38: "Sr", 39: "Y", 40: "Zr", 41: "Nb", 42: "Mo",
    43: "Tc", 44: "Ru", 45: "Rh", 46: "Pd", 47: "Ag", 48: "Cd", 49: "In",
    50: "Sn", 51: "Sb", 52: "Te", 53: "I", 54: "Xe", 55: "Cs", 56: "Ba",
    57: "La", 58: "Ce", 59: "Pr", 60: "Nd", 61: "Pm", 62: "Sm", 63: "Eu",
    64: "Gd", 65: "Tb", 66: "Dy", 67: "Ho", 68: "Er", 69: "Tm", 70: "Yb",
    71: "Lu", 72: "Hf", 73: "Ta", 74: "W", 75: "Re", 76: "Os", 77: "Ir",
    78: "Pt", 79: "Au", 80: "Hg", 81: "Tl", 82: "Pb", 83: "Bi", 84: "Po",
    85: "At", 86: "Rn", 87: "Fr", 88: "Ra", 89: "Ac", 90: "Th", 91: "Pa",
    92: "U", 93: "Np", 94: "Pu", 95: "Am", 96: "Cm", 97: "Bk", 98: "Cf",
    99: "Es", 100: "Fm", 101: "Md", 102: "No", 103: "Lr", 104: "Rf",
    105: "Db", 106: "Sg", 107: "Bh", 108: "Hs", 109: "Mt", 110: "Ds",
}
_SYM_TO_Z = {v: k for k, v in _Z_FALLBACK.items()}


def _parse_isotope_name(name: str) -> tuple[int, int, int]:
    """Parse 'Cs137', 'Tc99m', 'U238m1' into (A, Z, m)."""
    match = _NAME_RE.match(name.strip())
    if not match:
        raise ValueError(f"Cannot parse isotope name: {name!r}")
    sym, a_str, m_str = match.groups()
    sym = sym[0].upper() + (sym[1:].lower() if len(sym) > 1 else "")
    if sym not in _SYM_TO_Z:
        raise ValueError(f"Unknown element symbol: {sym!r}")
    A = int(a_str)
    Z = _SYM_TO_Z[sym]
    if m_str is None:
        m = 0
    elif m_str == "":
        m = 1
    else:
        m = int(m_str)
    return A, Z, m


class FlukaDecay:
    """FLUKA radioactive-decay generator and database.

    No kinematics are required.  The first instance triggers FLUKA's
    decay-table init (``CMSPPR`` / ``ZEROIN`` / ``NCDTRD`` / ``RDFLUO`` /
    ...).  If a ``Fluka`` instance has already initialised the library
    (via ``STPXYZ``), construction is a no-op.

    Two ``FlukaDecay`` instances in the same process raise
    ``RuntimeError`` (matches the FLUKA singleton).
    """

    def __init__(self, seed: int | None = None):
        global _INSTANCE_COUNT
        if _INSTANCE_COUNT > 0:
            raise RuntimeError(
                "FlukaDecay can only be instantiated once per process "
                "(FLUKA singleton)."
            )
        _INSTANCE_COUNT += 1
        from chromo.models import _fluka

        self._lib = _fluka
        _ensure_init(self._lib)
        if seed is not None:
            self._lib.init_rng_state(int(seed))
        self._catalog = None  # filled in Task 10

    # -- lookup ---------------------------------------------------------

    def lookup(self, *args) -> FlukaIsotope | None:
        """Return the FlukaIsotope for the given (A, Z, m) or name.

        Accepts:

        - ``lookup("Cs137")`` / ``lookup("Tc99m")`` / ``lookup("U238m1")``
        - ``lookup(A, Z, m)``
        - ``lookup(pdg_ion_id)``
        """
        A, Z, m = self._parse_arg(args)
        found, t12, exm, jsp, jpt = self._lib.chromo_dcy_lookup(A, Z, m)
        if not found:
            return None
        return FlukaIsotope(
            owner=self,
            A=A,
            Z=Z,
            m=m,
            t_half=float(t12),
            mass_excess=float(exm),
            symbol=_z_to_symbol(Z) if Z >= 0 else "n",
            j_spin=int(jsp),
            j_parity=int(jpt),
        )

    @staticmethod
    def _parse_arg(args) -> tuple[int, int, int]:
        if len(args) == 1:
            arg = args[0]
            if isinstance(arg, str):
                return _parse_isotope_name(arg)
            if isinstance(arg, int):
                # PDG ion code: 10LZZZAAAI
                if arg >= 1_000_000_000:
                    A = (arg // 10) % 1000
                    Z = (arg // 10_000) % 1000
                    m = (arg // 100_000_000) % 10
                    return A, Z, m
                raise ValueError(f"Cannot parse single-int arg {arg}")
        if len(args) == 3:
            return int(args[0]), int(args[1]), int(args[2])
        raise TypeError(
            "lookup() expects (name) or (A, Z, m); got " + repr(args)
        )
```

For the `_z_to_symbol` helper, simplify since the fallback table covers Z=0..110:

Replace the `_z_to_symbol` function body with just:
```python
def _z_to_symbol(z: int) -> str:
    """Element symbol for atomic number Z (or 'n' for free neutron)."""
    return _Z_FALLBACK.get(z, "?")
```

- [ ] **Step 4: Run tests**

```bash
python -m pytest tests/test_fluka_decay.py::test_fluka_decay_construct_and_lookup_cs137 tests/test_fluka_decay.py::test_fluka_decay_lookup_isomer tests/test_fluka_decay.py::test_fluka_decay_lookup_tuple_and_missing -v
```
Expected: 3 PASS.

- [ ] **Step 5: Commit**

```bash
git add src/chromo/models/fluka_decay.py tests/test_fluka_decay.py
git commit -m "feat(fluka): FlukaDecay.__init__ and lookup()"
```

---

## Task 10: `FlukaDecay.catalog()` + filters

**Files:**
- Modify: `src/chromo/models/fluka_decay.py`
- Modify: `tests/test_fluka_decay.py`

- [ ] **Step 1: Write failing tests**

Append to `tests/test_fluka_decay.py`:

```python
def _catalog_full_and_filter():
    from chromo.models.fluka_decay import FlukaDecay

    dcy = FlukaDecay()
    full = dcy.catalog()
    long_lived = dcy.catalog(t_half_min=1e10)
    short_only = dcy.catalog(t_half_max=1.0, t_half_min=1e-3)
    actinides = dcy.catalog(z_min=89, z_max=99)
    return len(full), len(long_lived), len(short_only), len(actinides)


def test_catalog_filters():
    n, n_long, n_short, n_act = run_in_separate_process(
        _catalog_full_and_filter
    )
    assert 4000 <= n <= 5500
    assert 200 <= n_long <= 600         # plenty of T1/2 > 1e10 s entries
    assert n_short > 0
    assert n_act > 0


def _catalog_contains_known():
    from chromo.models.fluka_decay import FlukaDecay

    dcy = FlukaDecay()
    full = dcy.catalog()
    azm = {(i.A, i.Z, i.m) for i in full}
    return ((238, 92, 0) in azm, (137, 55, 0) in azm,
            (14, 6, 0) in azm, (99, 43, 1) in azm)


def test_catalog_contains_known_isotopes():
    found = run_in_separate_process(_catalog_contains_known)
    assert all(found), found
```

- [ ] **Step 2: Run, verify failure**

```bash
python -m pytest tests/test_fluka_decay.py::test_catalog_filters -v
```
Expected: FAIL (`catalog()` not implemented).

- [ ] **Step 3: Implement `catalog()`**

Add inside `class FlukaDecay:` in `src/chromo/models/fluka_decay.py` (after `lookup`):

```python
    # -- catalogue ------------------------------------------------------

    def catalog(
        self,
        *,
        t_half_min: float | None = None,
        t_half_max: float | None = None,
        a_min: int | None = None,
        a_max: int | None = None,
        z_min: int | None = None,
        z_max: int | None = None,
    ) -> list[FlukaIsotope]:
        """Return all isotopes/isomers in FLUKA's decay-data table.

        First call materialises the full catalogue (~4 500 entries) by
        calling ``chromo_dcy_catalog``; subsequent calls reuse the cache.

        Parameters
        ----------
        t_half_min, t_half_max : float, optional
            Half-life range in seconds (inclusive).
        a_min, a_max, z_min, z_max : int, optional
            Mass / atomic-number ranges (inclusive).
        """
        if self._catalog is None:
            self._catalog = self._fetch_full_catalog()

        out = self._catalog
        if t_half_min is not None:
            out = [i for i in out if i.t_half >= t_half_min]
        if t_half_max is not None:
            out = [i for i in out if i.t_half <= t_half_max]
        if a_min is not None:
            out = [i for i in out if i.A >= a_min]
        if a_max is not None:
            out = [i for i in out if i.A <= a_max]
        if z_min is not None:
            out = [i for i in out if i.Z >= z_min]
        if z_max is not None:
            out = [i for i in out if i.Z <= z_max]
        return list(out)

    def _fetch_full_catalog(self) -> list[FlukaIsotope]:
        import numpy as np

        max_n = 5500
        a = np.zeros(max_n, dtype=np.int32)
        z = np.zeros(max_n, dtype=np.int32)
        m = np.zeros(max_n, dtype=np.int32)
        t12 = np.zeros(max_n, dtype=np.float64)
        exm = np.zeros(max_n, dtype=np.float64)
        jsp = np.zeros(max_n, dtype=np.int32)
        jpt = np.zeros(max_n, dtype=np.int32)
        n = self._lib.chromo_dcy_catalog(max_n, a, z, m,
                                         t12, exm, jsp, jpt)
        return [
            FlukaIsotope(
                owner=self,
                A=int(a[i]),
                Z=int(z[i]),
                m=int(m[i]),
                t_half=float(t12[i]),
                mass_excess=float(exm[i]),
                symbol=_z_to_symbol(int(z[i])),
                j_spin=int(jsp[i]),
                j_parity=int(jpt[i]),
            )
            for i in range(int(n))
        ]
```

- [ ] **Step 4: Run tests**

```bash
python -m pytest tests/test_fluka_decay.py::test_catalog_filters tests/test_fluka_decay.py::test_catalog_contains_known_isotopes -v
```
Expected: 2 PASS.

- [ ] **Step 5: Commit**

```bash
git add src/chromo/models/fluka_decay.py tests/test_fluka_decay.py
git commit -m "feat(fluka): FlukaDecay.catalog() with filters"
```

---

## Task 11: Lazy `FlukaIsotope.channels`

**Files:**
- Modify: `src/chromo/models/fluka_decay.py`
- Modify: `tests/test_fluka_decay.py`

- [ ] **Step 1: Write failing test**

Append to `tests/test_fluka_decay.py`:

```python
def _cs137_channels():
    from chromo.models.fluka_decay import FlukaDecay

    dcy = FlukaDecay()
    iso = dcy.lookup("Cs137")
    chs = iso.channels
    return len(chs), chs[0].name, chs[0].daughter_A, chs[0].daughter_Z


def test_isotope_channels_cs137():
    n, name, dA, dZ = run_in_separate_process(_cs137_channels)
    assert n >= 1
    assert name == "B-"
    assert (dA, dZ) == (137, 56)


def _u238_alpha():
    from chromo.models.fluka_decay import FlukaDecay

    dcy = FlukaDecay()
    iso = dcy.lookup("U238")
    return [(c.name, c.branching, c.daughter_A, c.daughter_Z)
            for c in iso.channels]


def test_isotope_channels_u238_has_alpha():
    chs = run_in_separate_process(_u238_alpha)
    assert any(c[0] == "alpha" and c[2] == 234 and c[3] == 90 for c in chs)
```

- [ ] **Step 2: Run, verify failure**

```bash
python -m pytest tests/test_fluka_decay.py::test_isotope_channels_cs137 -v
```
Expected: FAIL.

- [ ] **Step 3: Add `channels` property to `FlukaIsotope`**

Add inside `class FlukaIsotope:` in `src/chromo/models/fluka_decay.py`:

```python
    @property
    def channels(self) -> tuple[DecayChannel, ...]:
        """Decay channels (lazy fetch on first access)."""
        if self._channels is None:
            self._channels = self._owner._fetch_channels(
                self.A, self.Z, self.m
            )
        return self._channels
```

Add `_fetch_channels` to `FlukaDecay`:

```python
    # -- lazy backends used by FlukaIsotope -----------------------------

    def _fetch_channels(self, A: int, Z: int, m: int) -> tuple[DecayChannel, ...]:
        import numpy as np

        max_ch = 8
        kind = np.zeros(max_ch, dtype=np.int32)
        br = np.zeros(max_ch, dtype=np.float64)
        da = np.zeros(max_ch, dtype=np.int32)
        dz = np.zeros(max_ch, dtype=np.int32)
        dm = np.zeros(max_ch, dtype=np.int32)
        qv = np.zeros(max_ch, dtype=np.float64)
        n = self._lib.chromo_dcy_channels(
            A, Z, m, max_ch, kind, br, da, dz, dm, qv
        )
        return tuple(
            DecayChannel(
                name=_CHANNEL_NAMES.get(int(kind[i]), "other"),
                branching=float(br[i]),
                daughter_A=int(da[i]),
                daughter_Z=int(dz[i]),
                daughter_m=int(dm[i]),
                q_value=float(qv[i]),
            )
            for i in range(int(n))
        )
```

- [ ] **Step 4: Run tests**

```bash
python -m pytest tests/test_fluka_decay.py::test_isotope_channels_cs137 tests/test_fluka_decay.py::test_isotope_channels_u238_has_alpha -v
```
Expected: 2 PASS.

- [ ] **Step 5: Commit**

```bash
git add src/chromo/models/fluka_decay.py tests/test_fluka_decay.py
git commit -m "feat(fluka): FlukaIsotope.channels lazy property"
```

---

## Task 12: Lazy line-list properties + `FlukaIsotope.__str__`

**Files:**
- Modify: `src/chromo/models/fluka_decay.py`
- Modify: `tests/test_fluka_decay.py`

- [ ] **Step 1: Write failing test**

Append to `tests/test_fluka_decay.py`:

```python
def _cs137_gamma_and_str():
    from chromo.models.fluka_decay import FlukaDecay

    dcy = FlukaDecay()
    iso = dcy.lookup("Cs137")
    gammas = iso.gamma_lines
    s = str(iso)
    return [(round(g.energy, 5), round(g.branching, 4))
            for g in gammas], s


def test_gamma_lines_and_str():
    gammas, s = run_in_separate_process(_cs137_gamma_and_str)
    energies = [e for e, _b in gammas]
    assert any(abs(e - 0.66166) < 0.001 for e in energies), gammas
    # str(iso) should contain identifier + gamma energy
    assert "Cs137" in s
    assert "0.66" in s or "0.6617" in s
```

- [ ] **Step 2: Run, verify failure**

```bash
python -m pytest tests/test_fluka_decay.py::test_gamma_lines_and_str -v
```
Expected: FAIL (`gamma_lines` not implemented).

- [ ] **Step 3: Add line-list properties + `__str__`**

Add inside `class FlukaIsotope:` in `src/chromo/models/fluka_decay.py`:

```python
    def _lines(self, kind_int: int, kind_str: str) -> tuple[DecayLine, ...]:
        if kind_str in self._lines:
            return self._lines[kind_str]
        result = self._owner._fetch_lines(
            self.A, self.Z, self.m, kind_int
        )
        self._lines[kind_str] = result
        return result

    @property
    def gamma_lines(self) -> tuple[DecayLine, ...]:
        return self._lines(1, "gamma")

    @property
    def alpha_lines(self) -> tuple[DecayLine, ...]:
        return self._lines(2, "alpha")

    @property
    def ce_lines(self) -> tuple[DecayLine, ...]:
        return self._lines(3, "ce")

    @property
    def beta_spectra(self) -> tuple[DecayLine, ...]:
        return self._lines(4, "beta")

    def __str__(self) -> str:
        m_tag = "" if self.m == 0 else f"m{self.m}"
        head = (
            f"Isotope {self.symbol}{self.A}{m_tag}  "
            f"(A={self.A}, Z={self.Z}, m={self.m})\n"
            f"  T1/2     = {self.t_half:.3e} s\n"
            f"  ExMass   = {self.mass_excess:.4f} MeV\n"
            f"  J        = {self.j_spin / 2:.1f}  parity = "
            f"{'+' if self.j_parity > 0 else '-' if self.j_parity < 0 else '?'}"
        )
        rows = []
        for c in self.channels:
            d = (
                f"-> {c.daughter_A:3d}/{c.daughter_Z:3d}/m{c.daughter_m}"
                if c.daughter_A >= 0
                else "-> (no single daughter)"
            )
            rows.append(
                f"    {c.name:<5s}  BR={c.branching*100:7.3f}%  "
                f"{d}  Q={c.q_value:.4f} MeV"
            )
        chan = "\n  Channels:\n" + "\n".join(rows) if rows else ""

        def _fmt_lines(label, lines, n_max=10):
            if not lines:
                return ""
            out = [f"\n  {label} ({len(lines)}):"]
            for ln in lines[:n_max]:
                out.append(
                    f"    E={ln.energy:9.5f} MeV  "
                    f"BR={ln.branching*100:7.3f}%"
                )
            if len(lines) > n_max:
                out.append(f"    ... ({len(lines) - n_max} more)")
            return "\n".join(out)

        body = (
            chan
            + _fmt_lines("Gamma lines",  self.gamma_lines)
            + _fmt_lines("Alpha lines",  self.alpha_lines)
            + _fmt_lines("CE/Auger lines", self.ce_lines)
            + _fmt_lines("Beta+/- spectra", self.beta_spectra)
        )
        return head + body
```

Add `_fetch_lines` to `FlukaDecay`:

```python
    def _fetch_lines(self, A: int, Z: int, m: int,
                     kind: int) -> tuple[DecayLine, ...]:
        import numpy as np

        max_l = 256
        br = np.zeros(max_l, dtype=np.float64)
        e_mev = np.zeros(max_l, dtype=np.float64)
        nlev = np.zeros(max_l, dtype=np.int32)
        n = self._lib.chromo_dcy_lines(A, Z, m, kind, max_l,
                                       br, e_mev, nlev)
        return tuple(
            DecayLine(
                energy=float(e_mev[i]),
                branching=float(br[i]),
                end_level=int(nlev[i]),
                is_positron=(kind == 4 and float(e_mev[i]) > 0),
            )
            for i in range(int(n))
        )
```

- [ ] **Step 4: Run tests**

```bash
python -m pytest tests/test_fluka_decay.py::test_gamma_lines_and_str -v
```
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/chromo/models/fluka_decay.py tests/test_fluka_decay.py
git commit -m "feat(fluka): FlukaIsotope line-list properties + __str__"
```

---

## Task 13: `FlukaDecay.__call__` — inclusive sampling

**Files:**
- Modify: `src/chromo/models/fluka_decay.py`
- Modify: `tests/test_fluka_decay.py`

- [ ] **Step 1: Write failing test**

Append to `tests/test_fluka_decay.py`:

```python
def _sample_cs137_1000():
    from chromo.models.fluka_decay import FlukaDecay

    dcy = FlukaDecay(seed=42)
    n_events = 0
    n_with_electron = 0
    np_total = 0
    for ev in dcy("Cs137", n=1000):
        n_events += 1
        np_total += len(ev.pid)
        if 11 in ev.pid:               # PDG 11 = e-
            n_with_electron += 1
    return n_events, n_with_electron, np_total


def test_sample_cs137_1000_events():
    n, n_e, np_total = run_in_separate_process(_sample_cs137_1000)
    assert n == 1000
    # Cs-137 always β-: every event has at least one e-
    assert n_e == 1000
    # Mean NP ≈ 2.93 from prototype; allow generous range
    assert 2.0 < np_total / n < 4.0


def _sample_stable_raises():
    from chromo.models.fluka_decay import FlukaDecay

    dcy = FlukaDecay()
    try:
        list(dcy("Fe56", n=1))
    except ValueError:
        return "raised"
    return "did_not_raise"


def test_sample_stable_raises():
    assert run_in_separate_process(_sample_stable_raises) == "raised"
```

- [ ] **Step 2: Run, verify failure**

```bash
python -m pytest tests/test_fluka_decay.py::test_sample_cs137_1000_events -v
```
Expected: FAIL.

- [ ] **Step 3: Implement `__call__`**

Add inside `class FlukaDecay:`:

```python
    # -- inclusive event sampling --------------------------------------

    _STABLE_T_HALF = 1e38

    def __call__(self, parent, n: int):
        """Yield ``n`` correlated decay events for ``parent``.

        Parameters
        ----------
        parent : str | tuple | int
            Same forms accepted by ``lookup``.
        n : int
            Number of events to generate.

        Yields
        ------
        EventData
            One event per yield, with FLUKA's products in
            ``ev.pid / px / py / pz / en / mass``.
        """
        iso = self.lookup(*([parent] if not isinstance(parent, tuple)
                            else parent))
        if iso is None:
            raise ValueError(
                f"FLUKA has no decay data for parent {parent!r}."
            )
        if iso.t_half >= self._STABLE_T_HALF:
            raise ValueError(
                f"{parent!r} is stable in FLUKA's table "
                f"(T1/2 = {iso.t_half:.3e} s); cannot decay-sample."
            )

        for _ in range(int(n)):
            yield self._sample_one(iso.A, iso.Z, iso.m)

    def _sample_one(self, A: int, Z: int, m: int):
        """One SPDCEV call → one EventData via FlukaEvent extraction."""
        from chromo.models.fluka import FlukaEvent

        for _attempt in range(2):
            ok, _kdcy, _ilv = self._lib.chromo_dcy_sample(A, Z, m)
            if ok:
                return FlukaEvent(self._lib).copy()
        # Persistent failure for this isotope: emit empty event with
        # info log so the caller can detect/skip.
        from chromo.util import info as _info

        _info(0, f"chromo_dcy_sample failed for ({A},{Z},{m}) twice")
        return FlukaEvent(self._lib).copy()
```

Note: ``FlukaEvent`` already exists in ``chromo/models/fluka.py``; its
constructor reads NP / KPART / TKI / direction cosines from the FLUKA
common blocks via `chromo_fllhep` having populated HEPEVT.

- [ ] **Step 4: Run tests**

```bash
python -m pytest tests/test_fluka_decay.py::test_sample_cs137_1000_events tests/test_fluka_decay.py::test_sample_stable_raises -v
```
Expected: 2 PASS.

- [ ] **Step 5: Commit**

```bash
git add src/chromo/models/fluka_decay.py tests/test_fluka_decay.py
git commit -m "feat(fluka): FlukaDecay.__call__ inclusive sampling"
```

---

## Task 14: `STABLE_DEFAULT` set + `DecayChainHandler.expand`

**Files:**
- Modify: `src/chromo/models/fluka_decay.py`
- Modify: `tests/test_fluka_decay.py`

- [ ] **Step 1: Write failing test**

Append to `tests/test_fluka_decay.py`:

```python
def _stable_default_populated():
    from chromo.models.fluka_decay import FlukaDecay, STABLE_DEFAULT

    FlukaDecay()  # construction populates STABLE_DEFAULT
    return len(STABLE_DEFAULT), 22 in STABLE_DEFAULT  # 22 = photon


def test_stable_default_populated():
    n, has_photon = run_in_separate_process(_stable_default_populated)
    assert n > 200       # leptons + photons + ~250 stable nuclei
    assert has_photon


def _chain_cs137_handler():
    from chromo.models.fluka_decay import FlukaDecay, DecayChainHandler

    dcy = FlukaDecay(seed=7)
    handler = DecayChainHandler(dcy)
    n_events = 0
    n_with_2_gamma = 0
    for ev in dcy("Cs137", n=200):
        chained = handler.expand(ev)
        n_events += 1
        n_g = sum(1 for p in chained.pid if p == 22)
        if n_g >= 1:
            n_with_2_gamma += 1
    return n_events, n_with_2_gamma


def test_chain_handler_cs137_includes_gamma():
    n, n_g = run_in_separate_process(_chain_cs137_handler)
    assert n == 200
    # Cs137->Ba137m1->Ba137 + 0.661 MeV γ in ~94.7% of the chain decays.
    assert n_g >= 150
```

- [ ] **Step 2: Run, verify failure**

```bash
python -m pytest tests/test_fluka_decay.py::test_stable_default_populated -v
```
Expected: FAIL.

- [ ] **Step 3: Implement `STABLE_DEFAULT` and `DecayChainHandler`**

Add to `src/chromo/models/fluka_decay.py` (append at module bottom):

```python
# --------------------------------------------------------------------- #
# Stable-default set + DecayChainHandler                                #
# --------------------------------------------------------------------- #

# Lepton/photon/light-nucleon PDG ids that anchor the chain.
_BASE_STABLE_PDG = {
    11, -11,        # e-, e+
    12, -12,        # nu_e, anu_e
    13, -13,        # mu-, mu+
    14, -14,        # nu_mu, anu_mu
    16, -16,        # nu_tau, anu_tau
    22,             # gamma
    2212,           # proton
    2112,           # neutron
}


def _populate_stable_default(catalog) -> None:
    """Add to STABLE_DEFAULT every (A,Z,m) entry with T1/2 = 1e38."""
    STABLE_DEFAULT.clear()
    STABLE_DEFAULT.update(_BASE_STABLE_PDG)
    for iso in catalog:
        if iso.t_half >= 1e38:
            pdg = 1_000_000_000 + iso.Z * 10_000 + iso.A * 10 + iso.m
            STABLE_DEFAULT.add(pdg)


class DecayChainHandler:
    """Recursive decay-chain post-processor.

    For each event passed to ``expand()``, every product whose PDG id is
    *not* in ``final_state`` *and* corresponds to a FLUKA-decayable
    nuclide is sampled via ``SPDCEV`` and replaced with its products.
    Recursion stops when all products are in ``final_state`` or the
    isotope is stable / has no decay data.
    """

    def __init__(self, owner: FlukaDecay,
                 final_state: set[int] | None = None,
                 max_depth: int = 20,
                 on_max_depth: str = "raise"):
        self._owner = owner
        self.final_state = (
            set(STABLE_DEFAULT) if final_state is None else set(final_state)
        )
        self.max_depth = int(max_depth)
        if on_max_depth not in {"raise", "warn"}:
            raise ValueError(
                "on_max_depth must be 'raise' or 'warn', got "
                + repr(on_max_depth)
            )
        self.on_max_depth = on_max_depth

    def expand(self, event):
        """Expand chained decays until all products ∈ final_state."""
        import numpy as np
        import warnings

        owner = self._owner
        # Concatenated product arrays + parent indices.
        pid    = list(event.pid.tolist())
        px     = list(event.px.tolist())
        py     = list(event.py.tolist())
        pz     = list(event.pz.tolist())
        en     = list(event.en.tolist())
        mass   = list(event.mass.tolist())
        parents = [-1] * len(pid)

        # Active queue: indices of products that may need further decay.
        queue: list[tuple[int, int]] = [(i, 0) for i in range(len(pid))]

        while queue:
            i, depth = queue.pop()
            this_pid = pid[i]
            if this_pid in self.final_state:
                continue
            if abs(this_pid) < 1_000_000_000:
                # Standard particle but not in final_state: leave as-is
                # (FLUKA's decay tables don't cover it).
                continue
            A = (abs(this_pid) // 10) % 1000
            Z = (abs(this_pid) // 10_000) % 1000
            m = (abs(this_pid) // 100_000_000) % 10
            iso = owner.lookup(A, Z, m)
            if iso is None or iso.t_half >= 1e38:
                continue
            if depth >= self.max_depth:
                msg = (
                    f"DecayChainHandler hit max_depth={self.max_depth} "
                    f"on ({A},{Z},{m})"
                )
                if self.on_max_depth == "raise":
                    raise RuntimeError(msg)
                warnings.warn(msg)
                continue

            ok, _kdcy, _ilv = owner._lib.chromo_dcy_sample(A, Z, m)
            if not ok:
                continue
            from chromo.models.fluka import FlukaEvent
            sub = FlukaEvent(owner._lib).copy()

            base = len(pid)
            pid.extend(sub.pid.tolist())
            px.extend(sub.px.tolist())
            py.extend(sub.py.tolist())
            pz.extend(sub.pz.tolist())
            en.extend(sub.en.tolist())
            mass.extend(sub.mass.tolist())
            parents.extend([i] * len(sub.pid))
            for j in range(base, base + len(sub.pid)):
                queue.append((j, depth + 1))

            # Mark the parent product as "decayed away" by replacing its
            # pid with 0 (consumers can filter on |pid|>0 or rely on
            # parent indices).
            pid[i] = 0

        # Build a new EventData-shaped record by mutating a copy of the
        # input event in-place, since EventData has many fields.
        chained = event.copy()
        # Drop placeholders (pid==0) before storing.
        keep = [k for k, p in enumerate(pid) if p != 0]
        chained.pid    = np.array([pid[k]    for k in keep], dtype=np.int64)
        chained.px     = np.array([px[k]     for k in keep], dtype=np.float64)
        chained.py     = np.array([py[k]     for k in keep], dtype=np.float64)
        chained.pz     = np.array([pz[k]     for k in keep], dtype=np.float64)
        chained.en     = np.array([en[k]     for k in keep], dtype=np.float64)
        chained.mass   = np.array([mass[k]   for k in keep], dtype=np.float64)
        # Optional: re-map parent indices to the compacted arrays.
        index_map = {old: new for new, old in enumerate(keep)}
        chained.parents = np.array(
            [[index_map.get(parents[k], -1), -1] for k in keep],
            dtype=np.int32,
        )
        return chained
```

Modify `FlukaDecay._fetch_full_catalog` to populate `STABLE_DEFAULT` once:

```python
    def _fetch_full_catalog(self) -> list[FlukaIsotope]:
        import numpy as np

        max_n = 5500
        a = np.zeros(max_n, dtype=np.int32)
        z = np.zeros(max_n, dtype=np.int32)
        m = np.zeros(max_n, dtype=np.int32)
        t12 = np.zeros(max_n, dtype=np.float64)
        exm = np.zeros(max_n, dtype=np.float64)
        jsp = np.zeros(max_n, dtype=np.int32)
        jpt = np.zeros(max_n, dtype=np.int32)
        n = self._lib.chromo_dcy_catalog(max_n, a, z, m,
                                         t12, exm, jsp, jpt)
        cat = [
            FlukaIsotope(
                owner=self,
                A=int(a[i]),
                Z=int(z[i]),
                m=int(m[i]),
                t_half=float(t12[i]),
                mass_excess=float(exm[i]),
                symbol=_z_to_symbol(int(z[i])),
                j_spin=int(jsp[i]),
                j_parity=int(jpt[i]),
            )
            for i in range(int(n))
        ]
        if not STABLE_DEFAULT:
            _populate_stable_default(cat)
        return cat
```

Also: trigger catalog materialisation in `__init__` so `STABLE_DEFAULT`
is ready before any chain expansion. Modify `FlukaDecay.__init__` to
end with:

```python
        self._catalog = None
        # Eagerly materialise so STABLE_DEFAULT is populated for chain ops.
        self._catalog = self._fetch_full_catalog()
```

- [ ] **Step 4: Run tests**

```bash
python -m pytest tests/test_fluka_decay.py::test_stable_default_populated tests/test_fluka_decay.py::test_chain_handler_cs137_includes_gamma -v
```
Expected: 2 PASS.

- [ ] **Step 5: Commit**

```bash
git add src/chromo/models/fluka_decay.py tests/test_fluka_decay.py
git commit -m "feat(fluka): STABLE_DEFAULT + DecayChainHandler.expand"
```

---

## Task 15: `FlukaDecay.chain()` — chain sampling generator

**Files:**
- Modify: `src/chromo/models/fluka_decay.py`
- Modify: `tests/test_fluka_decay.py`

- [ ] **Step 1: Write failing test**

Append to `tests/test_fluka_decay.py`:

```python
def _chain_cs137_method():
    from chromo.models.fluka_decay import FlukaDecay

    dcy = FlukaDecay(seed=11)
    n_events = 0
    n_with_ba137_stable = 0
    for ev in dcy.chain("Cs137", n=200):
        n_events += 1
        # Ba137 g.s. PDG: 1000000000 + 56*10000 + 137*10 + 0
        ba137 = 1000000000 + 56 * 10000 + 137 * 10
        if ba137 in ev.pid:
            n_with_ba137_stable += 1
    return n_events, n_with_ba137_stable


def test_chain_cs137_terminates_on_ba137():
    n, n_ba = run_in_separate_process(_chain_cs137_method)
    assert n == 200
    # All Cs137 chains end on stable Ba137
    assert n_ba >= 195


def _chain_max_depth_raises():
    from chromo.models.fluka_decay import FlukaDecay

    dcy = FlukaDecay()
    try:
        for _ev in dcy.chain("Th232", n=1, max_depth=2):
            pass
    except RuntimeError:
        return "raised"
    return "no_raise"


def test_chain_max_depth_raises():
    assert run_in_separate_process(_chain_max_depth_raises) == "raised"
```

- [ ] **Step 2: Run, verify failure**

```bash
python -m pytest tests/test_fluka_decay.py::test_chain_cs137_terminates_on_ba137 -v
```
Expected: FAIL.

- [ ] **Step 3: Implement `chain()`**

Add inside `class FlukaDecay:`:

```python
    # -- decay chains ---------------------------------------------------

    def chain(self, parent, n: int,
              final_state: set[int] | None = None,
              max_depth: int = 20,
              on_max_depth: str = "raise"):
        """Yield ``n`` decay-chain events for ``parent``.

        Each event contains products from the initial decay plus all
        recursive decays; products are kept (final state) when their
        PDG id is in ``final_state`` (default: ``STABLE_DEFAULT``) or
        when FLUKA cannot decay them.

        Parameters
        ----------
        parent : str | tuple | int
        n : int
        final_state : set[int], optional
            PDG ids to keep as final products. Defaults to
            ``STABLE_DEFAULT``.
        max_depth : int
            Maximum chain depth before erroring/warning.
        on_max_depth : str
            ``"raise"`` (default) or ``"warn"``.
        """
        handler = DecayChainHandler(
            self, final_state=final_state,
            max_depth=max_depth, on_max_depth=on_max_depth,
        )
        for ev in self(parent, n):
            yield handler.expand(ev)
```

- [ ] **Step 4: Run tests**

```bash
python -m pytest tests/test_fluka_decay.py::test_chain_cs137_terminates_on_ba137 tests/test_fluka_decay.py::test_chain_max_depth_raises -v
```
Expected: 2 PASS.

- [ ] **Step 5: Commit**

```bash
git add src/chromo/models/fluka_decay.py tests/test_fluka_decay.py
git commit -m "feat(fluka): FlukaDecay.chain() generator"
```

---

## Task 16: `Fluka.post_event` integration (req 5)

**Files:**
- Modify: `src/chromo/models/fluka.py`
- Modify: `tests/test_fluka_decay.py`

- [ ] **Step 1: Write failing test**

Append to `tests/test_fluka_decay.py`:

```python
def _fluka_post_event():
    from chromo.kinematics import FixedTarget
    from chromo.models import Fluka
    from chromo.models.fluka_decay import (
        FlukaDecay, DecayChainHandler, STABLE_DEFAULT,
    )

    dcy = FlukaDecay(seed=1)
    handler = DecayChainHandler(dcy, final_state=STABLE_DEFAULT,
                                max_depth=20, on_max_depth="warn")
    fluka = Fluka(FixedTarget(100, "p", "O16"),
                  seed=1, post_event=handler.expand)

    n_events = 0
    n_undecayed_residuals = 0
    for ev in fluka(50):
        n_events += 1
        for pdg in ev.pid:
            if abs(int(pdg)) >= 1_000_000_000:
                A = (abs(int(pdg)) // 10) % 1000
                Z = (abs(int(pdg)) // 10_000) % 1000
                m = (abs(int(pdg)) // 100_000_000) % 10
                iso = dcy.lookup(A, Z, m)
                if (iso is not None
                        and iso.t_half < 1e38
                        and int(pdg) not in STABLE_DEFAULT):
                    n_undecayed_residuals += 1
    return n_events, n_undecayed_residuals


def test_fluka_post_event_chain():
    n, undecayed = run_in_separate_process(_fluka_post_event)
    assert n == 50
    assert undecayed == 0, f"{undecayed} unstable residuals leaked through"
```

- [ ] **Step 2: Run, verify failure**

```bash
python -m pytest tests/test_fluka_decay.py::test_fluka_post_event_chain -v
```
Expected: FAIL (`Fluka` doesn't accept `post_event` kwarg).

- [ ] **Step 3: Add `post_event` to `Fluka.__init__`**

Edit `src/chromo/models/fluka.py`:

Find the `__init__` signature:

```python
    def __init__(
        self,
        evt_kin,
        *,
        seed=None,
        targets=None,
        interaction_type=InteractionType.INELASTIC,
```

Add after `interaction_type=InteractionType.INELASTIC,`:

```python
        post_event=None,
```

Inside `__init__`, after the existing init body (before any other state
is set up — after `self._lib` is assigned), add:

```python
        self._post_event = post_event
```

Find the `_generate(self)` method in the same file. After the line that
constructs the `FlukaEvent` (search for `self._event_class(` or the
`yield`), wrap the yield like:

```python
        ev = self._event_class(self._lib).copy()
        if self._post_event is not None:
            ev = self._post_event(ev)
        return ev
```

(Replace whatever returns `FlukaEvent(self._lib).copy()` directly. If
`_generate` returns a single event per call, you may need to read the
existing structure and make the equivalent change — the rule is "after
the event is fully extracted, before it leaves the function, run
`post_event` if set".)

After `chromo_stpxyz` succeeds in `__init__`, mark the FLUKA init flag
so `FlukaDecay` knows decay tables are loaded:

```python
        from chromo.models.fluka_decay import _mark_fluka_init_done

        _mark_fluka_init_done()
```

- [ ] **Step 4: Run tests**

```bash
python -m pytest tests/test_fluka_decay.py::test_fluka_post_event_chain -v
```
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/chromo/models/fluka.py tests/test_fluka_decay.py
git commit -m "feat(fluka): Fluka(post_event=...) for decay-chain post-processing"
```

---

## Task 17: Public re-export

**Files:**
- Modify: `src/chromo/models/__init__.py`
- Modify: `tests/test_fluka_decay.py`

- [ ] **Step 1: Write failing test**

Append to `tests/test_fluka_decay.py`:

```python
def _public_export():
    import chromo.models as M

    return hasattr(M, "FlukaDecay") and hasattr(M, "FlukaIsotope")


def test_public_export():
    assert run_in_separate_process(_public_export) is True
```

- [ ] **Step 2: Run, verify failure**

```bash
python -m pytest tests/test_fluka_decay.py::test_public_export -v
```
Expected: FAIL.

- [ ] **Step 3: Add re-exports**

In `src/chromo/models/__init__.py`, find the section that re-exports
existing models (search for `from .fluka import Fluka`). Add nearby:

```python
from .fluka_decay import (
    DecayChainHandler,
    DecayChannel,
    DecayLine,
    FlukaDecay,
    FlukaIsotope,
    STABLE_DEFAULT,
)
```

Update the module's `__all__` list (if present) to include these names.

- [ ] **Step 4: Run tests, full suite**

```bash
python -m pytest tests/test_fluka_decay.py -v
```
Expected: ALL PASS (≥ 18 tests).

- [ ] **Step 5: Run pre-commit on the changed files**

```bash
pre-commit run --files src/chromo/models/fluka_decay.py src/chromo/models/fluka.py src/chromo/models/__init__.py tests/test_fluka_decay.py 2>&1 | tail -10
```
Expected: All hooks pass (or auto-fix and re-run).

- [ ] **Step 6: Commit**

```bash
git add src/chromo/models/__init__.py tests/test_fluka_decay.py
git commit -m "feat(fluka): re-export FlukaDecay from chromo.models"
```

---

## Task 18: Acceptance suite — Th-232 chain + Co-60 spectrum

**Files:**
- Modify: `tests/test_fluka_decay.py`

- [ ] **Step 1: Write final integration tests**

Append to `tests/test_fluka_decay.py`:

```python
def _th232_chain_terminates():
    from chromo.models.fluka_decay import FlukaDecay

    dcy = FlukaDecay(seed=2)
    events = list(dcy.chain("Th232", n=20, max_depth=20))
    # Final products should be stable Pb208 (and leptons / gammas).
    pb208 = 1000000000 + 82 * 10000 + 208 * 10
    return [pb208 in ev.pid for ev in events]


def test_th232_chain_terminates_on_pb208():
    flags = run_in_separate_process(_th232_chain_terminates)
    # Th232 series ultimately ends on Pb-208.
    assert sum(flags) >= 18, sum(flags)


def _co60_two_gammas():
    from chromo.models.fluka_decay import FlukaDecay

    dcy = FlukaDecay(seed=3)
    n_two_g = 0
    for ev in dcy("Co60", n=500):
        n_g = sum(1 for p in ev.pid if int(p) == 22)
        if n_g >= 2:
            n_two_g += 1
    return n_two_g


def test_co60_emits_two_gammas():
    n = run_in_separate_process(_co60_two_gammas)
    # Co-60 → Ni-60* → Ni-60 + 1.17 + 1.33 MeV in ~all decays.
    assert n >= 450, n
```

- [ ] **Step 2: Run final suite**

```bash
python -m pytest tests/test_fluka_decay.py -v
```
Expected: ALL PASS.

- [ ] **Step 3: Run a broader subset to catch regressions**

```bash
python -m pytest tests/test_fluka_decay.py tests/test_decay_handler.py -v -n 4 2>&1 | tail -20
```
Expected: no new failures vs main.

- [ ] **Step 4: Commit**

```bash
git add tests/test_fluka_decay.py
git commit -m "test(fluka): Th232 chain termination + Co60 two-gamma spectrum"
```

---

## Self-Review

**Spec coverage:**
- §1 Goal 1 (catalog, T1/2 < threshold): Tasks 3 (Fortran wrapper) + 10 (Python `catalog(t_half_min/max=)`).
- §1 Goal 2 (inspection): Tasks 4, 5 (wrappers) + 8, 11, 12 (`short()`, `channels`, line-list properties + `__str__`).
- §1 Goal 3 (inclusive sampling): Tasks 6 (`SPDCEV` wrapper) + 13 (`__call__`).
- §1 Goal 4 (chains): Tasks 14 (`STABLE_DEFAULT` + `DecayChainHandler.expand`) + 15 (`chain()` method).
- §1 Goal 5 (Fluka post_event chain): Task 16.
- §6 single-instantiation guard: Task 9 (init guard) + Task 16 (`_mark_fluka_init_done` from Fluka).
- §7 data flow walkthrough: covered by Task 14's chain test on Cs-137 → Ba-137.
- §8 error handling: invalid lookup (Task 9), stable-isotope sample raises (Task 13), max_depth raises (Task 15), 2 instances of FlukaDecay rejected (Task 9).
- §9 testing matrix: every named test (catalog filters, lookup, channel mode, sample, chain, Th232 termination, Co-60) has a matching task.

**Placeholder scan:** No "TBD", "TODO", "implement later" tokens. Each step has explicit code or commands.

**Type consistency:**
- `FlukaIsotope.__init__` signature matches usage in `_fetch_full_catalog` (Task 10) and `lookup` (Task 9). Both pass `owner=` kwarg-only and the named fields.
- `DecayChainHandler.expand(event)` matches the `post_event=callable(event) -> event` contract used in Task 16 (`Fluka(post_event=handler.expand)`).
- `_fetch_channels` / `_fetch_lines` signatures (Tasks 11, 12) match what `FlukaIsotope.channels` / `gamma_lines` etc. call.
- Fortran wrapper signatures match the Python ndarray-passing pattern (Tasks 1–6 are consistent in argument order: `(in_args..., out_arrays..., out_count)`).
- The PDG ion code expression `1_000_000_000 + Z * 10_000 + A * 10 + m` is used identically in Tasks 14 and 16.

---

## Execution Handoff

Plan complete and saved to `docs/superpowers/plans/2026-04-30-fluka-decay-interface.md`. Two execution options:

**1. Subagent-Driven (recommended)** — fresh subagent per task, two-stage review between tasks, fast iteration.

**2. Inline Execution** — execute tasks in this session using executing-plans, batch execution with checkpoints.

Which approach?
