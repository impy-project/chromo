# FLUKA radioactive-decay interface for chromo — design

Status: draft for review
Date: 2026-04-30
Author: A. Fedynitch (with Claude)

## 1. Goals

A new `FlukaDecay` Python class that exposes FLUKA 2025.1's radioactive-decay
machinery to chromo users, covering five requirements:

1. **Catalog query.** Walk FLUKA's nuclear-decay database and return all
   isotopes/isomers, with filtering by half-life, A, Z, decay mode.
2. **On-demand inspection.** Pretty-printed table of an isotope's
   properties (mass, decay channels, gamma/alpha/beta-spectrum line lists,
   total Q) — equivalent to `DCYPRN`'s console dump but driven from
   Python.
3. **Inclusive decay sampling.** Given an isotope, generate N correlated
   decay events as `EventData`-shaped records. FLUKA's branching ratios
   weight the channel mixture; users do not pick the channel.
4. **Decay chains.** A user supplies a "final-state" set of allowed
   particles. `FlukaDecay` decays the parent; for each product that is
   an unstable nuclide *not* in the final-state set, it samples a
   further decay; recursion continues until every product is in the
   set. The returned event contains the union of all products with
   correct parent links.
5. **Optional: chain post-processing on `Fluka` hadronic events.** The
   same chain machinery is applied after each `Fluka(...)` event so
   unstable hadronic residuals are decayed until all final products are
   in the user's allowed set. Single-process coexistence between
   `Fluka` and `FlukaDecay` is the prerequisite — see §6 / §11.

## 2. Out of scope

- Per-channel sampling. `SPDCEV` already weights channels by FLUKA's
  internal branching ratios; exposing a "sample only β−" mode is not
  needed.
- Modifying `Fluka`'s existing public methods. Req 5 is bolted on via
  an opt-in keyword (a callable post-processor), not by changing the
  default behaviour.
- Two simultaneous `FlukaDecay` instances. Single-instantiation per
  process applies (§6).
- The full FRDM mass table (~13 k entries). Only the ~4 500 nuclides
  with valid `ISMRCH` decay data are exposed. A future `mass_only=True`
  catalogue mode is left as a hook.

## 3. Architecture

```
src/chromo/models/
  fluka.py             (existing — add optional `post_decay` kwarg)
  fluka_decay.py       (NEW — FlukaDecay, FlukaIsotope, DecayChannel,
                        DecayLine, DecayChainHandler, _fluka_decay_init)
src/fortran/fluka/
  chromo_fluka.f       (+ 6 new wrapper subroutines)
src/chromo/models/__init__.py   (re-export FlukaDecay)
```

`FlukaDecay` is **not** an `MCRun` subclass — there is no projectile/target
kinematics, no cross-section concept. It is a sibling class with its own
lifecycle but sharing the FLUKA library's process globals via a
module-level init guard (§6).

`FlukaEvent` (already defined in `fluka.py`) is reused as the event type
returned by `FlukaDecay(...)` and `FlukaDecay.chain(...)`. This keeps
downstream tooling — HepMC export, ROOT writer, particle filtering —
working without changes.

## 4. Fortran wrappers

Six new subroutines added to `src/fortran/fluka/chromo_fluka.f`. All use
`INCLUDE '(DBLPRW)'` (the prototype proved that `DBLPRC` corrupts FHEAVY
across `SPDCEV` calls; `DBLPRW` pulls in the blank common where SPDCEV
writes line tables).

```fortran
SUBROUTINE CHROMO_DCY_INIT()
* Idempotent. Mirrors dcytst.f's main-program init:
*   CMSPPR / ZEROIN / flag setup / NCDTRD / KPIXSR / INCINI /
*   minimal EMF flag setup / RDFLUO
* Sets a module-level LDCYINI flag. No-op if STPXYZ has already run
* (chromo_stpxyz internally calls NCDTRD/RDFLUO via FLUKA's master init).

SUBROUTINE CHROMO_DCY_LOOKUP(IA, IZ, IM, IFOUND, T12, EXM, JSP, JPT)
* Thin ISMRCH wrapper. IFOUND = 1 if found (KISITP > 0), else 0.
* The element symbol is derived in Python from Z via the existing
* `particle` package (or `chromo.util`); no Fortran-side string return
* is needed (avoids Fortran string interop).

SUBROUTINE CHROMO_DCY_CATALOG(MAX_N, N_OUT,
&                             A_OUT, Z_OUT, M_OUT,
&                             T12_OUT, EXM_OUT,
&                             JSP_OUT, JPT_OUT)
* One-shot full walk over the conservative valley-of-stability Z-band
* per A (proven safe in the prototype, yields 4 478 entries).
* Probes ground state first; isomers (m=1..4) only probed if g.s. exists,
* avoiding ISMRCH bounds violations.
* MAX_N = caller-allocated output buffer size (allocate ≥ 5 000).

SUBROUTINE CHROMO_DCY_CHANNELS(IA, IZ, IM, MAX_CH, N_CH,
&                              KIND_CH, BR_CH,
&                              DA_CH, DZ_CH, DM_CH, Q_CH)
* Per-isotope channel data. Reads BRDECY/BRDISM, IDCNUC/IDCISM,
* IDCYDA/IDCYDZ, calls QRDDCY for the Q value.
* KIND_CH(:) integer codes (1=alpha, 2=B-, 3=B+, 4=EC, 5=IT, 6=SF,
* 7=B-N, 8=B+P, ...) — caller maps via a Python lookup table.
* DA/DZ/DM = -1 if no single daughter (e.g. SF).

SUBROUTINE CHROMO_DCY_LINES(IA, IZ, IM, KIND, MAX_L, N_L,
&                           BR_L, E_L, NLEV_L, LPOS_L)
* Per-isotope line list. KIND = 1 (γ), 2 (α), 3 (CE/Auger), 4 (β±).
* For β±: E_L is end-point, NLEV_L is daughter level index. The
* positron flag is captured separately in LPOS_L (1 = e+, 0 = e-);
* FLUKA's SIGGTT layout encodes it in the sign of <E> at offset +2,
* which we read before overwriting E_L with the unsigned endpoint.
* For γ/α/CE: LPOS_L is always 0.
* Reads via NGMLNS/KGMLNS, NCELNS/KCELNS, NALLNS/KALLNS, NBTSPC/KBTSPC
* and walks SIGGTT(:) (a REAL*4 module pointer, not a function).

SUBROUTINE CHROMO_DCY_SAMPLE(IA, IZ, IM, LSUCCS, KDCY_OUT, ILV_OUT)
* Zeroes NP/NP0/NPHEAV/NPEMF, calls SPDCEV with the stock argument
* tuple from dcytst.f, then FLLHEP (which itself resets NHEP=0 on
* entry, so no manual HEPEVT scrub is needed in the Python sampling
* loop). Products land in HEPEVT (f2py attribute `_fluka.hepevt`,
* not `_fluka.hepcmm` -- the include file is (HEPCMM) but the COMMON
* block name is /HEPEVT/) for chromo's existing extraction pipeline.
* KDCY_OUT: 1=alpha, 2=B-, 3=B+, 4=EC, 5=IT, 6=SF, -1=other.
* Note: KDCY_OUT is restricted to 1..6 even though chromo_dcy_channels
* returns 1..12 -- DCYFLG only has six fundamental-mode flags.
* ILV_OUT: daughter level (ILVDCY).
```

Six new entries appended to `meson.build`'s `fluka_syms` list:
`'chromo_dcy_init', 'chromo_dcy_lookup', 'chromo_dcy_catalog',
'chromo_dcy_channels', 'chromo_dcy_lines', 'chromo_dcy_sample'`.

## 5. Public Python API

### 5.1 Construction

```python
from chromo.models.fluka_decay import FlukaDecay

dcy = FlukaDecay(seed=42)
```

`seed` is forwarded to `init_rng_state` (existing chromo helper) so
sampled decays are reproducible. If a `Fluka` instance already exists in
the process, init is a no-op (it already loaded NCDTRD/RDFLUO via
STPXYZ); otherwise `chromo_dcy_init` runs the lighter dcytst.f-style
init. Either ordering works.

### 5.2 Catalogue (req 1)

```python
catalog = dcy.catalog()              # list[FlukaIsotope], ~4 500 entries
short = dcy.catalog(t_half_max=1.0)  # T1/2 < 1 s
long  = dcy.catalog(t_half_min=1e10) # T1/2 > 10^10 s
```

Filter kwargs (all optional, applied in Fortran-loaded full-list
materialisation step in Python):

| kwarg          | meaning                                     |
|----------------|---------------------------------------------|
| `t_half_min`   | keep only T1/2 ≥ this (seconds)             |
| `t_half_max`   | keep only T1/2 ≤ this (seconds)             |
| `a_min`/`a_max`| mass-number range                           |
| `z_min`/`z_max`| atomic-number range                         |
| `decay_modes`  | iterable of mode strings: `{"B-", "EC", …}` to keep only isotopes with any matching channel; triggers per-isotope `channels` fetch for filter evaluation |

The full catalog is fetched once on first `.catalog()` call and cached
on the `FlukaDecay` instance. Subsequent calls are filtered Python views.

### 5.3 Lookup and inspection (req 2)

```python
iso = dcy.lookup("Cs137")             # str: "Cs137", "Tc99m", "U238m1"
iso = dcy.lookup(137, 55, 0)          # tuple form (A, Z, m)
iso = dcy.lookup(pdg_id)              # PDG ion code
print(iso)                            # full DCYPRN-style table
print(iso.short())                    # one-line summary
```

Naming follows `chromo.util.process_particle` conventions. Trailing `m`
or `m1` denotes 1st isomer; `m2`, `m3` follow.

`FlukaIsotope` (regular class with `__slots__`, lazy attributes for the
heavy data):

```python
class FlukaIsotope:
    __slots__ = ("A", "Z", "m", "t_half", "mass_excess",
                 "symbol", "j_spin", "j_parity",
                 "_channels", "_lines", "_owner")
    # _lines: dict[str, tuple[DecayLine, ...]] keyed by kind
    #   ("gamma" | "alpha" | "ce" | "beta"); populated lazily per kind.
    # _channels: tuple[DecayChannel, ...] | None; populated on first access.
    # _owner: weakref-like back-pointer to the FlukaDecay instance whose
    #   Fortran library will service the lazy fetches.

    @property
    def channels(self) -> tuple[DecayChannel, ...]:
        # Lazy: calls chromo_dcy_channels on first access, caches.

    @property
    def gamma_lines(self) -> tuple[DecayLine, ...]: ...
    @property
    def alpha_lines(self) -> tuple[DecayLine, ...]: ...
    @property
    def ce_lines(self) -> tuple[DecayLine, ...]: ...
    @property
    def beta_spectra(self) -> tuple[DecayLine, ...]: ...

    def __str__(self) -> str:
        # DCYPRN-style multi-line table; lazy fetches all line lists.

    def short(self) -> str:
        # "Cs137 (Z=55, m=0): T1/2=9.49e+08 s, ExM=-86.55 MeV"

@dataclass(frozen=True, slots=True)
class DecayChannel:
    name: str                # "B-", "EC", "alpha", "IT", "SF", "B-N", ...
    branching: float         # 0..1
    daughter_A: int          # -1 if no single daughter (SF)
    daughter_Z: int
    daughter_m: int
    q_value: float           # MeV

@dataclass(frozen=True, slots=True)
class DecayLine:
    energy: float            # MeV
    branching: float         # 0..1 of the parent decay
    end_level: int           # for α / β±, daughter level index; 0 otherwise
    is_positron: bool        # for β only
```

The `__str__` implementation reproduces `DCYPRN`'s output: header line
with A/sym/Z/m/T1/2/ExMass; the channel table; CE / α / γ / β line
lists; total decay energy with and without neutrinos. Lazy lookups
ensure that `print(iso)` triggers the channel + line-list fetches only
on demand.

### 5.4 Inclusive decay sampling (req 3)

```python
for ev in dcy("Cs137", n=50_000):
    pids = ev.final_state().pid
```

`dcy.__call__(parent, n)` returns a generator yielding `FlukaEvent`
instances. Each call to `chromo_dcy_sample` populates HEPEVT through
`FLLHEP`; chromo's existing `FlukaEvent` extraction pipeline returns the
event. No initial-beam record is emitted (matches
`FlukaEvent._prepend_initial_beam = pass`).

### 5.5 Decay chains (req 4)

```python
from chromo.models.fluka_decay import STABLE_DEFAULT

final_state = STABLE_DEFAULT | {"Cs137"}   # extend or replace as needed
# STABLE_DEFAULT contains: e±, γ, ν_{e,μ,τ} and antiparticles, p, n,
# plus all FLUKA isotopes with T1/2 == 1e38 (computed once on first
# FlukaDecay instantiation from the cached catalog).

for ev in dcy.chain("Cs137", n=10_000,
                    final_state=final_state,
                    max_depth=20):
    # ev contains the union of products across all chained decays.
    # ev.parents links each decay product to the index of its parent
    # in the same event.
    ...
```

Semantics, exactly as the user specified:

> If a decay product is an unstable nuclide that is **not in
> `final_state`**, sample its decay and replace it with its products.
> Recurse until every product is either in `final_state` or cannot be
> decayed further by FLUKA (T1/2 ≥ 1e38 or no decay data).

Defaults:
- If `final_state` is omitted, the implementation uses `STABLE_DEFAULT`
  = `{e±, γ, ν, ν̄, p, n, all isotopes with T1/2 == 1e38}` so the chain
  runs to the lowest bound nucleus.
- `max_depth=20` (uranium series is 14, plenty of headroom). Exceeding
  it raises `RuntimeError`; configurable with `on_max_depth="warn"` to
  truncate instead.

A helper `STABLE_DEFAULT` is exported from `fluka_decay` so users can
extend it (e.g. add `"Cs137"` to halt at Cs137 even though it would
otherwise decay).

`DecayChainHandler` is the underlying class:

```python
class DecayChainHandler:
    def __init__(self, owner: FlukaDecay, final_state: set,
                 max_depth: int = 20, on_max_depth: str = "raise"): ...

    def expand(self, event: EventData) -> EventData:
        # Walks the event's product list; for each unstable nuclide
        # not in final_state, calls owner._sample_decay(A, Z, m), merges
        # the products into the event with correct parent indices.
```

`expand()` is what `dcy.chain(...)` calls per event. It's also the entry
point for req 5 (§5.6).

### 5.6 Optional: chain mode for Fluka hadronic events (req 5)

```python
from chromo.models import Fluka
from chromo.models.fluka_decay import FlukaDecay, DecayChainHandler

dcy   = FlukaDecay(seed=42)
chain = DecayChainHandler(dcy, final_state=final_state, max_depth=20)
fluka = Fluka(FixedTarget(100, "p", "O16"),
              seed=42,
              post_event=chain.expand)         # NEW kwarg

for ev in fluka(1_000):
    # ev is the hadronic event, with all unstable residuals decayed
    # via SPDCEV until products are in final_state.
    ...
```

Implementation: `Fluka.__init__` accepts an optional `post_event`
callable. After each `_generate()`, before the event is yielded, the
callable is applied: `event = post_event(event)`. If `post_event is
None`, behaviour is unchanged.

`DecayChainHandler.expand` re-enters Fortran via `chromo_dcy_sample` for
each unstable nuclear residual. This works because:

1. Both `chromo_evtxyz` (hadronic) and `chromo_dcy_sample` are linked
   into the same `_fluka` extension module — no separate library load.
2. `chromo_evtxyz` populates HEPEVT; `FlukaEvent.__init__` extracts the
   numpy arrays before returning, so subsequent HEPEVT writes by
   `chromo_dcy_sample` cannot corrupt the already-extracted arrays.
3. `STPXYZ` (run by `Fluka.__init__`) calls `NCDTRD`/`RDFLUO` as part
   of FLUKA's master init; the decay tables are already loaded when
   `_dcy_sample` runs.
4. `chromo_dcy_sample` zeroes its own NP/NP0/NPHEAV/NPEMF before each
   `SPDCEV`, so it never inherits dirty state from `chromo_evtxyz`.

Risks (to verify in the implementation phase):
- `SPDCEV`'s blank-common access (DBLPRW) must coexist with the rest of
  the FLUKA library's blank-common usage during a hadronic run. The
  chromo `_fluka` extension is built with the same linker flags as
  `dcytst`, so this is expected to work — but a smoke test that runs
  `Fluka` then `FlukaDecay.expand` on a real residual is a required
  acceptance test (§9).
- HEPEVT clobbering vs intermediate cross-section queries: any
  in-flight call (e.g. `cross_section()`) that depends on HEPEVT must
  not happen between event yield and chain expansion. Solved by
  running `expand` synchronously inside `_generate()` before yielding.

If the smoke test fails (one of the unforeseen interactions), req 5
falls back to a "post-process EventData in pure Python via repeated
FlukaDecay calls in a separate process", which is strictly less
efficient but guaranteed safe. Stated as a fallback in §11.

## 6. Single-instantiation, init guard

FLUKA's library uses Fortran globals; one process can only have its
master init done once. The constraint applies across `Fluka` and
`FlukaDecay` together.

Module-level guard in `fluka_decay.py`:

```python
_INITIALIZED = False    # set True after CHROMO_DCY_INIT or chromo_stpxyz

def _ensure_init():
    global _INITIALIZED
    if _INITIALIZED:
        return
    _lib.chromo_dcy_init()
    _INITIALIZED = True
```

`Fluka.__init__` sets `_INITIALIZED = True` after `chromo_stpxyz`
(STPXYZ subsumes `chromo_dcy_init`'s work). `FlukaDecay.__init__` calls
`_ensure_init()`.

Either ordering works. Two `FlukaDecay` instances in the same process
is rejected with a clear error (mirrors `Fluka`'s existing
single-instance guard).

## 7. Data flow — chain example (Cs-137 → final state)

User-visible sequence for `dcy.chain("Cs137", n=1, final_state={e-, e+, γ, ν, ν̄, stable_nuclides})`:

1. `_sample_decay(137, 55, 0)` → SPDCEV. Products in HEPEVT:
   `e- (β-), ν̄_e, Ba137_m1 (excited)`.
2. `expand()` scans: `e-` and `ν̄_e` are in final_state — keep.
   `Ba137_m1` is unstable, not in final_state — recurse.
3. `_sample_decay(137, 56, 1)` → SPDCEV. Products: `γ (0.66 MeV),
   Ba137 (g.s.)`.
4. `expand()` scans the new products: `γ` in final_state — keep.
   `Ba137` is stable → in `STABLE_DEFAULT` → keep.
5. Returned event particle list: `[e-, ν̄_e, γ, Ba137]`. Parents:
   `[parent=Cs137, parent=Cs137, parent=Ba137_m1, parent=Ba137_m1]`.

Implementation invariant: each `_sample_decay` call clobbers HEPEVT, so
the intermediate event arrays are extracted into Python *before* the
next call. The merge happens entirely in Python on `EventData` arrays.

## 8. Error handling

| Condition                                  | Behaviour                                |
|--------------------------------------------|------------------------------------------|
| `lookup("XXX")` invalid name               | `ValueError` with parsing detail         |
| `lookup(A, Z, m)` not in catalog           | returns `None`                           |
| `dcy("Fe56", N)` (stable nuclide)          | `ValueError("Fe56 has no decay data")`   |
| `chromo_dcy_sample` returns LSUCCS=False   | retry once; if still failed, skip event with `chromo.info` warning |
| `chain(..., max_depth=20)` exceeded        | `RuntimeError`; `on_max_depth="warn"` truncates instead |
| Two `FlukaDecay()` calls in same process   | `RuntimeError("FLUKA singleton...")`     |
| `decay_modes` filter invalid               | `ValueError` listing valid modes         |

## 9. Testing

`tests/test_fluka_decay.py`, all wrapped via
`tests/util.py::run_in_separate_process` to satisfy the FLUKA
single-instantiation constraint.

Catalogue:
- `test_catalog_size`: full catalog has between 4 000 and 5 000 entries.
- `test_catalog_known_isotopes`: U-238, Cs-137, C-14 each appear with
  T1/2 within 1 % of accepted values.
- `test_catalog_filters`: `t_half_min=1e10` returns ≥ 200; `a_min=200,
  a_max=240` shrinks proportionally; `decay_modes={"alpha"}` returns
  the alpha-emitter subset.

Lookup / inspection:
- `test_lookup_isomer`: `lookup("Tc99m")` returns m=1 with T1/2 ≈ 6 hr.
- `test_lookup_invalid`: `lookup("Foo")` → ValueError.
- `test_str_format`: `str(iso)` for Cs-137 contains "Beta-", "94.7",
  the 0.661 MeV gamma line, and the total decay energy.

Sampling:
- `test_decay_cs137`: 10 000 events; ≥99 % β−; ⟨γ multiplicity⟩ ≈ 0.94;
  ⟨E_β−⟩ ≈ 0.18 MeV.
- `test_decay_co60`: 10 000 events; ⟨γ multiplicity⟩ ≈ 2.0 with the
  1.17 and 1.33 MeV lines dominant.
- `test_decay_stable_raises`: `dcy("Fe56", 1)` → ValueError.
- `test_seed_reproducibility`: same seed yields identical first-event
  product list.

Chains:
- `test_chain_cs137`: every event ends with `Ba137 (g.s.) + e- + ν̄_e
  + 0 or 1 γ`. Parents recorded.
- `test_chain_th232`: 14-step uranium-thorium series terminates on
  Pb-208; depth ≤ 14.
- `test_chain_max_depth_raises`: artificial `max_depth=2` on Th232 →
  RuntimeError (default) or warn (`on_max_depth="warn"`).
- `test_chain_final_state_halts`: `final_state` includes `Ba137_m1` →
  chain stops there, only β− and ν̄_e are emitted.

Coexistence (req 5):
- `test_fluka_then_decay`: build `Fluka(FixedTarget(100,"p","O16"))`,
  then `FlukaDecay()`, then `Fluka` events still produce sensible
  cross-sections.
- `test_fluka_post_event_chain`: `Fluka(post_event=chain.expand)` runs
  100 events; final-state arrays contain only particles in
  `final_state`. Assertion: no isotope in the output has
  `t_half < 1e38` *and* is outside `final_state`.

## 10. Performance

- `catalog()` first call: ~ 4 500 ISMRCH calls (≈ 0.5 s end-to-end)
  plus filter walk in Python. Result cached on the instance.
- `iso.channels` first access: 1 Fortran call (~ µs).
- Sampling: same throughput as `dcytst.f` (≥ 10 k events/s for Cs-137).
- Chain expansion: dominated by recursion depth; ~ 14 recursive calls
  per event for Th232 series.

## 11. Future / risks

- **Req 5 fallback.** If `FlukaDecay.expand` cannot be safely invoked
  inside the same process as a running `Fluka`, the documented
  workaround is: pickle the EventData, send it to a worker subprocess
  that hosts `FlukaDecay`, expand there, return. Slower but bypasses
  any unforeseen common-block interaction.
- **Mass table / FRDM (~13 k entries).** Add `mass_only=True` to
  `catalog()` later if needed; entries would lack `channels` /
  line-lists. Pure read of `WAPS` and `Inwadm` tables.
- **Per-channel sampling.** If a future user wants "sample only β−",
  expose via a `channel="B-"` kwarg on `__call__` that pre-selects
  the channel inside the SPDCEV call (small Fortran change).
- **Decay-time Monte Carlo.** Currently only the products are returned;
  an optional `with_time=True` mode could also return sampled decay
  times from the exponential distribution per parent (purely Python).

## 12. Open questions

None blocking. The single design decision deferred to implementation is
whether `decay_modes` filter materialises `channels` for every isotope
(simple but slow) or maintains a Fortran-side bitmap (fast but more
work). Default to the simple Python-side path; revisit if perf is bad.
