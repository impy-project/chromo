# Robust FLUKA integration in chromo — design

Date: 2026-04-15
Status: design approved, pending user review of written spec
Branch: feature_fluka

## Goal

Provide a production-grade FLUKA model in chromo that supports the full physics
scope the user cares about: hadron–nucleon, hadron–nucleus, and nucleus–nucleus
collisions; photohadronic and photonuclear interactions; electromagnetic
dissociation (EMD); and extraction of nuclear remnants. The implementation must
link to the official FLUKA 2025.1 distribution and re-use FLUKA's own routines
for conversions and event generation wherever possible. The current draft on
`feature_fluka` is "kind of working" for a few hard-coded targets and has not
been validated against the broader physics surface; this design covers the
gap.

## Non-goals

- Weizsäcker-Williams equivalent-photon machinery for e± projectiles.
  `e+/e-` are not supported as projectiles; document and refuse at
  `_check_kinematics`.
- Dynamic material re-registration via subprocess restart. The user-supplied
  `targets=` kwarg at construction time is the single extension point.
- Custom output format for nuclear remnants. They appear as entries in the
  standard HEPEVT/`EventData` arrays with PDG nucleus ids
  (`10LZZZAAAI` encoding).
- Runtime mutation of FLUKA's decay settings. `_set_stable` remains a no-op
  that emits an info-level log.
- Bundling the FLUKA binary distribution with chromo. Users install FLUKA
  separately per its license; chromo build fails fast if `$FLUPRO` is unset.

## User contract

```python
from chromo.models import Fluka
from chromo.kinematics import FixedTarget, CenterOfMass
from chromo.util import CompositeTarget

# Default: hadron/photon/nucleus projectile on a registered nucleus.
gen = Fluka(FixedTarget(100, "p", "O16"), seed=42)
for event in gen(10):
    print(event.final_state_charged().pid)

# Extend the registered target list at construction.
gen = Fluka(FixedTarget(100, "p", "Cu63"), targets=["Cu63", "Au197"], seed=42)

# CompositeTarget (air) — components must resolve to registered nuclei.
air = CompositeTarget([("N14", 0.78), ("O16", 0.21), ("Ar40", 0.01)])
gen = Fluka(FixedTarget(100, "p", air), seed=42)

# AA collision.
gen = Fluka(FixedTarget(100 * 16, "O16", "Pb208"), seed=42)

# Photohadronic / photonuclear.
gen = Fluka(FixedTarget(0.3, "gamma", "Pb208"))   # Δ-resonance region

# EMD and interaction-mode selection.
from chromo.models.fluka import InteractionType
gen = Fluka(FixedTarget(100, "p", "Pb208"),
            interaction_type=InteractionType.INELA_ELA_EMD)
xs = gen.cross_section()
print(xs.inelastic, xs.elastic, xs.emd)
```

Frame handling: the generator runs internally in fixed-target. Users supplying
`CenterOfMass` kinematics get correct results via the standard chromo
base-class conversion path. Events are returned in the user's requested frame.

## Architecture

Five well-bounded units:

1. **FLUKA installation + meson build glue** — `scripts/install_fluka.sh`,
   `meson.build` `fluka` block.
2. **Fortran wrapper layer** — `src/fortran/fluka/chromo_fluka.f` (single
   file, extends the existing draft; ~500 LOC).
3. **Python wrapper class** — `src/chromo/models/fluka.py` (extends draft).
4. **Event class** — `FlukaEvent(MCEvent)` in the same file.
5. **Test suite** — `tests/test_fluka.py`.

### Existing code reused from the `feature_fluka` draft

Keep:
- All of `chromo_fluka.f` except the typo fix (see below). `init_rng_state`,
  `load_rng_state`, `save_rng_state`, `fluka_rand`, `icode_from_pdg`,
  `icode_from_pdg_arr`, `charge_from_pdg_arr`, `pdg_from_icode`,
  `random_direction`, `ICRVRCK`, `fluka_particle_scheme`, `CHROMO_STPXYZ`,
  `CHROMO_EVTXYZ`, `CHROMO_FLLHEP` stay as-is.
- The `fluka_all.a` custom-target in `meson.build` that combines the five
  FLUKA archives (`libflukahp.a`, `libdpmmvax.a`, `librqmdmvax.a`,
  `latestRQMD/librqmd.a`, `interface/libdpmjet*.a`).
- The `InteractionType` enum skeleton in `fluka.py` (extended to `IntEnum`,
  values unchanged).
- `FlukaEvent`'s `_get_charge` via `charge_from_pdg_arr`, plus the empty hooks
  `_history_zero_indexing`, `_prepend_initial_beam`, `_repair_initial_beam`.
- `_init_rng`, `_cleanup_fort`, the FLUKA `FLUPRO` env-var check in the Python
  constructor.
- The full-material/`_DEFAULT_MATERIALS` skeleton in `_init_fluka_materials`;
  generalise per Section "Material handling" below.

### Gaps and bugs to close

- **Bug.** `CHROMO_SGMXYZ` at `src/fortran/fluka/chromo_fluka.f:232` assigns to
  `CHROMO_SGMXY` (typo) — returns undefined value. Fix to
  `CHROMO_SGMXYZ = SGMXYZ(...)`.
- **Bug.** `meson.build` reads `meson.project_source_root() +
  '/../FLUKA/interface/dpmvers'` which hard-codes a path outside the repo.
  Replace with `$FLUPRO/interface/dpmvers` (`fluprod + '/interface/dpmvers'`).
- **Gap.** No ion-projectile path. `_fluka_pid` returns 0 for nuclei (because
  `MCIHAD` only covers -6..390). Must add PDG→FLUKA-ion-code conversion and
  route to EVTXYZ/SGMXYZ via FLUKA's own `PDGION`/`SETION`.
- **Gap.** No `CompositeTarget` per-component switching in `_set_kinematics`.
- **Gap.** No tests at all.
- **Gap.** No validation of energy bounds per projectile/target class. The
  draft's flat `_ekin_max = 20 TeV` is wrong for photonuclear and for AA.
- **Stale comment.** The reference `evtxyz.f` has `CALL FLABRT('EVTXYZ','EMD
  NOT YET IMPLEMENTED')`. FLUKA 2025.1 ships an updated EVTXYZ that supports
  EMD. chromo must link the library version, not include the reference `.f`
  file, and must not skip EMD.

## Python class structure

`src/chromo/models/fluka.py`, single `Fluka(MCRun)` class plus
`FlukaEvent(MCEvent)`. Class attributes:

```python
class Fluka(MCRun):
    _name = "FLUKA"
    _event_class = FlukaEvent
    _frame = EventFrame.FIXED_TARGET
    _version = "2025.1"
    _library_name = "_fluka"
    _projectiles = standard_projectiles | Nuclei() | {lp.photon.pdgid}
    _targets = Nuclei()
    _ecm_min = 0.01 * GeV                 # GDR region for γA
    _ekin_per_nucleon_max_hadron = 20 * TeV
    _ekin_per_nucleon_max_photon = 1 * TeV   # refined by investigation 6
```

Constructor:

```python
def __init__(self, evt_kin, *, seed=None,
             targets=None,                    # list of extra nuclei to register
             interaction_type=InteractionType.INELASTIC,
             transition_energy=None,          # None → FLUKA default
             transition_smearing=None,
             enable_quasielastic=False,
             rng_state_file=None):
    ...
```

Internal helpers:

- `_build_material_list(evt_kin, user_targets)` — returns an ordered,
  deduplicated list of nucleus PDG ids. Combines
  `_DEFAULT_MATERIALS` ∪ `user_targets` ∪ `components(evt_kin.p2)`. Free
  proton is represented as `2212` (same Z=1 as `H1`). Internally normalises
  all entries to canonical PDG nucleus ids via `particle.Particle.from_pdgid`.
- `_init_fluka_once(materials)` — calls `chromo_stpxyz` exactly once. Stores
  `{pdg → fluka_material_idx}` mapping as `self._materials_map`. Raises if
  called twice (STPXYZ would `FLABRT` anyway).
- `_fluka_projectile_code(pdg)` — delegates to the Fortran helper
  `pdg_to_proj_code` (see Fortran layer).
- `_set_kinematics(kin)` — resolves the target's FLUKA material index from
  `self._materials_map`. Drives from `_temporary_kinematics` so
  `CompositeTarget` per-component iteration works via the base class.
- `_cross_section(kin)` — calls `chromo_sgmxyz` once per channel requested by
  `self._interaction_type` (up to three calls), packs into
  `CrossSectionData(inelastic, elastic, emd, ...)`.
- `_generate()` — calls `chromo_evtxyz` with projectile code, target material
  index, `kin.ekin`, `ppm=0`, direction `(0, 0, 1)`, and the interaction-type
  flag. Returns True unless the flag indicates a FLUKA-internal rejection.
- `_check_kinematics(kin)` — delegates to base class for
  projectile/target/ecm_min; additionally enforces per-nucleon upper bounds
  (hadron vs photon); raises `KeyError` with remediation hint if the target
  is outside the registered material list.

`CrossSectionData` in `src/chromo/common.py` gets an `emd` field (default 0
for models that don't populate it). Audit existing models to ensure no
constructor-by-position break.

Default material list expanded to a pragmatic cosmic-ray + heavy-ion set:

```python
_DEFAULT_MATERIALS = [
    2212,    # free proton (Z=1, stand-in for H)
    "H1", "He4",
    "C12", "N14", "O16", "Ne20", "Ar40",
    "Fe56", "Cu63", "Ag108", "Au197", "Pb208",
]
```

## Material handling — extension strategy

1. **Construction time (always supported).** Union of `_DEFAULT_MATERIALS`,
   explicit `targets=...`, and components of `evt_kin.p2` if
   `CompositeTarget`. `_init_fluka_once` registers them all via a single
   `CHROMO_STPXYZ` call.
2. **After `_init_fluka_once` (investigation 3 required).** Probe
   whether FLUKA's lower-level material-registration routines (`SETMAT`,
   direct manipulation of `FLKMAT.ZTAR/AMSS/NMAT`) permit adding a new
   element without re-entering `STPXYZ`. If yes → expose
   `Fluka.register_target(pdg)` at runtime. If no → keep the constructor
   as the only extension point, and raise `KeyError` at `_set_kinematics`
   with this message:
   `"target {name} not initialised; pass targets=['{name}'] to the Fluka "
   "constructor."`
3. **Never re-call `STPXYZ`.** The reference source is explicit: a second
   call fires `FLABRT('STPXYZ', 'MULTIPLE CALLS!')`.

`CompositeTarget` support falls out automatically because the base class's
`_temporary_kinematics` iterates components and calls `_set_kinematics` per
component, and each component PDG is already in `self._materials_map` after
the initial registration.

## Projectile encoding

Single Fortran entry point `pdg_to_proj_code(pdg_id, proj_code)`:

```fortran
if (abs(pdg_id) .lt. 1000000000) then
   kflk = MCIHAD(pdg_id)
   if (pdg_id .eq. 22) then
      proj_code = 7                             ! photon external
   else
      proj_code = KPTOIP(kflk)                  ! FLUKA's own internal→external
   end if
else
   ! Nucleus PDG 10LZZZAAAI
   A = mod(pdg_id / 10,        1000)
   Z = mod(pdg_id / 10000,     1000)
   L = mod(pdg_id / 10000000,  10)
   proj_code = A*10 + Z*10000 + L*10000000 + 1000000000
end if
```

EVTXYZ and SGMXYZ already branch on `abs(proj_code) .ge. 1e9` and call
`PDGION`/`SETION` internally — we don't duplicate FLUKA's ion-setup logic.

## Fortran wrapper layer

Single file `src/fortran/fluka/chromo_fluka.f`. Additions to the existing
draft:

- `pdg_to_proj_code(pdg_id, proj_code)` — as above.
- `fluka_elem_properties(n_elements, z_out, a_out, mass_out)` — diagnostic:
  read `FLKMAT.ZTAR/AMSS` post-STPXYZ to let tests and the Python layer
  verify which elements are registered.
- `fluka_hepevt_summary(nhep_out, n_ejectiles, n_heavies, n_residuals)` —
  diagnostic counters derived from `HEPCMM.NHEP` plus a scan of the HEPEVT
  entries to categorise hadrons vs heavy fragments vs residual nucleus.
  Used by tests to verify remnant presence.

Fix:
- `CHROMO_SGMXYZ` typo (`CHROMO_SGMXY` → `CHROMO_SGMXYZ`).

`meson.build` `fluka_syms` list:
```
chromo_evtxyz, chromo_stpxyz, chromo_sgmxyz, chromo_fllhep,
pdg_to_proj_code, fluka_elem_properties, fluka_hepevt_summary,
icode_from_pdg, icode_from_pdg_arr, charge_from_pdg_arr, pdg_from_icode,
random_direction,
init_rng_state, load_rng_state, save_rng_state, fluka_rand,
fluka_particle_scheme
```
Removed from previous list: `fluka_key` (unused in Python). `chromo_fllhep`
stays.

## HEPEVT filling and nuclear remnants

`CHROMO_FLLHEP` already wraps FLUKA's `FLLHEP`. Implementation-time
investigation (item 1 below): confirm whether `FLLHEP` in FLUKA 2025.1's
`libflukahp.a` already populates HEPEVT with

- GENSTK ejectiles (standard particles),
- FHEAVY heavy fragments (d, t, ³He, α, evaporation residues A>4),
- RESNUC final residual nucleus.

If yes → done, `_get_charge` via `charge_from_pdg_arr` handles hadron
charges, and nuclei charges come from PDG decoding in Python
(`(pid // 10000) % 1000` for Z). If no → add `chromo_fill_remnants` that
reads the missing common blocks and appends HEPEVT entries with proper PDG
ion codes (`10LZZZAAAI` per 2019 PDG standard).

`FlukaEvent._get_charge` extends `charge_from_pdg_arr` to branch on
`|pdg| ≥ 10⁹` and compute Z directly from the PDG ion id.

## Frame handling

`_frame = EventFrame.FIXED_TARGET` declares the generator's native frame.
The chromo base class converts `CenterOfMass` inputs to equivalent
`FixedTarget` kinematics before calling `_set_kinematics` / `_generate`, and
transforms the returned event back to the user's requested frame. The only
thing our `Fluka` class must do is consume `kin.ekin` correctly (it's the
projectile kinetic energy in the target rest frame). Tests
(`test_xsec_cms_vs_ft_equivalent`, `test_event_cms_vs_ft_equivalent`)
exercise the round-trip.

## Energy bounds

- ecm_min: `0.01 * GeV` (covers GDR region for photonuclear).
- ekin upper per-nucleon:
  - hadron projectile: `20 * TeV`  (PEANUT + DPMJET default transition)
  - photon projectile: preliminary `1 * TeV`, refined by investigation 6.
- AA: upper bound applied per-nucleon (`kin.elab / max(A_proj, 1)`).
- FLUKA↔DPMJET transition uses FLUKA default (transition_energy=None).
  Users can override via `transition_energy` kwarg.

`_check_kinematics` enforces all of these with explicit error messages.

## Interaction-mode support (EMD included)

```python
class InteractionType(IntEnum):
    INELASTIC     = 1
    ELASTIC       = 10
    INELA_ELA     = 11
    EMD           = 100
    INELA_EMD     = 101
    ELA_EMD       = 110
    INELA_ELA_EMD = 111
```

Constructor default: `INELASTIC` (matches user expectation for standard
event generation). Users requesting EMD-only or combined modes get proper
behaviour routed through FLUKA's internal `ELDSEV` / EMD branch of EVTXYZ.

`cross_section()` calls `chromo_sgmxyz` up to three times (one per channel
the user requested) and packs results into `CrossSectionData`. The new
`emd` field on `CrossSectionData` is added for this purpose; all existing
models default it to zero.

## Install procedure

`scripts/install_fluka.sh`:

1. Check `$FLUPRO` is set; default suggestion `$HOME/devel/FLUKA`.
2. `mkdir -p $FLUPRO` if absent.
3. Idempotent short-circuit: if `$FLUPRO/libflukahp.a` already exists,
   skip to step 7 (verification only).
4. Locate both archives in `$HOME/devel/FLUKA-dev/`:
   `fluka2025.1-linux-gfor64bit-10.3-glibc2.32-AA.tar.gz` and
   `fluka2025.1-data.tar.gz`. Fail if missing.
5. `cd $FLUPRO && tar xzf <code> && tar xzf <data>`.
6. `cd $FLUPRO && make`.
7. Verify all five archives are present (`libflukahp.a`, `libdpmmvax.a`,
   `librqmdmvax.a`, `latestRQMD/librqmd.a`, `interface/libdpmjet*.a`).
8. Print shell `export FLUPRO=...` reminder.

Documented in a new "FLUKA" subsection of `CLAUDE.md`.

CI: opt-in. Fluka is not built in public CI unless a `FLUPRO` provider is
configured. Tests gated with `pytest.mark.skipif` on `_fluka` import.

## `meson.build` changes

- Replace `meson.project_source_root() + '/../FLUKA/interface/dpmvers'` with
  `fluprod + '/interface/dpmvers'`.
- Add presence checks for each of the five input libraries before
  `custom_target('fluka_all', ...)`; error with a clear message pointing to
  `scripts/install_fluka.sh` if missing.
- Update `fluka_syms` per the Fortran layer section.

## Test plan

File: `tests/test_fluka.py`. Structure mirrors `tests/test_pythia8_cascade.py`
and `tests/test_pythia8_angantyr.py`:
- `pytestmark = pytest.mark.skipif(not _fluka_available, reason=...)`.
- Each test runs in a separate subprocess via `run_in_separate_process`
  (FLUKA is single-instantiation per process).
- Fixtures cached with `@lru_cache(maxsize=1)`.

**Cross sections:**
- `test_xsec_p_A_sweep` — p+{H1, N14, O16, Ar40, Fe56, Pb208} at 100 GeV:
  100 mb < σ_inel < 3000 mb; monotone in A.
- `test_xsec_pi_N14` — π⁺ + N14 at 100 GeV: 200–600 mb.
- `test_xsec_p_composite_air` — CompositeTarget (air); σ within envelope
  of its components.
- `test_xsec_AA_O_O` — ¹⁶O+¹⁶O at 100 GeV/nucleon: 500–2000 mb.
- `test_xsec_cms_vs_ft_equivalent` — CMS and FT with equivalent kinematics
  give identical cross sections.

**Photohadronic (γ+p):**
- `test_xsec_gamma_p` — γ+p at 1/10/100 GeV: σ_inel in PDG-compatible range.
- `test_generate_gamma_p` — 1 event; hadrons in final state; Q_fs = 0.

**Photonuclear (γ+A, multiple regimes):**
- `test_xsec_gamma_Pb_GDR` — γ+Pb208 at 15 MeV: σ large (GDR peak).
- `test_xsec_gamma_Pb_delta` — γ+Pb208 at 300 MeV: Δ-resonance.
- `test_xsec_gamma_Pb_DIS` — γ+Pb208 at 10 GeV: DIS regime.
- `test_generate_gamma_Pb_at_delta` — 1 event; residual A close to 207.
- `test_photonuclear_energy_bounds` — γ+Pb at 10 MeV/10 GeV/100 GeV/1 TeV;
  either succeeds or raises `ValueError`. Output of investigation 6.

**EMD:**
- `test_xsec_emd_p_Pb208` — p+Pb208 with `INELA_EMD`; both σ.inel and σ.emd
  > 0; σ.emd < σ.inel.
- `test_xsec_emd_AA_O_Pb` — O16+Pb208 at 100 GeV/nucleon; large EMD
  contribution (Z² UPC enhancement).
- `test_generate_emd_event` — EMD-only event; small multiplicity; projectile
  survives; a few nucleons emitted from target.

**Events:**
- `test_generate_p_O16` — standard pA at 100 GeV ekin.
- `test_generate_pi_Fe56` — π+Fe56.
- `test_generate_AA_O_O` — O+O event; n_final > 10.
- `test_generate_composite_air` — events distribute across air components.

**Conservation (parametrised):**
- Baryon number: `sum(B(pid)) == B_in`.
- Charge: `sum(Q(pid)) == Q_proj + Q_targ`.
- Energy/momentum: `sum(p4) == p4_in` to tolerance 1e-3 · E_beam.

**Nuclear remnants:**
- `test_remnant_present_hA` — p+Pb208 at 10 GeV: at least one nucleus in
  `event.pid` with A close to target A - few.
- `test_remnant_present_AA` — O+Au197: both projectile and target
  remnants present.
- `test_heavy_fragments` — p+Fe56: d/t/α fragments occasionally appear.

**Infrastructure:**
- `test_charge_from_pdg` — `event.charge == reference_charge(event.pid)`.
- `test_rng_state_roundtrip` — save → generate N → load → generate N →
  identical event arrays.
- `test_registered_targets` — default list registered; constructor
  `targets=` extension works.
- `test_unregistered_target_error` — raising with clear remediation text.

**Energy bounds:**
- `test_below_ecm_min` — ValueError.
- `test_above_ekin_max_hN` — ValueError.
- `test_above_ekin_per_nucleon_max_AA` — ValueError.

Total: ~32 tests. Also add `Fluka` to `tests/test_generators.py` and
`tests/test_cross_sections.py` model matrices if they iterate over models.

## Implementation-time investigations

These are explicit, scoped steps in the implementation plan — not design
decisions.

1. **`FLLHEP` content.** Run one hA event; inspect HEPEVT; confirm
   whether it already contains FHEAVY and RESNUC entries. Branches the
   Fortran layer plan (add `chromo_fill_remnants` only if needed).
2. **Projectile-code convention.** Verify `pdg_to_proj_code` output via
   one proton event. Confirm the downstream FLUKA routine accepts
   `KPTOIP(MCIHAD(pdg))` and not raw `MCIHAD(pdg)`.
3. **Runtime material extension.** Probe whether a lower-level FLUKA
   routine can add an element to FLKMAT after STPXYZ. If yes → expose
   `Fluka.register_target(pdg)`.
4. **Ion decode correctness.** Smoke-test `pdg_to_proj_code` for d, ³He,
   ⁴He, ¹⁶O, ²⁰⁸Pb against direct `PDGION` call output.
5. **Frame round-trip.** Statistical equivalence between CMS input and
   equivalent FT input for the same collision.
6. **Photonuclear upper validity.** Sweep γ+Pb at 10 MeV to 1 TeV; set
   `_ekin_per_nucleon_max_photon` to the confirmed upper limit; document.
7. **EMD smoke.** Send IFLXYZ=101 (INELA_EMD) through CHROMO_EVTXYZ;
   confirm no FLABRT; validate EMD cross section is non-zero for charged
   projectile on heavy target.

## Error handling

- Unregistered target → `KeyError` with text:
  `"target {name} not initialised; pass targets=['{name}'] to the Fluka "
  "constructor."`
- Energy out of range → `ValueError` from `_check_kinematics`.
- `e+/e-` projectile → `ValueError` from `_check_kinematics` (not in
  `_projectiles`).
- Missing `$FLUPRO` at import time → actionable message pointing to
  `scripts/install_fluka.sh`.
- FLABRT from FLUKA is fatal (`STOP`). Don't attempt to catch; prevent by
  validating inputs before calling into FLUKA.

## Files changed or added

Added:
- `docs/superpowers/specs/2026-04-15-fluka-robust-integration-design.md`
  (this spec).
- `scripts/install_fluka.sh`.
- `tests/test_fluka.py`.

Modified:
- `src/chromo/models/fluka.py` — materials, projectile encoding, EMD,
  frame, bounds, composite.
- `src/chromo/models/__init__.py` — already exports `Fluka`.
- `src/fortran/fluka/chromo_fluka.f` — SGMXY typo fix, three new
  helpers, symbol list.
- `src/chromo/common.py` — add `emd` field on `CrossSectionData`.
- `meson.build` — dpmvers path, library presence checks, symbols.
- `CLAUDE.md` — FLUKA install + usage section.

No changes to other model files beyond the `CrossSectionData` default
(zero-emd).
