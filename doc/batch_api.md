# Batch generation with varying initial states (issue #188)

Minimal design note for `MCRun.generate_batch(states, *, seed=None)`.

## Motivation

The usual chromo loop `for e in model(nevents)` amortises the cost of setting
the initial state over many events. In air-shower-like use cases (CRPROPA,
hybrid shower MCs) the initial state changes every event, and model setup
dominates the runtime. Issue #188 proposed pushing the loop into Fortran and
returning awkward arrays; the maintainers noted that a cheaper, model-side
win is a *sorted input stack*: generate all configurations grouped by equal
initial state, so the generator switches beams as rarely as possible.

## Design

`generate_batch` is implemented once in `MCRun` (Python layer) rather than
per model:

- Input is a sequence of `EventKinematicsBase` states (the existing API
  object for an initial state; particle-list/4-vector input is representable
  as `EventKinematicsWithRestframe(p1, p2, beam=(pz1, pz2))`).
- Output is `list[EventData]` in input order (picklable, same element type
  as the loop API; awkward arrays are a possible later optimisation but not
  required for correctness).
- States are grouped by value into first-appearance order (dict grouping,
  deterministic across processes) so that *equal states run consecutively*
  ("sorted input stack"); the generator state is switched once per group via
  the ordinary `generator.kinematics = ...` setter, then the group is drained
  with the ordinary `generator(n)` call. Rejection, composite-target sampling
  and decay validation behave exactly as in `__call__`.
- Optional `seed` resets the RNG; for Pythia8 the seed is re-applied in
  `_set_kinematics`, since Pythia re-reads `Random:seed` on every `init()`
  and `generate_batch` re-initialises per group.

Because grouping happens internally, callers do not have to sort anything,
and behaviour is identical to the one-by-one loop except for the (documented)
RNG stream layout: duplicates of a state share one initialisation and
continue the RNG stream.

## Per-model support

Any model whose `generator.kinematics = ...` switching works after
construction can use `generate_batch`; it inherits all model-specific
restrictions of that setter:

- **Pythia8** (hN/γγ/γN): full support; re-inits per group (cheap). Tested.
- **DPMJET-III / PHOJET**: work, but all batch energies must stay below the
  initialization energy and below the initialized A of both sides
  (initialize the model at the highest energy of the batch, cf. issue #242).
  Nuclear-target batches also re-run the Glauber MC per switch. Tested.
- **UrQMD**: works; note that switching to a nuclear target recomputes the
  Glauber-like cross section table and that UrQMD keeps additional hidden
  state across kinematics switches, so batch results are only reproducible
  for a fixed seed, not identical to one-by-one runs across different
  states. Tested.
- **SIBYLL / QGSJet / EPOS**: kinematics switching is supported by the
  wrapper, so `generate_batch` works, but the beam setup rewrites several
  Fortran module variables per switch; no speed-up claim is made.
- **Pythia8Cascade / Pythia8Angantyr**: supported in principle (kinematics
  switching exists), untested; Angantyr nuclear modes re-init tables per
  switch.

## What this does *not* do (future work)

- No Fortran-side loop (point 1 of #188): the Python overhead per event
  remains, but the number of *expensive model state switches* now scales
  with the number of distinct states, matching the sorted-stack suggestion.
- No awkward-array output; `EventData.copy()` per event means one
  numpy-array copy per particle, same as `final_state()` usage today.
- A packed-table output (`EventStack`) can be added later without changing
  this signature.
