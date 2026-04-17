# Decay Handling

chromo uses a consistent definition of "final state" across all generators: particles with lifetimes greater than 30 ps are considered stable (the ALICE convention). Short-lived resonances are decayed by the generator.

## Default Behavior

By default, each generator handles decays according to its own physics model. The `final_state()` method returns particles with `status == 1`, which are the prompt, long-lived particles.

## Pythia8 Decay Handler

For generators that do not decay all short-lived particles natively (e.g., some produce undecayed charm or strange hadrons), chromo can use Pythia 8 as an auxiliary decay engine:

```python
from chromo.decay_handler import Pythia8DecayHandler
from chromo.constants import long_lived

handler = Pythia8DecayHandler(stable_pids=long_lived, seed=42)
```

The `Pythia8DecayHandler` is automatically activated for generators that need it. You typically do not need to configure it manually.

!!! note
    Some generators (e.g., FLUKA) have their own global decay model that is not runtime-configurable. For FLUKA, `_set_stable` is a no-op.

## Controlling Stable Particles

Some models support `set_stable()` to control which particles are treated as stable during generation:

```python
generator.set_stable(pid, stable=True)   # make particle stable
generator.set_stable(pid, stable=False)  # allow particle to decay
```

Not all models support this. Check the [Model Overview](../models/overview.md) for per-model capabilities.

## Reproducibility Note

When the Pythia8DecayHandler is active, its internal RNG is independent from the main generator's RNG. For fully reproducible event records (e.g., with FLUKA), deactivate it:

```python
generator._activate_decay_handler(on=False)
```

See the [Decay Handler notebook](../examples/decayhandler.ipynb) for a detailed example.
