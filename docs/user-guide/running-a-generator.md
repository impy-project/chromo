# Running a Generator

## The Generator Lifecycle

1. **Import** a model class from `chromo.models`
2. **Create** an instance with kinematics (and optional seed)
3. **Iterate** over events using the call syntax

```python
import chromo

kin = chromo.kinematics.CenterOfMass(100 * chromo.constants.GeV, "proton", "proton")
generator = chromo.models.Sibyll23d(kin, seed=42)

for event in generator(1000):
    # process each event
    pass
```

The `seed` parameter controls the random number generator for reproducibility.

## Single-Instantiation Constraint

!!! warning "One model per process"
    Most event generators are written in Fortran and use **global state** (COMMON blocks). This means you can only create **one instance of each model family per Python process**. Attempting to create a second instance will raise an `AssertionError`.

This is a fundamental limitation of the underlying Fortran codes, not a chromo bug. It means:

- You **cannot** create two `Sibyll23d` instances in the same process
- You **can** create one `Sibyll23d` and one `EposLHC` in the same process (different model families)
- You **can** reuse a single instance with different kinematics by assigning to `generator.kinematics`

## Running Multiple Models

To compare models, run each in a separate process:

```python
import multiprocessing as mp
import chromo

def run_model(model_class):
    kin = chromo.kinematics.CenterOfMass(
        100 * chromo.constants.GeV, "proton", "proton"
    )
    gen = model_class(kin)
    results = []
    for event in gen(100):
        fs = event.final_state_charged()
        results.append(len(fs))
    return model_class.__name__, results

# Run in separate processes
with mp.Pool(3) as pool:
    for name, counts in pool.map(run_model, [
        chromo.models.Sibyll23d,
        chromo.models.EposLHC,
        chromo.models.QGSJetII04,
    ]):
        print(f"{name}: mean multiplicity = {sum(counts)/len(counts):.1f}")
```

Alternatively, use `subprocess` for complete isolation:

```python
import subprocess, sys

script = '''
import chromo
gen = chromo.models.Sibyll23d(
    chromo.kinematics.CenterOfMass(100 * chromo.constants.GeV, "proton", "proton")
)
for event in gen(10):
    print(len(event.final_state()))
'''

result = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True)
print(result.stdout)
```

## Event Iteration

The generator call returns a Python generator (iterator):

```python
# Generate exactly 1000 events
for event in generator(1000):
    pass

# Access the total number of generated events
print(generator.nevents)
```

If the underlying model rejects an event (e.g., failed kinematics), chromo retries automatically. After 1000 consecutive failures, it raises a `RuntimeError`.

## Random Seeds

```python
# Fixed seed for reproducibility
gen = chromo.models.Sibyll23d(kin, seed=42)

# Access the seed
print(gen.seed)

# Access/restore full RNG state
state = gen.random_state
# ... generate some events ...
gen.random_state = state  # restore to earlier state
```
