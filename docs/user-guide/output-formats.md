# Output Formats

chromo supports multiple output formats for generated events.

## HepMC3

The recommended format for interoperability. Requires `pyhepmc`:

```bash
pip install pyhepmc
```

### From Python

```python
event.to_hepmc3()  # returns a pyhepmc.GenEvent
```

### From the CLI

```bash
# Write HepMC3 text output
chromo -m Sibyll23d -n 1000 -o events.hepmc

# Write gzip-compressed HepMC3
chromo -m Sibyll23d -n 1000 -o events.hepmc.gz
```

HepMC output can be piped directly into [RIVET](https://rivet.hepforge.org/) and other tools.

## ROOT

Write events to ROOT files via `uproot` (no ROOT installation needed):

```bash
pip install uproot awkward
```

```bash
chromo -m Sibyll23d -n 1000 -o events.root
```

See the [ROOT Output example notebook](../examples/write_root.ipynb) for details.

## SVG Visualization

Visualize individual events as SVG images (requires `pyhepmc`):

```python
hepmc_event = event.to_hepmc3()
# In Jupyter, the event renders automatically as HTML/SVG
hepmc_event
```

See the [HepMC IO & Visualization notebook](../examples/hepmc_io_and_visualization.ipynb) for examples.

## Command-Line Interface

The CLI is designed to be familiar for users of [CRMC](https://gitlab.iap.kit.edu/AirShowerPhysics/crmc):

```bash
chromo --help

# Basic usage
chromo -m Sibyll23d -n 1000 -o output.hepmc \
    --ecm 13000 -p 2212 -P 2212

# Fixed-target mode
chromo -m EposLHC -n 100 -o output.hepmc \
    --elab 1000000 -p 2212 -P 1000070140
```

If `chromo` is not on your PATH, use `python -m chromo` instead.
