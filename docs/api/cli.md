# chromo.cli

Command-line interface compatible with [CRMC](https://gitlab.iap.kit.edu/AirShowerPhysics/crmc).

## Usage

```bash
chromo [options]
```

If `chromo` is not on your PATH, use `python -m chromo` instead.

## Options

| Flag | Description |
|------|-------------|
| `-m MODEL` | Model name (e.g., `Sibyll23d`, `EposLHC`, `QGSJetII04`) |
| `-n N` | Number of events to generate |
| `-o FILE` | Output file (`.hepmc`, `.hepmc.gz`, `.root`) |
| `--ecm E` | Center-of-mass energy in GeV |
| `--elab E` | Lab-frame energy in GeV |
| `--plab P` | Lab-frame momentum in GeV |
| `-p PID` | Projectile PDG ID |
| `-P PID` | Target PDG ID |
| `-s SEED` | Random seed |
| `--list-models` | List available models |

## Examples

```bash
# Proton-proton at 13 TeV CMS, 1000 events, HepMC output
chromo -m Sibyll23d -n 1000 -o events.hepmc --ecm 13000 -p 2212 -P 2212

# Fixed-target: 1 PeV proton on nitrogen
chromo -m EposLHC -n 100 -o events.hepmc --elab 1000000 -p 2212 -P 1000070140

# Gzip-compressed output
chromo -m QGSJetII04 -n 500 -o events.hepmc.gz --ecm 100 -p 2212 -P 2212
```
