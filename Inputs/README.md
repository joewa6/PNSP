# Inputs Directory

**Human-curated, version-controlled input data for pipeline runs.**

## Principles

- **Inputs are source data** - Hand-edited, reviewed, git-tracked
- **Treat as immutable** - Don't generate outputs here
- **Small, curated files** - For large datasets, store pointers + checksums
- **Profile-based configuration** - YAML profiles for different run configurations

## Structure

```
Inputs/
├── stageA/
│   ├── ligands/              # Molecule CSV files
│   │   └── *.csv             # Name,SMILES format
│   └── nanoparticle/         # NP profiles
│       ├── NP_profile_c50.yaml
│       ├── NP_c50/           # Canonical NP files for -50e charge
│       │   ├── NP_correct.lib
│       │   ├── NP_charged.frcmod
│       │   ├── NP_neut.frcmod
│       │   └── NP_dry.pdb
│       └── (future profiles)
└── stageB/
    └── sim_profiles/         # MD simulation profiles
        └── *.yaml            # Simulation parameters
```

## Stage A Inputs

### Ligands

CSV files with molecule definitions:
- Required columns: `Name`, `SMILES`
- One molecule per row
- Example: `PET_screening_v1.csv`

### Nanoparticle Profiles

YAML configuration files defining NP charge states:
- `NP_profile_c50.yaml` - Standard -50e charge state
- Each profile points to force field files in subdirectory

## Stage B Inputs

### Simulation Profiles

YAML files defining MD simulation parameters:
- Temperature, pressure, duration
- Equilibration protocol
- Production run settings

## Usage

Reference inputs from pipeline stages:

```bash
# Stage A
python Orchestra/stageA_builder/src/build.py \
  --csv Inputs/stageA/ligands/PET_screening_v1.csv \
  --np-profile Inputs/stageA/nanoparticle/NP_profile_c50.yaml \
  --run-root Runs/PET_SCREENING/

# Stage B  
python Orchestra/stageB_simulator/src/simulate.py \
  --run-dir Runs/PET_SCREENING/<run_id> \
  --sim-profile Inputs/stageB/sim_profiles/prod_10ns.yaml
```

## What Does NOT Belong Here

- ❌ Generated outputs (goes in `Runs/`)
- ❌ Test run data (goes in `Runs/__SMOKE__/` or `tests/fixtures/`)
- ❌ Logs, temporary files (goes in `Runs/`)
- ❌ Large binary datasets (store pointer + checksum instead)
- ❌ Stage code (goes in `Orchestra/`)
