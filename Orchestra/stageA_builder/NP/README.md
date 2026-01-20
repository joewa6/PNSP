# NP Build Scripts

**Legacy directory containing tleap build scripts for nanoparticle preparation.**

## Purpose

This directory contains the original tleap input scripts used to create the canonical nanoparticle force field files that are now stored in `Inputs/stageA/nanoparticle/`.

## Contents

- `make_dry.in` - tleap script to generate dry NP structure
- `check.in` - tleap script to verify NP structure
- `check_charge.i` - Script to calculate total charge
- `NP_*.lib`, `NP_*.frcmod`, `NP_*.pdb` - Canonical force field files (duplicated from Inputs/)

## Migration Note

The **canonical nanoparticle files** have been moved to:
```
Inputs/stageA/nanoparticle/NP_c50/
├── NP_correct.lib
├── NP_charged.frcmod
├── NP_neut.frcmod
└── NP_dry.pdb
```

With configuration in:
```
Inputs/stageA/nanoparticle/NP_profile_c50.yaml
```

## Usage

The build.py script now reads NP files from `Inputs/` via the `--np-profile` argument, not from this directory.

This directory is kept for:
1. **Provenance** - Shows how canonical files were originally created
2. **Regeneration** - If NP parameters need to be rebuilt with different charge states

## What Does NOT Belong Here

- ❌ Generated artifacts (leap.log, NP_solv.*, prod.rst) - REMOVED
- ❌ Old charge variants (OLD_-34/, OLD_-50/) - REMOVED  
- ❌ Run outputs - Go to `Runs/`
