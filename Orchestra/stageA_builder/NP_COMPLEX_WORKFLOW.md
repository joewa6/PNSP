# NP Complex Building Workflow

## Overview

The `build_test.py` script now includes complete NP complex assembly with Epik microstate enumeration.

## Full Pipeline

```
SMILES Input
    ↓
[1] Epik Microstate Enumeration (pH 7.4±2.0)
    ↓
[2] Ligand Parameterization (GAFF2, BCC charges)
    ↓
[3] NP Complex Assembly
    ├── Load NP structure from library
    ├── Calculate bounding box
    ├── Pack 10 ligand copies with Packmol
    └── Combine NP + ligands
    ↓
[4] System Solvation & Finalization
    ├── Solvate in TIP3P water (20Å buffer)
    ├── Add neutralizing Na+ ions
    └── Generate AMBER topology
    ↓
Final Output: system.prmtop, system.inpcrd
```

## Command

```bash
python src/build_test.py --smiles "CCO" --output /path/to/output
```

## Output Files

### Primary Outputs
- `system.prmtop` - AMBER topology (NP + ligands + water + ions)
- `system.inpcrd` - AMBER coordinates
- `run_manifest.json` - Metadata and provenance

### Ligand Files
- `lig_XXXXXXXX.sdf` - Epik best microstate (SDF format)
- `lig_XXXXXXXX.mol2` - Epik best microstate (MOL2 format)
- `LIG.mol2` - GAFF2-typed ligand
- `LIG.frcmod` - Missing force field parameters
- `LIG.lib` - AMBER library file

### NP Complex Files
- `NP_clean.pdb` - Cleaned NP structure from library
- `combined.pdb` - NP + ligands from Packmol
- `ligands_only.pdb` - Extracted ligand coordinates
- `system_solvated.pdb` - Final solvated system

### Intermediate Files
- `packmol.inp` - Packmol input script
- `tleap_*.in` - tleap input scripts
- `leap.log` - tleap execution log

## NP Complex Assembly Details

### Step 1: Create NP_clean.pdb

```bash
tleap -f tleap_np_clean.in
```

Handles two NP library formats:
- **Format A**: Top-level `NP` unit → `saveamberparm NP NP_clean.pdb`
- **Format B**: Multiple units (F04, F06, etc.) → `list` to find first unit

### Step 2: Calculate Bounding Box

```bash
grep "^ATOM\|^HETATM" NP_clean.pdb | awk '{print $6, $7, $8}'
```

Extracts X/Y/Z coordinates, finds min/max, adds 20Å buffer for ligands.

### Step 3: Packmol

```
tolerance 2.0
filetype pdb
output combined.pdb

structure NP_clean.pdb
  number 1
  fixed xmin ymin zmin 0. 0. 0.
end structure

structure LIG.pdb
  number 10
  inside box xmin ymin zmin xmax ymax zmax
end structure
```

Places 10 ligand copies randomly in box around fixed NP.

### Step 4: Final tleap

```
source leaprc.protein.ff19SB
source leaprc.gaff2
source leaprc.water.tip3p
loadamberparams frcmod.ionsjc_tip3p

loadamberparams NP_charged.frcmod
loadoff NP_correct.lib
loadamberparams LIG.frcmod
loadoff LIG.lib

SYS = loadpdb ligands_only.pdb
combine {NP SYS}

solvatebox SYS TIP3PBOX 20.0
addions SYS Na+ 0

saveamberparm SYS system.prmtop system.inpcrd
savepdb SYS system_solvated.pdb
quit
```

## Fallback Hierarchy

1. **NP Complex** (preferred)
   - NP files found in `../../NP` or `../NP` or `../PET_NP`
   - Full pipeline: Epik → ligand prep → NP assembly → solvation

2. **Ligand-Only System** (fallback #1)
   - NP files not found
   - Builds solvated ligand system without NP

3. **Placeholder Files** (fallback #2)
   - Epik unavailable or critical error
   - Creates minimal manifest with placeholder topology/coords

## System Requirements

- **Schrödinger Suite**: epikx, structconvert
- **AmberTools**: antechamber, parmchk2, tleap
- **Packmol**: Molecular packing
- **OpenBabel**: obabel for SDF→MOL2 conversion
- **RDKit**: Molecular handling (Python)

## NP Directory Structure

Expected layout:
```
NP/
├── NP_correct.lib      # AMBER library with NP residues
├── NP_charged.frcmod   # Force field modifications (charged)
├── NP_neut.frcmod      # Force field modifications (neutral)
└── NP_dry.pdb          # Reference structure (optional)
```

## Typical System Size

For a simple ligand like ethanol (CCO):
- NP: ~23,000 atoms
- 10 ligands: 90 atoms (9 atoms × 10 copies)
- Water: ~27,500 molecules (~82,500 atoms)
- Ions: ~50 Na+ for neutralization
- **Total: ~106,000 atoms**

## Performance

Single ligand build time: ~2-3 minutes
- Epik enumeration: 30s
- Ligand parameterization: 20s
- NP preparation: 10s
- Packmol: 30s
- Solvation: 60s
