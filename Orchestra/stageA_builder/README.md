# Stage A Builder

Production-ready builder for creating NP-ligand complex systems with Epik microstate enumeration.

## Files

- **`build.py`** - **PRODUCTION** Main builder with complete Epik + NP complex pipeline

## Usage

### Single Ligand Build with NP Complex

```bash
python src/build.py --smiles "CCO" --output /path/to/output
```

**What it does:**
1. Uses Epik to enumerate protonation/tautomeric states
2. Selects the lowest-penalty microstate
3. Generates 3D coordinates
4. Prepares AMBER topology with antechamber/parmchk2/tleap
5. **Builds complete NP complex:**
   - Loads nanoparticle structure from NP library
   - Packs 10 ligand copies around NP using Packmol
   - Solvates in TIP3P water box
   - Adds neutralizing ions
   - Generates final AMBER system files
6. Creates run manifest with metadata

**Output:**
- `system.prmtop` - AMBER topology file (NP + ligands + water + ions)
- `system.inpcrd` - AMBER coordinates file
- `run_manifest.json` - Metadata and provenance
- `LIG.mol2`, `LIG.frcmod`, `LIG.lib` - Ligand parameter files
- `NP_clean.pdb` - Cleaned NP structure
- `combined.pdb` - NP + ligands from Packmol
- `ligands_only.pdb` - Extracted ligand coordinates

**Fallback behavior:**
- If NP files not found: builds ligand-only system
- If Epik unavailable: creates placeholder files

### Batch Processing Mode (CSV)

```bash
python src/build.py --csv molecules.csv --output-dir ./batch_output
```

**CSV Format:**
```csv
Name,SMILES
Ethanol,CCO
Diuron,CC(C)Nc1nc(Cl)nc(NC(C)C)n1
Caffeine,CN1C=NC2=C1C(=O)N(C(=O)N2C)C
```

**What it does:**
- Processes multiple molecules from CSV file
- Creates separate directory for each molecule
- Handles errors gracefully (continues on failure)
- Generates detailed batch summary report

**Output Structure:**
```
batch_output/
├── Ethanol/
│   ├── system.prmtop
│   ├── system.inpcrd
│   └── run_manifest.json
├── Diuron/
│   ├── system.prmtop
│   ├── system.inpcrd
│   └── run_manifest.json
├── Caffeine/
│   ├── system.prmtop
│   ├── system.inpcrd
│   └── run_manifest.json
└── batch_summary.json
```

**Batch Summary (`batch_summary.json`):**
- Total molecules processed
- Success/failure/skipped counts
- Detailed error messages for failures
- Output directory paths
- Timestamp and git commit info

### Dry Run

```bash
python src/build.py --smiles "CCO" --output /path/to/output --dry-run
```

Creates manifest only, skips topology/coordinate generation.

## Requirements

- Python 3.11+
- RDKit
- Schrödinger Suite (epikx, structconvert)
- AmberTools (antechamber, parmchk2, tleap)
- Packmol
- OpenBabel
- AmberTools (antechamber, parmchk2, tleap)
- OpenBabel (obabel)

## Architecture

The refactoring separated two distinct systems:

1. **Single-ligand builder** (`build_test.py`) - Simple SMILES → topology workflow
2. **Batch complex builder** (`epik_builder.py`) - CSV → multiple NP complexes

Both share the same Epik-based ligand preparation pipeline.
