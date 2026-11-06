# PNSP: PET Nanoplastic Contaminant Screening

PNSP is the computational workflow used to investigate the association of
environmentally relevant small-molecule contaminants with an atomistic
polyethylene terephthalate nanoplastic (PET-NPP) model.

Associated paper: *Molecular Determinants of Contaminant Association with PET
Nanoplastics Revealed by Atomistic Simulations*,
[doi:10.1021/acs.jcim.6c02127](https://doi.org/10.1021/acs.jcim.6c02127).

The workflow supports:

- contaminant preparation and GAFF2/AM1-BCC parameterization;
- construction of PET-NPP/contaminant simulation systems;
- preparation and execution of the AMBER molecular dynamics protocol;
- calculation of copy-resolved surface-contact occupancy and contact-event durations;
- aggregation of replicated simulations into contaminant-level statistics; and
- physicochemical descriptor analysis using ridge regression and generalized additive models.

The repository also contains the processed numerical data and scripts used to
reproduce the figures and tables reported in the associated manuscript. Raw
molecular dynamics trajectories are not included.

## Installation

Create the Conda environment:

```bash
conda env create -f environment.yml
conda activate pnsp
```

The workflow uses AmberTools for system preparation and analysis. Production
simulations in the study were performed with AMBER `pmemd.cuda_SPFP`.

Contaminant microstate selection in the manuscript used Schrödinger Epik at
pH 7.4 ± 2.0. Schrödinger software is not distributed with this repository.

## Quick start

A minimal example containing pyrene, phenol, and ethylene glycol is provided in
`examples/`.

```bash
bash examples/run_example.sh
```

This prepares contaminant parameters, constructs a ten-copy contaminant/PET-NPP
system, and generates the AMBER input files.

The resulting simulation can then be run using:

```bash
bash example_output/phenol/rep1/run_md.sh
```

## Core workflow

### 1. Prepare contaminants

```bash
python pnsp/build_contaminant.py \
    examples/example_contaminants.csv \
    --output build/ligands \
    --microstate-source input
```

For Epik-based microstate selection:

```bash
python pnsp/build_contaminant.py \
    contaminants.csv \
    --output build/ligands \
    --microstate-source epik \
    --ph 7.4 \
    --pht 2.0
```

### 2. Build the PET-NPP/contaminant system

```bash
python pnsp/build_system.py \
    --pet examples/example_pet \
    --ligand build/ligands/phenol \
    --output build/systems/phenol/rep1 \
    --copies 10 \
    --replica rep1 \
    --seed 314159
```

The system is solvated in TIP3P water and neutralized with counterions.

### 3. Prepare or run molecular dynamics

Prepare the 15 ns protocol:

```bash
python pnsp/run_md.py \
    build/systems/phenol/rep1 \
    --production-ns 15
```

Run directly when a compatible AMBER engine is available:

```bash
python pnsp/run_md.py \
    build/systems/phenol/rep1 \
    --production-ns 15 \
    --engine pmemd.cuda_SPFP \
    --execute
```

### 4. Calculate occupancy and contact events

```bash
python pnsp/analyse_occupancy.py \
    --topology build/systems/phenol/rep1/system.prmtop \
    --trajectory build/systems/phenol/rep1/prod.nc \
    --system phenol \
    --replica rep1 \
    --output results/phenol/rep1
```

A contaminant copy is classified as surface-associated when its minimum
heavy-atom distance to the PET-NPP is below 4.0 Å.

### 5. Aggregate the screen

```bash
python analysis/make_screen_summary.py \
    results \
    --expected-replicas 3 \
    --expected-copies 10 \
    --output screen_summary
```

### 6. Descriptor analysis

```bash
python analysis/make_descriptor_models.py \
    --summary screen_summary/ranked_screen_summary.csv \
    --compounds contaminants.csv \
    --output descriptor_models
```

## Paper reproduction

The `paper/` directory contains:

- processed numerical data underlying the manuscript figures and tables;
- scripts for regenerating the reported numerical figures and tables; and
- a checksum manifest for the released numerical datasets.

See [`paper/README.md`](paper/README.md) for reproduction instructions.

Representative all-atom production snapshots of phenol and ammonium systems
are provided in [`examples/snapshots/`](examples/snapshots/README.md), together
with matching topologies, ligand parameters, original MD input files, and
provenance. Solvent and counterions are retained; full trajectories are not
distributed in this repository.

## Execution metadata

Each preparation, MD and trajectory-analysis stage writes a JSON manifest with
its settings and file checksums. The `provenance` record identifies the stage,
contaminant, replica (where applicable), UTC recording time, executing Git
commit, worktree state, and input/output paths. Pass `--replica` when building
a system; MD preparation inherits its identifiers from the system manifest.
Unavailable identifiers or Git metadata are recorded as `null`, not inferred.
These records describe new executions, not reconstructed historical runs.
Generated manifests contain local paths; review them before sharing.

Run the small CPU-only tests with `python -m unittest discover -s tests -v`.

## Repository structure

```text
PNSP/
├── pnsp/        Core simulation and analysis workflow
├── analysis/    Screen-level statistical and descriptor analyses
├── examples/    Minimal example inputs
├── paper/       Processed manuscript data and reproduction scripts
├── environment.yml
└── LICENSE
```

## Scope

PNSP is intended for comparative screening of contaminant surface association
within the PET-NPP model used in the associated study. The calculated
occupancies are finite-trajectory comparative measurements and should not be
interpreted as equilibrium adsorption constants or absolute kinetic rate
constants.

## License

MIT License. See [`LICENSE`](LICENSE).
