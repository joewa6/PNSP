# Stage Pipeline

Each stage is a **pure transformation** with defined inputs/outputs.

---

## Stage A: Builder

**Purpose**: Generate simulation input files from parameter sets.

**Inputs**:
- Parameter specification (temperature, pressure, composition, etc.)
- Force field selection
- Topology templates

**Outputs**:
- `run_manifest.json` (created with `run_id`, parameters, builder metadata)
- Topology files (`.prmtop`, `.inpcrd`, or equivalent)
- Input configurations (`.in` files for GROMACS/Amber/LAMMPS)

**Side effects**:
- Creates working directory under `stageA_builder/builds/{run_id}/`
- Logs build process

**Contract**:
```json
{
  "run_id": "unique-identifier",
  "created_at": "ISO-8601-timestamp",
  "parameters": { ... },
  "stageA": {
    "status": "complete",
    "builder_version": "1.2.3",
    "topology_path": "relative/path/to/topology",
    "completed_at": "ISO-8601-timestamp"
  }
}
```

---

## Stage B: Simulator

**Purpose**: Execute molecular dynamics simulation.

**Inputs**:
- `run_manifest.json` (from Stage A, with topology references)
- Simulation engine config (local vs HPC, resource requests)

**Outputs**:
- Updated `run_manifest.json` with simulation metadata
- Trajectory files (`.nc`, `.xtc`, `.dcd`)
- Energy/log files (`.out`, `.log`, `.mdout`)
- Checkpoint/restart files

**Side effects**:
- Creates scratch/work directory during execution
- Submits jobs to SLURM (if HPC mode)
- Writes large trajectory files to designated output location

**Contract** (manifest update):
```json
{
  "stageB": {
    "status": "complete",
    "engine": "amber22",
    "environment": "hpc-cluster-01",
    "job_id": "slurm-12345",
    "started_at": "...",
    "completed_at": "...",
    "trajectory_path": "scratch/output/prod.nc",
    "total_steps": 5000000
  }
}
```

---

## Stage C: Archive

**Purpose**: Freeze evidence to immutable storage and create pointer views.

**Inputs**:
- `run_manifest.json` (from Stage B, with trajectory paths)
- Access to WD4/WD8/HPC archive mounts

**Outputs**:
- Frozen archive directory on WD4/WD8 (timestamped, immutable)
- Updated `run_manifest.json` with archive pointers
- Symlinks under `stageC_archive/views/runs/{run_id}/` pointing to frozen data

**Side effects**:
- Copies trajectories/logs/topologies to cold storage
- Creates symlink forest for browsing
- **Never** modifies frozen data after write

**Contract** (manifest update):
```json
{
  "stageC": {
    "status": "archived",
    "archive_location": "WD8/MSE_LAB/runs/2026-01-13/run-abc123/",
    "archived_at": "...",
    "frozen_files": [
      "prod.nc",
      "topology.prmtop",
      "prod.out"
    ]
  }
}
```

**Rebuild capability**:
```bash
python stageC_rebuild_views.py --from-archive WD8/MSE_LAB
# Reconstructs entire views/runs/ forest from frozen manifests
```

---

## Stage D: Reduce

**Purpose**: Compute metrics, aggregations, and analysis tables from frozen runs.

**Inputs**:
- Frozen evidence (via `stageC_archive/views/runs/`)
- Analysis specification (what metrics to extract)
- `run_manifest.json` for metadata

**Outputs**:
- Metric tables (Parquet, CSV, HDF5)
- Summary statistics
- Optional plots/visualizations
- Updated `run_manifest.json` with reduction metadata

**Side effects**:
- Caches intermediate results in `stageD_reduce/cache/`
- Writes final outputs to `stageD_reduce/outputs/`
- **Never** modifies frozen data (read-only access)

**Contract** (manifest update):
```json
{
  "stageD": {
    "status": "complete",
    "metrics_computed": ["rmsd", "rg", "energy"],
    "output_table": "outputs/metrics_2026-01-13.parquet",
    "computed_at": "..."
  }
}
```

---

## Stage transitions (managed by Composer)

```
  NEW
   ↓
[Stage A] → BUILT
   ↓
[Stage B] → SIMULATED
   ↓
[Stage C] → ARCHIVED
   ↓
[Stage D] → REDUCED
```

Only Composer decides when a run can transition.

Stages are **pure workers** — they take a manifest, do work, return an updated manifest, and exit.

---

## Key invariants

1. **Stages never mutate Composer state directly**
2. **Manifests are append-only** (stages add sections, never delete)
3. **Frozen evidence is immutable** (Stage C writes once, never modifies)
4. **Stages can be re-run** (idempotent where possible, or fail-fast if not)
5. **No stage depends on another stage's code** (only on manifest schema)
