# Evidence Policy

## Core principle

> **Evidence is sacred. Control is versioned.**

---

## What is evidence?

Evidence is any file that represents **irreversible experimental/computational work**:

- Molecular dynamics trajectories (`.nc`, `.xtc`, `.dcd`, `.trr`)
- Topology/structure files (`.prmtop`, `.pdb`, `.top`, `.gro`)
- Simulation logs and output (`.out`, `.log`, `.mdout`, `.mdinfo`)
- Restart/checkpoint files (`.rst`, `.rst7`, `.crd`)
- Analysis outputs (`.parquet`, `.csv`, `.h5`, metrics tables)
- Intermediate scratch (build artifacts, temporary work directories)

---

## Evidence NEVER goes in GitHub

Why:

1. **Size**: Trajectories are GBs to TBs. GitHub repos become unmaintainable.
2. **Semantics**: Evidence is not "versioned" — it's **frozen at creation**.
3. **Pollution**: Committing data obscures the actual control plane changes.
4. **Irreversibility**: You can't "fix" a trajectory by editing it. You re-run the simulation.

GitHub stores the **recipe** (how to create evidence).  
WD4/WD8/HPC store the **evidence** (the result of running that recipe).

---

## Where evidence lives

### During execution (ephemeral)

- **Stage A builds**: `Orchestra/stageA_builder/builds/{run_id}/`
- **Stage B scratch**: `Orchestra/stageB_simulator/scratch/{run_id}/`

These are **temporary**. They are `.gitignore`d. They are deleted after archival.

---

### After archival (permanent)

- **WD4**: High-speed NVMe archive (recent runs, hot access)
- **WD8**: Long-term spinning disk archive (completed campaigns)
- **HPC scratch**: Temporary HPC storage (purged after retrieval)

Frozen structure:
```
WD8/MSE_LAB/
  runs/
    2026-01-13/
      run-abc123/
        run_manifest.json
        topology.prmtop
        prod.nc
        prod.out
```

Once written, **never modified**.

---

### Pointer views (symlinks)

`Orchestra/stageC_archive/views/runs/{run_id}/` contains **symlinks** to frozen evidence.

Example:
```
views/runs/run-abc123/
  run_manifest.json → /mnt/WD8/MSE_LAB/runs/2026-01-13/run-abc123/run_manifest.json
  prod.nc → /mnt/WD8/MSE_LAB/runs/2026-01-13/run-abc123/prod.nc
```

The symlinks are **not tracked by git** (`.gitignore` blocks `views/runs/`).

They are **rebuilt on demand**:
```bash
python stageC_rebuild_views.py
```

This design lets you:
- Browse runs as if they're local
- Keep GitHub clean
- Reconstruct the view forest from manifests alone

---

## What IS tracked by git

Only the **control plane**:

- Stage code (builders, simulators, archivers, reducers)
- Manifest schemas (`run_schema.yaml`)
- Composer ledger logic
- Configuration templates
- Documentation

These are the **laws of physics** for your lab.

Evidence is the **fossils** created by those laws.

---

## The `.gitignore` enforcement

The `.gitignore` file is **load-bearing infrastructure**.

It blocks:
1. All common MD file formats
2. Stage scratch directories
3. Archive views
4. Analysis outputs and caches

**Never edit `.gitignore` to allow evidence files.**

If you find yourself wanting to commit a `.nc` or `.prmtop` file, you are doing it wrong.

Instead:
- Commit the **template** or **builder script** that creates it
- Let Stage A/B generate it during execution
- Let Stage C archive it to WD4/WD8

---

## GitHub Actions enforcement

The `no-evidence` workflow hard-fails if evidence files are committed.

It checks for:
- File extensions (`.nc`, `.prmtop`, `.parquet`, etc.)
- Blocked directories (`views/runs/`, `cache/`, `outputs/`)

This is your **safety net** when you're tired and reckless.

---

## Rebuild guarantee

At any time, you can:

1. Clone `joewa6/MSE_LAB` from GitHub
2. Mount WD4/WD8/HPC archives
3. Run `stageC_rebuild_views.py`
4. Get back the complete symlink forest

No data is duplicated. No history is lost.

**GitHub + WD4/WD8 = complete lab state.**

---

## Summary

| Type | Location | Git tracked? | Mutable? |
|------|----------|--------------|----------|
| Stage code | `Orchestra/stageA_builder/` etc. | ✅ YES | ✅ YES (via commits) |
| Templates | `Orchestra/*/templates/` | ✅ YES | ✅ YES |
| Build scratch | `Orchestra/stageA_builder/builds/` | ❌ NO | ✅ YES (then deleted) |
| Frozen evidence | `WD4/WD8/MSE_LAB/runs/` | ❌ NO | ❌ NO (immutable) |
| Pointer views | `Orchestra/stageC_archive/views/runs/` | ❌ NO | ✅ YES (rebuilt on demand) |
| Analysis cache | `Orchestra/stageD_reduce/cache/` | ❌ NO | ✅ YES (ephemeral) |

---

**If it's data, it's frozen externally.**  
**If it's code, it's versioned in GitHub.**

No exceptions.
