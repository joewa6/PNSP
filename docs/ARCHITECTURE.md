# MSE_LAB Architecture

## Core principles

> **GitHub = control plane**  
> **WD4/WD8/HPC = freezer**

Everything in this repository must be **rebuildable from scratch** except the frozen evidence living on external storage.

---

## System components

### Composer (control plane)

The state machine and orchestration layer:

- **State**: Ledger of all runs, their lineage, and current status
- **Rules**: What can transition, when, and how
- **CLI**: Human interface for declaring intent and querying state

Composer **never** touches data files directly. It only:
- Reads/writes manifests
- Updates ledger state
- Invokes Orchestra stages

---

### Orchestra (data plane)

The workers that create and transform evidence:

- **Stage A (builder)**: Constructs simulation inputs
- **Stage B (simulator)**: Executes MD runs (local/HPC)
- **Stage C (archive)**: Freezes evidence + creates pointer views
- **Stage D (reduce)**: Generates metrics/tables from frozen runs

Each stage:
1. Reads a `run_manifest.json`
2. Performs work
3. Updates the manifest with new metadata
4. Exits cleanly

---

## Data flow

```
Composer ledger
    ↓
Stage A: creates inputs → run_manifest.json
    ↓
Stage B: runs simulation → updates manifest
    ↓
Stage C: freezes to WD4/WD8 → creates symlink views
    ↓
Stage D: reads frozen evidence → emits metrics
    ↓
Composer: updates ledger with results
```

---

## Why Composer and Orchestra are separate

**Separation of concerns**:

- Composer = **what** to do (intent, state, ordering)
- Orchestra = **how** to do it (execution, environments, resources)

This prevents:
- Stage scripts from mutating ledger state directly
- Tight coupling between control logic and execution environments
- Accidental data commits (stages never write to Composer/)

---

## The contract: `run_manifest.json`

Every run is defined by a single manifest file containing:

- **Identity**: `run_id`, `lineage`, `created_at`
- **Inputs**: parameter set, topology files, force fields
- **Stage metadata**: What each stage did, when, and where
- **Pointers**: Paths to frozen evidence (never the data itself)

The manifest is the **immutable record** of what happened to a run.

Stages append to it. Composer reads it. Nobody deletes it.

---

## What lives where

| Location | Purpose | Tracked by git? |
|----------|---------|----------------|
| `Composer/` | State, rules, CLI | ✅ YES |
| `Orchestra/stageA_builder/` | Build scripts + templates | ✅ YES |
| `Orchestra/stageB_simulator/` | Runner wrappers (local/HPC) | ✅ YES |
| `Orchestra/stageC_archive/` | Archival logic + view builders | ✅ YES |
| `Orchestra/stageD_reduce/` | Analysis scripts | ✅ YES |
| `Orchestra/stageC_archive/views/runs/` | Symlink forest (evidence pointers) | ❌ NO |
| `Orchestra/stageD_reduce/cache/` | Intermediate tables | ❌ NO |
| `WD4/WD8/HPC archives` | Frozen trajectories, logs, topologies | ❌ NO |

---

## Rebuild guarantee

At any point, you can:

1. Clone this repo
2. Point it at WD4/WD8/HPC (where frozen evidence lives)
3. Run `stageC_rebuild_views.py`
4. Get back the full symlink forest under `views/runs/`

No data is lost. The **control plane** (this repo) + **freezer** (external storage) = complete lab state.
