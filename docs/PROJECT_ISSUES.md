# MSE_LAB Construction Blueprint

**Executable project board specification for 3-day lab bootstrap**

Generated: 2026-01-13  
Repository: `joewa6/MSE_LAB`  
Branch strategy: Work on `dev`, merge to `main` only when stages actually work

---

## 📊 PROJECT BOARD STRUCTURE

### Columns:

1. **Backlog** — Not started, blocked, or future work
2. **Day 1** — Active work for Day 1 milestone
3. **Day 2** — Active work for Day 2 milestone  
4. **Day 3** — Active work for Day 3 milestone
5. **Review** — Awaiting validation/testing
6. **Done** — Acceptance criteria met, committed to `dev`

---

## 🎯 MILESTONES

### DAY 1: Wrap A & B + Prove One End-to-End Run
**Goal**: One ligand runs end-to-end locally and appears in Stage C view

### DAY 2: Full Spine Integration
**Goal**: Multi-run support + legacy WD8 indexing + archive pointers

### DAY 3: Reduce + HPC Skeleton + Hardening
**Goal**: Metrics extraction + HPC wrapper + production docs

---

## 📋 ISSUES (Execution Order)

---

## DAY 1 ISSUES

### Issue #1: Port Stage A builder into Orchestra

**Column**: Day 1  
**Labels**: `capability`, `stageA`, `day1`, `critical`  
**Dependencies**: None (blocking all other work)

**Goal**: Transplant legacy builder into `Orchestra/stageA_builder/src/` with clean entrypoint

**Required files/dirs**:
```
Orchestra/stageA_builder/
  src/
    build.py              # main entrypoint
    templates/            # force field / topology templates
  README.md               # usage + requirements
```

**Acceptance criteria**:
- [ ] `python Orchestra/stageA_builder/src/build.py --smiles "..." --output ./test_build/` succeeds
- [ ] Emits `run_manifest.json` with schema:
  ```json
  {
    "project": "PET_SCREENING",
    "ligand_norm": "...",
    "smiles": "...",
    "stageA": {
      "status": "complete",
      "builder_version": "...",
      "git_commit": "...",
      "timestamp": "ISO-8601"
    },
    "paths": {
      "topology": "relative/path/to/.prmtop",
      "coords": "relative/path/to/.inpcrd"
    }
  }
  ```
- [ ] Generates `.prmtop` + `.inpcrd` (or equivalent for GROMACS/LAMMPS)
- [ ] No hardcoded absolute paths
- [ ] Works with Python 3.11+ in fresh conda env

**Done when**:
```bash
cd Orchestra/stageA_builder
python src/build.py --smiles "CC(C)Nc1nc(Cl)nc(NC(C)C)n1" --output /tmp/test_diuron
ls /tmp/test_diuron/run_manifest.json
cat /tmp/test_diuron/run_manifest.json | jq .stageA.status
# outputs: "complete"
```

---

### Issue #2: Port Stage B local simulator into Orchestra

**Column**: Day 1  
**Labels**: `capability`, `stageB`, `day1`, `critical`  
**Dependencies**: #1 (needs topology files from Stage A)

**Goal**: Wrap legacy MD runner into `Orchestra/stageB_simulator/local/src/` with manifest-driven execution

**Required files/dirs**:
```
Orchestra/stageB_simulator/
  local/
    src/
      run.py            # main entrypoint
      settings/
        local_tiny_v1.yaml    # short test config (100 steps)
    README.md
  templates/
    amber.in.j2         # optional: Jinja2 template for .in files
```

**Acceptance criteria**:
- [ ] `python Orchestra/stageB_simulator/local/src/run.py --manifest ./run_manifest.json --settings local_tiny_v1` succeeds
- [ ] Reads topology paths from manifest `paths.topology` and `paths.coords`
- [ ] Updates `run_manifest.json` with:
  ```json
  {
    "stageB": {
      "status": "complete",
      "engine": "amber22",
      "settings_id": "local_tiny_v1",
      "started_at": "...",
      "completed_at": "...",
      "total_steps": 100
    },
    "paths": {
      "trajectory": "relative/path/to/prod.nc",
      "log": "relative/path/to/prod.out"
    }
  }
  ```
- [ ] Generates trajectory file (`.nc` or `.xtc`)
- [ ] Does NOT commit trajectories to git (`.gitignore` blocks them)
- [ ] Logs to stdout + optional file

**Done when**:
```bash
# After Stage A completes
cd /tmp/test_diuron
python ~/WORK/MSE_LAB/Orchestra/stageB_simulator/local/src/run.py \
  --manifest run_manifest.json \
  --settings local_tiny_v1

ls prod.nc prod.out
cat run_manifest.json | jq .stageB.status
# outputs: "complete"
```

---

### Issue #3: Create Stage C pointer-based archival

**Column**: Day 1  
**Labels**: `capability`, `stageC`, `day1`, `critical`  
**Dependencies**: #2 (needs completed runs)

**Goal**: Archive completed runs to local toy archive + create pointer views

**Required files/dirs**:
```
Orchestra/stageC_archive/
  src/
    archive.py          # freezes run to archive root
    rebuild_views.py    # rebuilds views from pointers
  pointers/             # .ptr files (JSON with archive locations)
  views/
    runs/               # symlink forest (gitignored)
```

**Acceptance criteria**:
- [ ] `python Orchestra/stageC_archive/src/archive.py --manifest /path/to/run_manifest.json --archive-root /tmp/toy_archive` succeeds
- [ ] Copies entire run directory to archive root with structure:
  ```
  /tmp/toy_archive/
    runs/
      2026-01-13/
        {run_id}/
          run_manifest.json
          topology.prmtop
          prod.nc
          prod.out
  ```
- [ ] Creates pointer file at `Orchestra/stageC_archive/pointers/{run_id}.ptr`:
  ```json
  {
    "run_id": "...",
    "archive_root": "/tmp/toy_archive",
    "archive_path": "runs/2026-01-13/{run_id}",
    "archived_at": "ISO-8601",
    "origin": "v2_local"
  }
  ```
- [ ] Updates manifest with `stageC.status = "archived"`
- [ ] `python Orchestra/stageC_archive/src/rebuild_views.py` creates symlinks:
  ```
  Orchestra/stageC_archive/views/runs/{run_id}/
    run_manifest.json → /tmp/toy_archive/runs/2026-01-13/{run_id}/run_manifest.json
    prod.nc → /tmp/toy_archive/runs/2026-01-13/{run_id}/prod.nc
  ```
- [ ] `views/runs/` is gitignored (no symlinks tracked)

**Done when**:
```bash
python Orchestra/stageC_archive/src/archive.py \
  --manifest /tmp/test_diuron/run_manifest.json \
  --archive-root /tmp/toy_archive

ls Orchestra/stageC_archive/pointers/*.ptr
python Orchestra/stageC_archive/src/rebuild_views.py
ls -l Orchestra/stageC_archive/views/runs/*/run_manifest.json
# shows symlink to /tmp/toy_archive/...
```

---

### Issue #4: End-to-end smoke test (A → B → C)

**Column**: Day 1  
**Labels**: `validation`, `day1`, `critical`  
**Dependencies**: #1, #2, #3 (all stages must exist)

**Goal**: Prove one ligand can flow through all stages without manual intervention

**Required files/dirs**:
```
scripts/
  smoke_day1.sh       # automated end-to-end test
```

**Acceptance criteria**:
- [ ] Script runs Stage A → B → C for one SMILES
- [ ] Generates unique `run_id` (timestamp + ligand hash)
- [ ] Archives to `/tmp/toy_archive`
- [ ] Rebuilds view
- [ ] Verifies symlink exists and points to correct archive location
- [ ] Exits 0 on success, non-zero on failure
- [ ] No manual steps required

**Done when**:
```bash
bash scripts/smoke_day1.sh
echo $?
# outputs: 0

ls Orchestra/stageC_archive/views/runs/
# shows one run directory with symlinks
```

---

### Issue #5: Create `run_manifest.json` schema definition

**Column**: Day 1  
**Labels**: `docs`, `day1`  
**Dependencies**: #1, #2, #3 (schema extracted from working stages)

**Goal**: Formalize the manifest contract so stages can be built independently

**Required files/dirs**:
```
docs/
  run_manifest_schema.yaml    # JSON Schema or Pydantic-style spec
```

**Acceptance criteria**:
- [ ] Documents all required fields: `project`, `ligand_norm`, `smiles`, `stageA`, `stageB`, `stageC`, `paths`
- [ ] Specifies types, enums, constraints
- [ ] Includes example manifest (valid JSON)
- [ ] States versioning policy (manifest schema is immutable, add new fields, never delete)

**Done when**:
```bash
cat docs/run_manifest_schema.yaml
# shows complete schema with descriptions
```

---

## DAY 2 ISSUES

### Issue #6: Stage A batch mode (read `test.csv`)

**Column**: Day 2  
**Labels**: `capability`, `stageA`, `day2`  
**Dependencies**: #1 (extends Stage A)

**Goal**: Build multiple ligands from CSV in one invocation

**Required files/dirs**:
```
test.csv                # example: ligand_norm, smiles, charge
Orchestra/stageA_builder/src/
  batch_build.py        # reads CSV, calls build.py per row
```

**Acceptance criteria**:
- [ ] `python Orchestra/stageA_builder/src/batch_build.py --csv test.csv --output ./builds/` succeeds
- [ ] Creates one subdirectory per ligand under `./builds/{ligand_norm}/`
- [ ] Each contains `run_manifest.json` + topology files
- [ ] Logs progress (N/M complete)
- [ ] Handles errors gracefully (logs failure, continues)

**Done when**:
```bash
python Orchestra/stageA_builder/src/batch_build.py --csv test.csv --output /tmp/batch_test
ls /tmp/batch_test/*/run_manifest.json | wc -l
# outputs: (number of rows in test.csv)
```

---

### Issue #7: Stage C legacy WD8 indexer

**Column**: Day 2  
**Labels**: `capability`, `stageC`, `day2`, `critical`  
**Dependencies**: #3 (extends Stage C)

**Goal**: Index existing WD8 archives and create pointer files for legacy runs

**Required files/dirs**:
```
Orchestra/stageC_archive/src/
  index_legacy.py       # scans WD8, emits .ptr files
```

**Acceptance criteria**:
- [ ] `python Orchestra/stageC_archive/src/index_legacy.py --root /mnt/WD8/PET_SCREENING --origin legacy_wd8` succeeds
- [ ] Discovers all `run_manifest.json` files under WD8 tree
- [ ] Creates `.ptr` files in `pointers/` with:
  ```json
  {
    "run_id": "...",
    "archive_root": "/mnt/WD8/PET_SCREENING",
    "archive_path": "relative/path/to/run",
    "origin": "legacy_wd8",
    "indexed_at": "ISO-8601"
  }
  ```
- [ ] Handles missing manifests gracefully (logs warning, skips)
- [ ] Does NOT copy data (only creates pointers)
- [ ] `rebuild_views.py` includes legacy runs in symlink forest

**Done when**:
```bash
python Orchestra/stageC_archive/src/index_legacy.py \
  --root /mnt/WD8/PET_SCREENING \
  --origin legacy_wd8

ls Orchestra/stageC_archive/pointers/*.ptr | wc -l
# shows total count (legacy + v2)

python Orchestra/stageC_archive/src/rebuild_views.py
ls Orchestra/stageC_archive/views/runs/ | head
# shows mix of legacy and v2 runs
```

---

### Issue #8: Stage B rsync archival after simulation

**Column**: Day 2  
**Labels**: `capability`, `stageB`, `day2`  
**Dependencies**: #2, #3 (Stage B must trigger Stage C)

**Goal**: Automatically archive completed runs to designated root after Stage B finishes

**Required files/dirs**:
```
Orchestra/stageB_simulator/local/src/
  run.py (modified)     # calls archive after completion
```

**Acceptance criteria**:
- [ ] Stage B accepts `--archive-root` flag
- [ ] On successful completion, invokes `stageC_archive/src/archive.py`
- [ ] Run manifest updated with `stageC` metadata
- [ ] Pointer file created automatically
- [ ] User can optionally skip archival with `--no-archive`

**Done when**:
```bash
python Orchestra/stageB_simulator/local/src/run.py \
  --manifest /tmp/test_diuron/run_manifest.json \
  --settings local_tiny_v1 \
  --archive-root /tmp/toy_archive

ls Orchestra/stageC_archive/pointers/*.ptr
cat /tmp/test_diuron/run_manifest.json | jq .stageC.status
# outputs: "archived"
```

---

### Issue #9: Composer ledger scan (multi-archive support)

**Column**: Day 2  
**Labels**: `capability`, `composer`, `day2`, `critical`  
**Dependencies**: #7 (needs pointer files for legacy + v2)

**Goal**: Composer can read all pointer files and build unified run ledger

**Required files/dirs**:
```
Composer/
  state/
    ledger.py           # reads pointers, builds in-memory table
  cli/
    list_runs.py        # CLI to query ledger
```

**Acceptance criteria**:
- [ ] `python Composer/cli/list_runs.py` succeeds
- [ ] Scans all `.ptr` files in `Orchestra/stageC_archive/pointers/`
- [ ] Outputs table with columns: `run_id`, `ligand_norm`, `origin`, `status`, `archive_path`
- [ ] Distinguishes `legacy_wd8` vs `v2_local` in `origin` column
- [ ] Supports filtering: `--origin legacy_wd8` or `--status archived`
- [ ] Reads manifests on-demand (via symlink view)

**Done when**:
```bash
python Composer/cli/list_runs.py
# outputs table showing all runs (legacy + new)

python Composer/cli/list_runs.py --origin v2_local | wc -l
# shows count of new runs only
```

---

### Issue #10: Composer `origin` column logic

**Column**: Day 2  
**Labels**: `capability`, `composer`, `day2`  
**Dependencies**: #9 (extends ledger)

**Goal**: Ledger tracks provenance (legacy vs v2) for all runs

**Acceptance criteria**:
- [ ] Pointer files include `"origin"` field
- [ ] Ledger exposes `origin` in queries
- [ ] `origin` values: `legacy_wd8`, `v2_local`, `v2_hpc` (future)
- [ ] Documented in `docs/STAGES.md` under Stage C

**Done when**:
```bash
python Composer/cli/list_runs.py --format json | jq '.[].origin' | sort | uniq
# outputs:
# "legacy_wd8"
# "v2_local"
```

---

## DAY 3 ISSUES

### Issue #11: Stage D metric extraction (first reducer)

**Column**: Day 3  
**Labels**: `capability`, `stageD`, `day3`, `critical`  
**Dependencies**: #7, #9 (needs unified view of archived runs)

**Goal**: Extract one metric (e.g., final energy) from archived trajectories

**Required files/dirs**:
```
Orchestra/stageD_reduce/
  src/
    extract_metrics.py    # reads from views, emits metrics
  outputs/
    run_stats.json        # or .parquet
```

**Acceptance criteria**:
- [ ] `python Orchestra/stageD_reduce/src/extract_metrics.py --runs-dir Orchestra/stageC_archive/views/runs/` succeeds
- [ ] Reads `prod.out` or `prod.nc` from symlink view
- [ ] Extracts one metric (final energy, RMSD, or total steps)
- [ ] Writes to `outputs/run_stats.json`:
  ```json
  [
    {
      "run_id": "...",
      "ligand_norm": "...",
      "metric_name": "final_energy",
      "value": -12345.67,
      "unit": "kcal/mol",
      "computed_at": "ISO-8601"
    }
  ]
  ```
- [ ] Handles missing files gracefully (logs warning, skips)
- [ ] Does NOT modify archived data (read-only)

**Done when**:
```bash
python Orchestra/stageD_reduce/src/extract_metrics.py \
  --runs-dir Orchestra/stageC_archive/views/runs/

cat Orchestra/stageD_reduce/outputs/run_stats.json | jq length
# shows number of runs processed
```

---

### Issue #12: Composer rank table from Stage D metrics

**Column**: Day 3  
**Labels**: `capability`, `composer`, `day3`  
**Dependencies**: #11 (needs metrics)

**Goal**: Composer CLI can display ranked list of runs by metric

**Required files/dirs**:
```
Composer/cli/
  rank_runs.py          # joins ledger + metrics, sorts
```

**Acceptance criteria**:
- [ ] `python Composer/cli/rank_runs.py --metric final_energy` succeeds
- [ ] Reads ledger (from pointers)
- [ ] Reads metrics (from `Orchestra/stageD_reduce/outputs/run_stats.json`)
- [ ] Joins on `run_id`
- [ ] Outputs sorted table: `rank`, `run_id`, `ligand_norm`, `final_energy`, `origin`
- [ ] Supports ascending/descending order

**Done when**:
```bash
python Composer/cli/rank_runs.py --metric final_energy --order ascending
# outputs ranked table with best (lowest) energy first
```

---

### Issue #13: Stage B HPC skeleton

**Column**: Day 3  
**Labels**: `capability`, `stageB`, `day3`  
**Dependencies**: #2 (parallel to local runner)

**Goal**: Create HPC submission wrapper (does not need to run, just structure)

**Required files/dirs**:
```
Orchestra/stageB_simulator/
  hpc/
    src/
      submit.py         # generates SLURM script, submits job
      templates/
        slurm.sh.j2     # Jinja2 template
    README.md           # usage + cluster-specific notes
```

**Acceptance criteria**:
- [ ] `python Orchestra/stageB_simulator/hpc/src/submit.py --manifest ./run_manifest.json --settings hpc_prod_v1 --dry-run` succeeds
- [ ] Generates SLURM script with:
  - Job name (from `run_id`)
  - Resource requests (nodes, time, partition)
  - Module loads (Amber, Python)
  - Execution command
- [ ] Does NOT actually submit in `--dry-run` mode
- [ ] Documents required environment (module names, paths)
- [ ] Placeholder for rsync-back logic (commented out)

**Done when**:
```bash
python Orchestra/stageB_simulator/hpc/src/submit.py \
  --manifest /tmp/test/run_manifest.json \
  --settings hpc_prod_v1 \
  --dry-run

cat slurm_*.sh
# shows valid SLURM script with placeholders
```

---

### Issue #14: Repository README with quickstart

**Column**: Day 3  
**Labels**: `docs`, `day3`  
**Dependencies**: All critical path issues (proves system works)

**Goal**: Document how to use MSE_LAB from clone to first run

**Required files/dirs**:
```
README.md (root)
```

**Acceptance criteria**:
- [ ] Describes system purpose (no fluff)
- [ ] Quickstart:
  ```bash
  git clone git@github.com:joewa6/MSE_LAB.git
  cd MSE_LAB
  conda env create -f environment.yml
  conda activate mse_lab
  bash scripts/smoke_day1.sh
  ```
- [ ] Links to `docs/ARCHITECTURE.md`, `docs/STAGES.md`, `docs/EVIDENCE_POLICY.md`
- [ ] States branch policy (`main` = stable, `dev` = work)
- [ ] States what is NOT in repo (WD4/WD8 data)

**Done when**:
```bash
cat README.md | head -20
# shows concise purpose + quickstart
```

---

### Issue #15: Pin Python environment (`environment.yml`)

**Column**: Day 3  
**Labels**: `infra`, `day3`  
**Dependencies**: All stage code (to know dependencies)

**Goal**: Reproducible conda environment for all stages

**Required files/dirs**:
```
environment.yml
```

**Acceptance criteria**:
- [ ] Lists all required packages:
  - Python 3.11
  - AmberTools (if needed)
  - MDAnalysis / mdtraj (if used in Stage D)
  - Jinja2 (for templates)
  - PyYAML, pandas, numpy
- [ ] Specifies channels (conda-forge, bioconda)
- [ ] Tested on clean conda install:
  ```bash
  conda env create -f environment.yml
  conda activate mse_lab
  bash scripts/smoke_day1.sh
  # succeeds
  ```

**Done when**:
```bash
conda env create -f environment.yml -n mse_test
conda activate mse_test
python -c "import yaml, pandas, mdtraj"
echo $?
# outputs: 0
```

---

### Issue #16: GitHub branch protection rules

**Column**: Day 3  
**Labels**: `infra`, `day3`  
**Dependencies**: None (can be done anytime)

**Goal**: Protect `main` from accidental commits

**Acceptance criteria**:
- [ ] On GitHub repo settings:
  - Protect `main` branch
  - Require pull request before merging
  - Require status checks: `ci`, `no-evidence`
  - No direct pushes to `main`
- [ ] Documented in `docs/GITHUB_WORKFLOW.md`

**Done when**:
- [ ] Attempt to `git push origin main` fails with protection message
- [ ] Merge from `dev` → `main` requires PR + green CI

---

## 🏷️ LABEL TAXONOMY

| Label | Meaning |
|-------|---------|
| `capability` | Adds new lab capability (always commit-worthy) |
| `validation` | Proves existing capabilities work |
| `docs` | Documentation only |
| `infra` | GitHub/CI/environment setup |
| `stageA` | Builder-related |
| `stageB` | Simulator-related |
| `stageC` | Archive-related |
| `stageD` | Reducer-related |
| `composer` | Ledger/CLI-related |
| `day1` | Must complete by end of Day 1 |
| `day2` | Must complete by end of Day 2 |
| `day3` | Must complete by end of Day 3 |
| `critical` | Blocking other work |

---

## 🔗 DEPENDENCY GRAPH

```
Day 1:
  #1 (Stage A port) → #2 (Stage B port) → #3 (Stage C archive) → #4 (smoke test)
                   ↘ #5 (manifest schema) ↗

Day 2:
  #1 → #6 (batch builds)
  #3 → #7 (legacy indexer) → #9 (Composer ledger) → #10 (origin column)
  #2 + #3 → #8 (B→C archival)

Day 3:
  #7 + #9 → #11 (Stage D metrics) → #12 (Composer rank)
  #2 → #13 (HPC skeleton)
  All → #14 (README)
  All → #15 (environment.yml)
  (independent) → #16 (branch protection)
```

---

## 📦 KANBAN INITIAL PLACEMENT

| Column | Issues |
|--------|--------|
| **Day 1** | #1, #2, #3, #4, #5 |
| **Day 2** | #6, #7, #8, #9, #10 |
| **Day 3** | #11, #12, #13, #14, #15, #16 |
| **Review** | (empty initially) |
| **Done** | (empty initially) |

---

## ✅ COMPLETION CRITERIA (ALL ISSUES DONE)

When all issues are in **Done** column:

```bash
# Full pipeline test
echo "ligand_norm,smiles,charge" > final_test.csv
echo "diuron,CC(C)Nc1nc(Cl)nc(NC(C)C)n1,-50" >> final_test.csv

python Orchestra/stageA_builder/src/batch_build.py --csv final_test.csv --output /tmp/final_builds
cd /tmp/final_builds/diuron
python ~/WORK/MSE_LAB/Orchestra/stageB_simulator/local/src/run.py \
  --manifest run_manifest.json \
  --settings local_tiny_v1 \
  --archive-root /tmp/final_archive

python ~/WORK/MSE_LAB/Orchestra/stageC_archive/src/rebuild_views.py
python ~/WORK/MSE_LAB/Orchestra/stageD_reduce/src/extract_metrics.py \
  --runs-dir ~/WORK/MSE_LAB/Orchestra/stageC_archive/views/runs/

python ~/WORK/MSE_LAB/Composer/cli/rank_runs.py --metric final_energy
# Shows ranked table with diuron + any legacy runs
```

If above succeeds → **MSE_LAB is operational**

---

## 🚀 NEXT STEPS (AFTER BLUEPRINT COMMIT)

1. Commit this file to `dev`:
   ```bash
   git add docs/PROJECT_ISSUES.md
   git commit -m "Blueprint: executable project board for 3-day bootstrap"
   git push
   ```

2. Create GitHub issues (manually or via script):
   - Use `.github/ISSUE_TEMPLATE/stage_capability.yml` for capability issues
   - Copy acceptance criteria from this document
   - Apply labels as specified

3. Create GitHub Project Board:
   - Board name: "MSE_LAB Bootstrap"
   - Add columns: Backlog, Day 1, Day 2, Day 3, Review, Done
   - Assign issues to columns per table above

4. Start execution:
   - Move #1 to "Day 1" column
   - Begin porting Stage A
   - Mark complete when acceptance criteria met
   - Commit to `dev` with message: "Stage A: [brief description]"

---

**This document is the executable specification.**  
**Every issue is a forcing function.**  
**The lab exists when all issues are Done.**
