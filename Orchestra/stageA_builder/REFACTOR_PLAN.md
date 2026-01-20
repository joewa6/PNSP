# Stage A Refactoring Plan: Inputs/Runs Architecture

## Changes Required

### 1. Update parse_args() - Add New CLI Arguments

**New arguments:**
```python
--project          # Project name (e.g., "PET_SCREENING")
--run-root         # Path to Runs/ project directory  
--np-profile       # Path to NP YAML config (e.g., Inputs/stageA/nanoparticle/NP_profile_c50.yaml)
```

**Keep existing:**
```python
--smiles / --csv   # Input mode (mutually exclusive)
--dry-run          # Testing flag
```

**Remove:**
```python
--output           # Replaced by --run-root + auto-generated run_id
--output-dir       # Replaced by --run-root + auto-generated run_id
--templates        # Not needed anymore (templates in stage code directory)
```

### 2. Create NP Profile Loader

**Function:** `load_np_profile(yaml_path: Path) -> dict`

Load and validate NP profile YAML:
```yaml
profile_id: c50
description: -50e charged gold nanoparticle
charge: -50
files:
  lib: NP_c50/NP_correct.lib
  charged_frcmod: NP_c50/NP_charged.frcmod
  neut_frcmod: NP_c50/NP_neut.frcmod
  dry_pdb: NP_c50/NP_dry.pdb
```

Returns resolved absolute paths.

### 3. Generate run_id

**Function:** `generate_run_id(project: str, np_profile_id: str, ligand_set: str, git_commit: str) -> str`

Format:
```
<timestamp>__<charge_state>__<ligand_set>__git-<commit>

Example:
2026-01-14T16-05-12Z__c50__example_ligands__git-a1b2c3d
```

### 4. Update _build_np_complex()

**Change signature:**
```python
def _build_np_complex(
    out_dir: pathlib.Path,
    np_files: dict[str, pathlib.Path],  # NEW: from NP profile
    ligand_copies: int = 10
) -> None:
```

**Changes:**
- Remove hardcoded NP directory search
- Use `np_files["lib"]`, `np_files["dry_pdb"]`, etc.
- Use absolute paths (no relative path tricks)

### 5. Update run_builder()

**Add parameter:** `np_files: dict[str, pathlib.Path] | None = None`

Pass to `_build_np_complex()`.

### 6. Update process_csv_batch()

**Change signature:**
```python
def process_csv_batch(
    csv_path: pathlib.Path,
    run_dir: pathlib.Path,        # NEW: full Runs/<project>/<run_id>/
    np_files: dict[str, pathlib.Path],  # NEW: from NP profile
    git_commit: str,
    dry_run: bool = False,
) -> dict[str, Any]:
```

**Changes:**
- Write outputs to `run_dir / "stageA" / "systems" / <ligand_id>/`
- Update manifest structure to match Runs/README.md spec
- Add `inputs_snapshot/` copy of CSV and NP profile

### 7. Update main()

**New flow:**

```python
def main(argv: list[str]) -> int:
    args = parse_args(argv)
    
    # Load NP profile
    np_profile_path = Path(args.np_profile).resolve()
    np_profile = load_np_profile(np_profile_path)
    np_files = resolve_np_files(np_profile_path.parent, np_profile)
    
    # Get git commit
    repo_root = find_repo_root(Path(__file__).parent)
    git_commit = get_git_commit(repo_root)
    
    # Generate run_id
    if args.csv:
        ligand_set = Path(args.csv).stem  # e.g., "example_ligands"
    else:
        ligand_set = "single"
    
    run_id = generate_run_id(
        project=args.project,
        np_profile_id=np_profile["profile_id"],
        ligand_set=ligand_set,
        git_commit=git_commit
    )
    
    # Create run directory
    run_dir = Path(args.run_root) / run_id
    safe_mkdir(run_dir)
    safe_mkdir(run_dir / "inputs_snapshot")
    safe_mkdir(run_dir / "stageA" / "systems")
    
    # Snapshot inputs
    shutil.copy(np_profile_path, run_dir / "inputs_snapshot" / "np_profile.yaml")
    if args.csv:
        shutil.copy(args.csv, run_dir / "inputs_snapshot" / "ligands.csv")
    
    # Initialize run_manifest.json
    manifest = {
        "project": args.project,
        "run_id": run_id,
        "created": datetime.utcnow().isoformat() + "Z",
        "git_commit": git_commit,
        "inputs": {
            "np_profile": "inputs_snapshot/np_profile.yaml"
        },
        "stageA": {"status": "in_progress"},
        "stageB": {"status": "pending"},
        "stageC": {"status": "pending"},
        "stageD": {"status": "pending"}
    }
    
    if args.csv:
        manifest["inputs"]["ligands_csv"] = "inputs_snapshot/ligands.csv"
        summary = process_csv_batch(
            csv_path=Path(args.csv),
            run_dir=run_dir,
            np_files=np_files,
            git_commit=git_commit,
            dry_run=args.dry_run
        )
        # Update manifest with results
        manifest["stageA"]["status"] = "complete"
        manifest["stageA"]["systems"] = summary["systems"]
        
    else:  # Single SMILES mode
        # Similar but for single molecule
        pass
    
    # Write manifest
    manifest_path = run_dir / "run_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2))
    
    return 0
```

### 8. Update Documentation

**Files to update:**
- `Orchestra/stageA_builder/README.md` - New CLI usage examples
- `Orchestra/stageA_builder/CHANGELOG.md` - Document architectural change

## Testing Plan

1. Test with example CSV:
   ```bash
   python Orchestra/stageA_builder/src/build.py \
     --csv Inputs/stageA/ligands/example_ligands.csv \
     --np-profile Inputs/stageA/nanoparticle/NP_profile_c50.yaml \
     --project TEST \
     --run-root Runs/__SMOKE__/
   ```

2. Verify Runs/ structure matches spec

3. Verify manifest format

4. Test single SMILES mode

## Implementation Order

1. ✅ Create Inputs/ structure
2. ✅ Create Runs/ structure  
3. ✅ Create NP profile YAML
4. ✅ Move canonical NP files
5. ✅ Update .gitignore
6. 🔄 Refactor build.py (NEXT)
   - Add load_np_profile()
   - Add generate_run_id()
   - Update parse_args()
   - Update _build_np_complex()
   - Update run_builder()
   - Update process_csv_batch()
   - Update main()
7. Update README.md
8. Update CHANGELOG.md
9. Test with smoke test
10. Commit to git
