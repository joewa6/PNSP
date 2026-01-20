# Integration Tests

## Purpose
End-to-end tests validating complete pipelines (input → Stage A → Stage B → ...).

---

## stageA_batch_run_20012026

**Test Date:** 2026-01-20  
**Git Commit:** cc8b307  
**Environment:** Cheminf conda  

### Test Questions

1. **Does stageA run and return files we want?**
   - ✅ YES: All 5 molecules processed successfully
   - Output: `system.prmtop` + `system.inpcrd` per molecule
   - Files contain full solvated NP-ligand complex systems (~16MB topology + ~3MB coords)
   - Each run generates metadata in `run_manifest.json`

2. **Does it catch wrong environment?**
   - ✅ YES: With base conda (missing rdkit):
     - Previous test: 5/5 failed with clear error: `"Current profile does not contain rdkit,epik"`
     - No silent fallback to incomplete systems
     - Batch processing reported all failures in `batch_summary.json`
   
### Test Data

**Input:** `/home/jwallace-iit.local/WORK/MSE_LAB/Inputs/stageA/ligands/example_ligands.csv`
- 5 molecules: Ethanol, Diuron, Aspirin, Caffeine, Benzene
- Various molecular complexity (3-24 atoms)

**Output Structure:**
```
tests/integration/stageA_batch_run_20012026/
├── Aspirin/           → system.prmtop, system.inpcrd, run_manifest.json, intermediates
├── Benzene/           → system.prmtop, system.inpcrd, run_manifest.json, intermediates
├── Caffeine/          → system.prmtop, system.inpcrd, run_manifest.json, intermediates
├── Diuron/            → system.prmtop, system.inpcrd, run_manifest.json, intermediates
├── Ethanol/           → system.prmtop, system.inpcrd, run_manifest.json, intermediates
└── batch_summary.json → Batch report (5 successful, 0 failed)
```

### Key Validations Passed

✅ Epik microstate enumeration works  
✅ Ligand parameterization (AM1-BCC) works  
✅ NP complex assembly with 10 ligand copies works  
✅ Full solvation with TIP3P + counter ions works  
✅ Environment checking catches missing dependencies  
✅ Batch processing completes all molecules  
✅ Output files have correct sizes (16MB prmtop, 3MB inpcrd per system)  
✅ Batch summary JSON correctly reports all 5 successes  

### Next Steps

1. **Detailed error catching** (not yet tested):
   - Wrong molecule format (e.g., invalid SMILES)
   - Invalid charge specifications
   - Isomer enumeration edge cases
   - File I/O errors
   - Memory constraints

2. **Stage B integration**: Test that these outputs feed into simulation pipeline

