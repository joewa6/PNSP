# Representative solvated systems

`phenol/` and `ammonium/` contain the first and last saved frames from the
canonical 15 ns production trajectory of replica 1. Each contains ten
contaminant copies and retains all PET, TIP3P water and counterion atoms.
Ammonium is labelled `Ammonia` in the historical archive; its recorded input
SMILES is not a claim about the final parameterized protonation state.

Load each `frame_*.nc` with the matching `system.prmtop`. These are individual
snapshots, not consecutive samples for residence-time analysis. Frame indices,
stored times, production-relative times, box dimensions, composition and SHA-256 checksums are recorded
in `manifest.json`. Original build timestamps and Git revisions are distinguished
from snapshot-export timestamps.

`MDIN/` contains the original simulation inputs; these expect the historical
restart filenames and are provided as protocol records, not as a standalone
restart chain. Use the main workflow to prepare a new simulation. The ligand
library, force-field parameters and structures accompany each example; PET
library and parameter files are in `examples/example_pet/`.

To reproduce an export from an available archive:

```bash
python examples/export_snapshot.py --archive-run /path/to/Phenol/rep1 \
  --system Phenol --output /path/to/new_snapshot_directory
```
