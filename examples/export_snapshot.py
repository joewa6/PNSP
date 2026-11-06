#!/usr/bin/env python3
"""Export the first and last production frames, with all solvent and ions."""
import argparse
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

import MDAnalysis as mda
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive-run", type=Path, required=True)
    parser.add_argument("--system", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source, output = args.archive_run, args.output
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite {output}")
    historical = json.loads((source / "run_manifest.json").read_text())
    universe = mda.Universe(str(source / "system.prmtop"), str(source / "prod.nc"))
    if len(universe.trajectory) != 1500:
        raise ValueError("Expected 1,500 frames from the 15 ns production trajectory")
    if not np.isclose(universe.trajectory.dt, 10.):
        raise ValueError("Expected 10 ps production frame spacing")
    counts = {name: int(sum(universe.residues.resnames == name))
              for name in np.unique(universe.residues.resnames)}
    if counts.get("LIG") != 10 or counts.get("WAT", 0) == 0:
        raise ValueError(f"Unexpected ligand/solvent composition: {counts}")
    output.mkdir(parents=True)
    for name in ("system.prmtop", "LIG.lib", "LIG.frcmod", "LIG.mol2", "LIG.pdb"):
        shutil.copy2(source / name, output / name)
    shutil.copytree(source / "MDIN", output / "MDIN")
    frames = []
    for index in (0, 1499):
        ts = universe.trajectory[index]
        if not np.isfinite(ts.positions).all() or ts.dimensions is None:
            raise ValueError("Invalid snapshot coordinates or box")
        name = f"frame_{index:04d}.nc"
        with mda.Writer(str(output / name), n_atoms=len(universe.atoms)) as writer:
            writer.write(universe.atoms)
        check = mda.Universe(str(output / "system.prmtop"), str(output / name))
        np.testing.assert_allclose(check.atoms.positions, ts.positions, atol=1e-4)
        np.testing.assert_allclose(check.dimensions, ts.dimensions, atol=1e-4)
        frames.append(dict(file=name, source_frame_index=index, time_ps=float(ts.time),
                           production_elapsed_ps=(index + 1) * 10.,
                           box_A_degrees=ts.dimensions.tolist()))
        check.trajectory.close()
    universe.trajectory.close()
    hashes = {str(p.relative_to(output)): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in sorted(output.rglob("*")) if p.is_file()}
    manifest = dict(
        record_type="archival_snapshot_export",
        exported_at_utc=datetime.now(timezone.utc).isoformat(),
        contaminant_id=args.system,
        replica_id=f"rep{historical['replica']['index']}",
        historical_build=historical["stageA"],
        historical_contaminant_id=historical["ligand_norm"],
        historical_input_smiles=historical["smiles"],
        packmol_seed=historical["replica"].get("packmol_seed"),
        source_archive_id="NEW_203_15NS",
        source_run=f"{args.system}/rep{historical['replica']['index']}",
        source_trajectory="prod.nc", source_frames=1500,
        atoms=len(universe.atoms), residue_counts=counts,
        topology_net_charge_e=float(universe.atoms.charges.sum()),
        snapshots=frames, sha256=hashes,
        note="Historical build metadata is preserved from run_manifest.json. Export time is not simulation time. Snapshots retain every atom and are not a continuous trajectory.",
    )
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Validated {args.system}: {len(universe.atoms)} atoms, {counts}")


if __name__ == "__main__":
    main()
