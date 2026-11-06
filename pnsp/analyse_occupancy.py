#!/usr/bin/env python3
"""Calculate copy-level PET contact occupancy and contiguous contact events."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

if __package__:
    from .provenance import execution_provenance
else:
    from provenance import execution_provenance

import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.lib.distances import capped_distance


DEFAULT_PET_RESNAMES = ("F01", "F02", "F03", "F04", "F05", "F06")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def heavy_indices(atoms: mda.core.groups.AtomGroup) -> np.ndarray:
    try:
        masses = np.asarray(atoms.masses, dtype=float)
        selected = atoms.indices[masses > 1.5]
    except (AttributeError, mda.exceptions.NoDataError):
        selected = atoms.indices[
            np.array([not name.upper().startswith("H") for name in atoms.names])
        ]
    if len(selected) == 0:
        raise RuntimeError("Heavy-atom selection is empty")
    return selected


def events(mask: np.ndarray) -> list[tuple[int, int]]:
    padded = np.r_[False, mask.astype(bool), False].astype(np.int8)
    changes = np.diff(padded)
    starts = np.flatnonzero(changes == 1)
    ends = np.flatnonzero(changes == -1)
    return list(zip(starts, ends))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", required=True, type=Path)
    parser.add_argument("--trajectory", required=True, type=Path, nargs="+")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--system", default=None)
    parser.add_argument("--replica", default=None)
    parser.add_argument("--ligand-resname", default="LIG")
    parser.add_argument("--pet-resnames", nargs="+", default=list(DEFAULT_PET_RESNAMES))
    parser.add_argument("--cutoff-A", type=float, default=4.0)
    parser.add_argument(
        "--frame-spacing-ps",
        type=float,
        default=10.0,
        help=(
            "Saved-frame interval in ps (default: 10, as used in the paper). "
            "This is explicit because AMBER NetCDF metadata may report the MD "
            "integration timestep rather than the coordinate-write interval."
        ),
    )
    parser.add_argument("--expected-frames", type=int, default=None)
    parser.add_argument("--write-timeseries", action="store_true")
    args = parser.parse_args()

    topology = args.topology.resolve()
    trajectories = [path.resolve() for path in args.trajectory]
    for path in [topology, *trajectories]:
        if not path.is_file():
            raise FileNotFoundError(path)
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    system = args.system or topology.parent.parent.name
    replica = args.replica or topology.parent.name

    universe = mda.Universe(str(topology), *map(str, trajectories))
    pet_residues = [residue for residue in universe.residues if residue.resname in args.pet_resnames]
    ligands = [residue for residue in universe.residues if residue.resname == args.ligand_resname]
    if not pet_residues or not ligands:
        raise RuntimeError(
            f"Selections are empty: PET residues={len(pet_residues)}, LIG residues={len(ligands)}"
        )
    pet_atoms = universe.atoms[
        heavy_indices(universe.atoms[np.concatenate([res.atoms.indices for res in pet_residues])])
    ]
    ligand_heavy_indices: list[int] = []
    atom_to_copy: list[int] = []
    for copy_index, residue in enumerate(ligands):
        indices = heavy_indices(residue.atoms)
        ligand_heavy_indices.extend(indices.tolist())
        atom_to_copy.extend([copy_index] * len(indices))
    ligand_atoms = universe.atoms[np.asarray(ligand_heavy_indices, dtype=int)]
    atom_to_copy_array = np.asarray(atom_to_copy, dtype=int)

    n_frames = len(universe.trajectory)
    if args.expected_frames is not None and n_frames != args.expected_frames:
        raise RuntimeError(f"Expected {args.expected_frames} frames; trajectory has {n_frames}")
    dt_ps = args.frame_spacing_ps
    if not np.isfinite(dt_ps) or dt_ps <= 0:
        raise ValueError("Frame spacing is unavailable; provide --frame-spacing-ps")

    bound = np.zeros((n_frames, len(ligands)), dtype=bool)
    for frame_index, timestep in enumerate(universe.trajectory):
        dimensions = np.asarray(timestep.dimensions, dtype=float)
        if dimensions.size < 6 or not np.all(np.isfinite(dimensions[:3])) or np.any(dimensions[:3] <= 0):
            raise RuntimeError(f"Frame {frame_index} has no valid periodic box")
        pairs, distances = capped_distance(
            ligand_atoms.positions,
            pet_atoms.positions,
            max_cutoff=args.cutoff_A,
            box=dimensions,
            return_distances=True,
        )
        # capped_distance includes its upper bound; the manuscript criterion is strict.
        pairs = pairs[distances < args.cutoff_A]
        if len(pairs):
            bound[frame_index, np.unique(atom_to_copy_array[pairs[:, 0]])] = True

    copy_rows: list[dict[str, object]] = []
    event_rows: list[dict[str, object]] = []
    for copy_index in range(len(ligands)):
        mask = bound[:, copy_index]
        copy_events = events(mask)
        durations = np.asarray([(end - start) * dt_ps for start, end in copy_events], dtype=float)
        for event_index, (start, end) in enumerate(copy_events, start=1):
            event_rows.append(
                {
                    "system": system,
                    "replica": replica,
                    "copy_id": copy_index + 1,
                    "event_id": event_index,
                    "start_frame": start,
                    "end_frame_inclusive": end - 1,
                    "duration_ns": (end - start) * dt_ps / 1000.0,
                    "left_censored": bool(start == 0),
                    "right_censored": bool(end == n_frames),
                }
            )
        copy_rows.append(
            {
                "system": system,
                "replica": replica,
                "copy_id": copy_index + 1,
                "n_frames": n_frames,
                "frame_spacing_ps": dt_ps,
                "bound_frames": int(mask.sum()),
                "occupancy": float(mask.mean()),
                "n_events": len(copy_events),
                "total_bound_ns": float(mask.sum() * dt_ps / 1000.0),
                "mean_event_duration_ns": (
                    float(durations.mean() / 1000.0) if len(durations) else np.nan
                ),
                "maximum_event_duration_ns": (
                    float(durations.max() / 1000.0) if len(durations) else np.nan
                ),
            }
        )

    copies = pd.DataFrame(copy_rows)
    event_columns = [
        "system", "replica", "copy_id", "event_id", "start_frame",
        "end_frame_inclusive", "duration_ns", "left_censored", "right_censored",
    ]
    event_frame = pd.DataFrame(event_rows, columns=event_columns)
    copies.to_csv(output / "copy_metrics.csv", index=False)
    event_frame.to_csv(output / "contact_events.csv", index=False)
    if args.write_timeseries:
        timeseries = pd.DataFrame(
            {
                "frame": np.repeat(np.arange(n_frames), len(ligands)),
                "time_ns": np.repeat(np.arange(n_frames) * dt_ps / 1000.0, len(ligands)),
                "copy_id": np.tile(np.arange(1, len(ligands) + 1), n_frames),
                "bound": bound.reshape(-1),
            }
        )
        timeseries.to_csv(output / "bound_timeseries.csv.gz", index=False)

    manifest = {
        "provenance": execution_provenance(
            "analyse_occupancy", system, replica, [topology, *trajectories], output,
        ),
        "system": system,
        "replica": replica,
        "topology": str(topology),
        "topology_sha256": sha256(topology),
        "trajectories": [str(path) for path in trajectories],
        "trajectory_sha256": {str(path): sha256(path) for path in trajectories},
        "pet_resnames": args.pet_resnames,
        "pet_residue_count": len(pet_residues),
        "ligand_resname": args.ligand_resname,
        "ligand_copy_count": len(ligands),
        "atom_selection": "PET and ligand atoms with mass >1.5 Da (hydrogens excluded)",
        "contact_definition": f"minimum-image heavy-atom distance < {args.cutoff_A} A",
        "cutoff_A": args.cutoff_A,
        "frames": n_frames,
        "frame_spacing_ps": dt_ps,
        "residence_definition": "maximal uninterrupted runs of bound frames; duration = frames x frame spacing",
        "boundary_policy": "boundary events retained for descriptive mean duration and labelled as censored",
    }
    (output / "analysis_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        f"[complete] {system}/{replica}: {len(ligands)} copies, {n_frames} frames, "
        f"{len(event_frame)} events"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
