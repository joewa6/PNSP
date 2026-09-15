#!/usr/bin/env python3
"""Assemble, solvate, and neutralize one PET/contaminant replica."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import shutil
import subprocess
from pathlib import Path

import numpy as np


PET_FILES = ("NP_correct.lib", "NP_charged.frcmod", "NP_neut.frcmod", "NP_dry.pdb")
LIGAND_FILES = ("LIG.lib", "LIG.frcmod", "LIG.pdb", "ligand_manifest.json")


def executable(name: str) -> str:
    path = shutil.which(name)
    if path is None:
        raise FileNotFoundError(f"Required executable is not on PATH: {name}")
    return path


def run(command: list[str], cwd: Path, stdin: Path | None = None) -> None:
    print("[run]", " ".join(command), flush=True)
    if stdin is None:
        subprocess.run(command, cwd=cwd, check=True)
    else:
        with stdin.open() as handle:
            subprocess.run(command, cwd=cwd, stdin=handle, check=True)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def pdb_atoms(path: Path) -> list[tuple[str, tuple[str, str, str], np.ndarray]]:
    atoms = []
    with path.open() as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            residue = line[17:20].strip()
            key = (line[21:22].strip() or "-", line[22:26].strip(), line[26:27].strip())
            xyz = np.array(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                dtype=float,
            )
            atoms.append((residue, key, xyz))
    return atoms


def topology_charge(path: Path) -> float:
    values: list[float] = []
    reading = False
    with path.open() as handle:
        for line in handle:
            if line.startswith("%FLAG CHARGE"):
                reading = True
                continue
            if reading and line.startswith("%FORMAT"):
                continue
            if reading and line.startswith("%FLAG"):
                break
            if reading:
                values.extend(float(value) for value in line.split())
    if not values:
        raise RuntimeError(f"Could not read CHARGE from {path}")
    return float(sum(values) / 18.2223)


def coordinate_box(path: Path) -> list[float]:
    lines = path.read_text().splitlines()
    atom_count = int(lines[1].split()[0])
    values = [float(value) for line in lines[2:] for value in line.split()]
    box = values[3 * atom_count : 3 * atom_count + 6]
    if len(box) != 6:
        raise RuntimeError(f"Could not read periodic box from {path}")
    return box


def minimum_ligand_pet_distance(path: Path, expected_copies: int) -> tuple[float, int]:
    atoms = pdb_atoms(path)
    pet = np.vstack([xyz for residue, _, xyz in atoms if residue != "LIG"])
    ligand_groups: dict[tuple[str, str, str], list[np.ndarray]] = {}
    for residue, key, xyz in atoms:
        if residue == "LIG":
            ligand_groups.setdefault(key, []).append(xyz)
    if len(ligand_groups) != expected_copies:
        raise RuntimeError(
            f"Packmol produced {len(ligand_groups)} LIG residues; expected {expected_copies}"
        )
    minimum = math.inf
    for coordinates in ligand_groups.values():
        ligand = np.vstack(coordinates)
        for start in range(0, len(pet), 1000):
            delta = ligand[:, None, :] - pet[None, start : start + 1000, :]
            minimum = min(minimum, float(np.sqrt(np.sum(delta * delta, axis=2)).min()))
    return minimum, len(ligand_groups)


def require_files(directory: Path, filenames: tuple[str, ...]) -> None:
    missing = [filename for filename in filenames if not (directory / filename).is_file()]
    if missing:
        raise FileNotFoundError(f"Missing files under {directory}: {missing}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pet", required=True, type=Path, help="Directory containing PET assets")
    parser.add_argument("--ligand", required=True, type=Path, help="Prepared ligand directory")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--copies", type=int, default=10)
    parser.add_argument("--seed", type=int, required=True, help="Explicit positive Packmol seed")
    parser.add_argument(
        "--pet-unit", default="NP", help="Unit name stored in NP_correct.lib"
    )
    parser.add_argument("--pet-charge", type=int, default=-50)
    parser.add_argument("--solvent-buffer", type=float, default=11.0, help="Isotropic TIP3P buffer in Angstrom")
    parser.add_argument("--packmol-tolerance", type=float, default=2.0, help="Minimum initial all-atom separation in Angstrom")
    args = parser.parse_args()
    if args.copies < 1 or args.seed < 1:
        raise ValueError("--copies and --seed must be positive")

    pet = args.pet.resolve()
    ligand = args.ligand.resolve()
    output = args.output.resolve()
    require_files(pet, PET_FILES)
    require_files(ligand, LIGAND_FILES)
    ligand_manifest = json.loads((ligand / "ligand_manifest.json").read_text())
    ligand_charge = int(ligand_manifest["formal_charge"])
    output.mkdir(parents=True, exist_ok=True)

    tleap = executable("tleap")
    packmol = executable("packmol")
    clean_input = output / "tleap_pet_clean.in"
    clean_input.write_text(
        f"loadoff {pet / 'NP_correct.lib'}\n"
        f"PET = copy {args.pet_unit}\n"
        "savepdb PET PET_clean.pdb\n"
        "quit\n",
        encoding="utf-8",
    )
    run([tleap, "-f", clean_input.name], output)

    pet_coordinates = np.vstack([xyz for _, _, xyz in pdb_atoms(output / "PET_clean.pdb")])
    lower = pet_coordinates.min(axis=0)
    upper = pet_coordinates.max(axis=0)
    packmol_input = output / "packmol.inp"
    packmol_input.write_text(
        f"tolerance {args.packmol_tolerance}\n"
        f"seed {args.seed}\n"
        "filetype pdb\n"
        "output combined.pdb\n\n"
        "structure PET_clean.pdb\n"
        "  number 1\n"
        "  fixed 0. 0. 0. 0. 0. 0.\n"
        "end structure\n\n"
        f"structure {ligand / 'LIG.pdb'}\n"
        f"  number {args.copies}\n"
        f"  inside box {' '.join(f'{value:.3f}' for value in np.r_[lower, upper])}\n"
        "end structure\n",
        encoding="utf-8",
    )
    run([packmol], output, stdin=packmol_input)
    minimum, observed_copies = minimum_ligand_pet_distance(
        output / "combined.pdb", args.copies
    )
    if minimum < args.packmol_tolerance - 0.002:
        raise RuntimeError(
            f"Initial ligand/PET distance {minimum:.3f} A violates Packmol tolerance"
        )

    with (output / "combined.pdb").open() as source, (output / "ligands_only.pdb").open("w") as target:
        for line in source:
            if line.startswith(("ATOM", "HETATM")) and line[17:20].strip() == "LIG":
                target.write(line)

    estimated_charge = args.pet_charge + args.copies * ligand_charge
    ion_line = ""
    if estimated_charge < 0:
        ion_line = "addionsrand SYS Na+ 0\n"
    elif estimated_charge > 0:
        ion_line = "addionsrand SYS Cl- 0\n"
    final_input = output / "tleap_system.in"
    final_input.write_text(
        "source leaprc.gaff2\n"
        "source leaprc.water.tip3p\n"
        f"loadamberparams {pet / 'NP_charged.frcmod'}\n"
        f"loadamberparams {pet / 'NP_neut.frcmod'}\n"
        f"loadoff {pet / 'NP_correct.lib'}\n"
        f"PET = copy {args.pet_unit}\n"
        f"loadamberparams {ligand / 'LIG.frcmod'}\n"
        f"loadoff {ligand / 'LIG.lib'}\n"
        "LIGS = loadpdb ligands_only.pdb\n"
        "SYS = combine {PET LIGS}\n"
        f"solvatebox SYS TIP3PBOX {args.solvent_buffer:.3f} iso\n"
        f"{ion_line}"
        "saveamberparm SYS system.prmtop system.inpcrd\n"
        "savepdb SYS system_solvated.pdb\n"
        "quit\n",
        encoding="utf-8",
    )
    run([tleap, "-f", final_input.name], output)

    products = ["system.prmtop", "system.inpcrd", "system_solvated.pdb"]
    require_files(output, tuple(products))
    final_charge = topology_charge(output / "system.prmtop")
    # AMBER stores partial charges at finite precision in the topology.  A formally
    # neutral system can consequently sum to a small residual (the canonical PET
    # example gives approximately -0.01 e after counterion addition).
    if abs(final_charge) > 0.05:
        raise RuntimeError(
            f"Neutralized topology has residual charge {final_charge:.6f} e "
            "(allowed rounding tolerance: 0.05 e)"
        )
    box = coordinate_box(output / "system.inpcrd")
    manifest = {
        "pet_charge_e": args.pet_charge,
        "ligand_charge_e": ligand_charge,
        "ligand_copies": observed_copies,
        "packmol_seed": args.seed,
        "packmol_tolerance_A": args.packmol_tolerance,
        "initial_minimum_ligand_pet_distance_A": minimum,
        "solvent_model": "TIP3P",
        "solvation_rule": "tleap solvatebox with isotropic solvent buffer",
        "solvent_buffer_A": args.solvent_buffer,
        "neutralization_rule": "counterions only; addionsrand to net charge zero",
        "charge_validation_tolerance_e": 0.05,
        "final_topology_charge_e": final_charge,
        "initial_box_A": box,
        "files": {filename: sha256(output / filename) for filename in products},
    }
    (output / "system_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        f"[complete] {observed_copies} copies; min initial distance={minimum:.3f} A; "
        f"box={box[0]:.2f} x {box[1]:.2f} x {box[2]:.2f} A"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
