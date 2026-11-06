#!/usr/bin/env python3
"""Prepare GAFF2/AM1-BCC ligand parameters from SMILES.

The paper workflow selected the dominant Epik microstate at pH 7.4 before
parameterization.  For an entirely open-source rerun, a preselected
microstate may instead be supplied directly as SMILES in the input CSV.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

if __package__:
    from .provenance import execution_provenance
else:
    from provenance import execution_provenance

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem


PENALTY_KEYS = (
    "r_epik_State_Penalty",
    "r_epik_Ionization_Penalty",
    "r_epik_Charging_Adjusted_Penalty",
)


def run(command: list[str], cwd: Path) -> None:
    print("[run]", " ".join(command), flush=True)
    subprocess.run(command, cwd=cwd, check=True)


def executable(name: str) -> str:
    path = shutil.which(name)
    if path is None:
        raise FileNotFoundError(f"Required executable is not on PATH: {name}")
    return path


def safe_name(value: str) -> str:
    name = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip()).strip("_.")
    if not name:
        raise ValueError(f"Invalid empty compound name derived from {value!r}")
    return name


def formal_charge(molecule: Chem.Mol) -> int:
    return int(sum(atom.GetFormalCharge() for atom in molecule.GetAtoms()))


def epik_penalty(molecule: Chem.Mol) -> float:
    properties = molecule.GetPropsAsDict()
    for key in PENALTY_KEYS:
        if key in properties:
            return float(properties[key])
    for key, value in properties.items():
        if "epik" in key.lower() and "penalt" in key.lower():
            return float(value)
    return float("inf")


def select_epik_microstate(smiles: str, name: str, ph: float, pht: float) -> Chem.Mol:
    structconvert = executable("structconvert")
    epikx = executable("epikx")
    with tempfile.TemporaryDirectory(prefix=f"pnsp-epik-{name}-") as temporary:
        work = Path(temporary)
        (work / "input.smi").write_text(f"{smiles} {name}\n", encoding="utf-8")
        run([structconvert, "input.smi", "input.mae"], work)
        run(
            [
                epikx,
                "-WAIT",
                "-ph",
                str(ph),
                "-pht",
                str(pht),
                "-ms",
                "16",
                "-q_lo",
                "-2",
                "-q_hi",
                "2",
                "input.mae",
                "output.mae",
            ],
            work,
        )
        candidates = sorted(work.rglob("*.mae*"), key=lambda path: path.stat().st_mtime)
        candidates = [path for path in candidates if path.name != "input.mae"]
        if not candidates:
            raise RuntimeError("Epik produced no Maestro output")
        run([structconvert, str(candidates[-1]), "states.sdf"], work)
        states = [mol for mol in Chem.SDMolSupplier(str(work / "states.sdf"), removeHs=False) if mol]
        if not states:
            raise RuntimeError("No valid Epik states could be read")
        states.sort(key=lambda mol: (epik_penalty(mol), abs(formal_charge(mol))))
        return Chem.Mol(states[0])


def embed_selected_state(molecule: Chem.Mol, seed_text: str) -> Chem.Mol:
    smiles = Chem.MolToSmiles(Chem.RemoveHs(molecule), isomericSmiles=True)
    embedded = Chem.AddHs(Chem.MolFromSmiles(smiles))
    seed = int(hashlib.sha256(seed_text.encode()).hexdigest()[:8], 16) % 2_000_000_000
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = seed
    if AllChem.EmbedMolecule(embedded, parameters) != 0:
        raise RuntimeError(f"RDKit could not embed selected state: {smiles}")
    if AllChem.MMFFHasAllMoleculeParams(embedded):
        AllChem.MMFFOptimizeMolecule(embedded)
    else:
        AllChem.UFFOptimizeMolecule(embedded)
    return embedded


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def prepare_one(
    row: pd.Series,
    output_root: Path,
    microstate_source: str,
    ph: float,
    pht: float,
    source_csv: Path | None = None,
) -> None:
    name = safe_name(str(row["name"]))
    input_smiles = str(row["smiles"]).strip()
    molecule = Chem.MolFromSmiles(input_smiles)
    if molecule is None:
        raise ValueError(f"RDKit could not parse {name}: {input_smiles}")

    if microstate_source == "epik":
        selected = select_epik_microstate(input_smiles, name, ph, pht)
    else:
        selected = molecule
    selected_smiles = Chem.MolToSmiles(Chem.RemoveHs(selected), isomericSmiles=True)
    selected_charge = formal_charge(selected)
    if "formal_charge" in row and pd.notna(row["formal_charge"]):
        declared_charge = int(row["formal_charge"])
        if declared_charge != selected_charge:
            raise ValueError(
                f"{name}: declared formal_charge={declared_charge}, but selected "
                f"SMILES has formal charge {selected_charge}"
            )

    output = output_root / name
    output.mkdir(parents=True, exist_ok=True)
    embedded = embed_selected_state(selected, selected_smiles)
    embedded.SetProp("_Name", name)
    sdf = output / "selected_microstate.sdf"
    writer = Chem.SDWriter(str(sdf))
    writer.write(embedded)
    writer.close()

    obabel = executable("obabel")
    antechamber = executable("antechamber")
    parmchk2 = executable("parmchk2")
    tleap = executable("tleap")
    run([obabel, sdf.name, "-O", "selected_microstate.mol2"], output)
    run(
        [
            antechamber,
            "-i",
            "selected_microstate.mol2",
            "-fi",
            "mol2",
            "-o",
            "LIG.mol2",
            "-fo",
            "mol2",
            "-c",
            "bcc",
            "-at",
            "gaff2",
            "-rn",
            "LIG",
            "-nc",
            str(selected_charge),
            "-s",
            "2",
            "-pf",
            "y",
        ],
        output,
    )
    run([parmchk2, "-i", "LIG.mol2", "-f", "mol2", "-o", "LIG.frcmod"], output)
    leap_input = output / "tleap_ligand.in"
    leap_input.write_text(
        "source leaprc.gaff2\n"
        "loadamberparams LIG.frcmod\n"
        "LIG = loadmol2 LIG.mol2\n"
        "saveoff LIG LIG.lib\n"
        "savepdb LIG LIG.pdb\n"
        "quit\n",
        encoding="utf-8",
    )
    run([tleap, "-f", leap_input.name], output)

    products = ["LIG.mol2", "LIG.frcmod", "LIG.lib", "LIG.pdb"]
    missing = [filename for filename in products if not (output / filename).is_file()]
    if missing:
        raise RuntimeError(f"{name}: missing Amber products: {missing}")
    manifest = {
        "provenance": execution_provenance(
            "build_contaminant", name, None,
            [source_csv] if source_csv is not None else [], output,
        ),
        "name": name,
        "input_smiles": input_smiles,
        "selected_smiles": selected_smiles,
        "formal_charge": selected_charge,
        "microstate_source": microstate_source,
        "epik_ph": ph if microstate_source == "epik" else None,
        "epik_pht": pht if microstate_source == "epik" else None,
        "force_field": "GAFF2",
        "charge_method": "AM1-BCC",
        "files": {filename: sha256(output / filename) for filename in products},
        "executables": {
            key: executable(key)
            for key in ("obabel", "antechamber", "parmchk2", "tleap")
        },
    }
    (output / "ligand_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"[complete] {name}: charge {selected_charge:+d}, {selected_smiles}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path, help="CSV with name, smiles, and optional formal_charge")
    parser.add_argument("--output", type=Path, default=Path("build/ligands"))
    parser.add_argument(
        "--microstate-source",
        choices=("input", "epik"),
        default="input",
        help="Use the supplied SMILES as the selected state, or run Epik as in the paper",
    )
    parser.add_argument("--ph", type=float, default=7.4)
    parser.add_argument("--pht", type=float, default=2.0)
    args = parser.parse_args()

    table = pd.read_csv(args.csv)
    required = {"name", "smiles"}
    missing = required - set(table.columns)
    if missing:
        raise ValueError(f"Input CSV is missing columns: {sorted(missing)}")
    names = table["name"].map(lambda value: safe_name(str(value)))
    if names.duplicated().any():
        raise ValueError("Compound names are not unique after filename sanitization")
    args.output.mkdir(parents=True, exist_ok=True)
    for _, row in table.iterrows():
        prepare_one(row, args.output, args.microstate_source, args.ph, args.pht, args.csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
