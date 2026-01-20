#!/usr/bin/env python3
"""
Stage A builder entrypoint - CLEAN REFACTORED VERSION

This is the simple scaffold builder that creates placeholder topology files.
For full Epik-based microstate enumeration and complex building, use epik_builder.py

Contract:
  - CLI: python src/build_test.py --smiles "..." --output /path/to/out
  - Produces: out/run_manifest.json, out/topology.prmtop, out/coords.inpcrd
"""

from __future__ import annotations

import argparse
import datetime as _dt
import json
import os
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile
import yaml
from dataclasses import dataclass
from typing import Any, List, Optional, Tuple

PROJECT_NAME = "PET_SCREENING"
BUILDER_VERSION = "0.1.0"

# Epik defaults
DEFAULT_PH = 7.4
DEFAULT_PHT = 2.0
DEFAULT_MAX_STATES = 16
DEFAULT_Q_LO, DEFAULT_Q_HI = -2, 2


def iso_now_utc() -> str:
    return _dt.datetime.now(tz=_dt.timezone.utc).isoformat(timespec="seconds")


def safe_mkdir(p: pathlib.Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def relpath_str(path: pathlib.Path, start: pathlib.Path) -> str:
    return os.path.relpath(path.resolve(), start.resolve())


def get_git_commit(repo_root: pathlib.Path) -> str:
    try:
        out = subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=str(repo_root),
            stderr=subprocess.DEVNULL,
            text=True,
        ).strip()
        if re.fullmatch(r"[0-9a-f]{7,40}", out):
            return out
    except Exception:
        pass
    return "unknown"


def find_repo_root(start: pathlib.Path) -> pathlib.Path:
    cur = start.resolve()
    for _ in range(25):
        if (cur / ".git").exists():
            return cur
        if cur.parent == cur:
            break
        cur = cur.parent
    return start.resolve()


def load_np_profile(yaml_path: pathlib.Path) -> dict[str, Any]:
    """Load and validate NP profile YAML configuration."""
    if not yaml_path.exists():
        raise FileNotFoundError(f"NP profile not found: {yaml_path}")
    
    with open(yaml_path, 'r') as f:
        profile = yaml.safe_load(f)
    
    # Validate required fields
    required_fields = ['profile_id', 'charge', 'files']
    for field in required_fields:
        if field not in profile:
            raise ValueError(f"NP profile missing required field: {field}")
    
    required_files = ['lib', 'charged_frcmod', 'neut_frcmod', 'dry_pdb']
    for file_key in required_files:
        if file_key not in profile['files']:
            raise ValueError(f"NP profile missing required file: {file_key}")
    
    return profile


def resolve_np_files(np_profile_dir: pathlib.Path, profile: dict[str, Any]) -> dict[str, pathlib.Path]:
    """Resolve NP file paths from profile configuration to absolute paths."""
    np_files = {}
    for key, rel_path in profile['files'].items():
        abs_path = (np_profile_dir / rel_path).resolve()
        if not abs_path.exists():
            raise FileNotFoundError(f"NP file not found: {abs_path}")
        np_files[key] = abs_path
    return np_files


def resolve_default_np_files(stageA_root: pathlib.Path) -> dict[str, pathlib.Path] | None:
    """Return built-in NP files if present; otherwise None."""
    np_dir = stageA_root / "NP"
    req = {
        "lib": np_dir / "NP_correct.lib",
        "charged_frcmod": np_dir / "NP_charged.frcmod",
        "neut_frcmod": np_dir / "NP_neut.frcmod",
        "dry_pdb": np_dir / "NP_dry.pdb",
    }
    if all(p.exists() for p in req.values()):
        return {k: v.resolve() for k, v in req.items()}
    return None


def generate_run_id(
    project: str,
    np_profile_id: str,
    ligand_set: str,
    git_commit: str
) -> str:
    """
    Generate unique run identifier.
    
    Format: <timestamp>__<charge_state>__<ligand_set>__git-<commit>
    Example: 2026-01-14T16-05-12Z__c50__example_ligands__git-a1b2c3d
    """
    timestamp = _dt.datetime.now(tz=_dt.timezone.utc).strftime("%Y-%m-%dT%H-%M-%SZ")
    git_short = git_commit[:7] if len(git_commit) >= 7 else git_commit
    
    # Sanitize ligand_set name (remove special chars, spaces)
    ligand_set_clean = re.sub(r'[^a-zA-Z0-9_-]', '_', ligand_set)
    
    return f"{timestamp}__{np_profile_id}__{ligand_set_clean}__git-{git_short}"


def normalize_smiles(smiles: str) -> str:
    s = smiles.strip()
    if not s:
        raise ValueError("Empty SMILES after stripping whitespace.")
    try:
        from rdkit import Chem  # type: ignore
        mol = Chem.MolFromSmiles(s)
        if mol is None:
            raise ValueError("RDKit could not parse SMILES.")
        s_can = Chem.MolToSmiles(mol, isomericSmiles=True)
        return s_can
    except ImportError:
        return s
    except Exception as e:
        raise ValueError(f"Invalid SMILES (RDKit parse failed): {e}") from e


def ligand_id_from_smiles(smiles_norm: str) -> str:
    import hashlib
    h = hashlib.sha256(smiles_norm.encode("utf-8")).hexdigest()
    return f"lig_{h[:12]}"


# ============================================================================
# Epik Pipeline Functions
# ============================================================================

def runp(args: List[str], cwd: Optional[pathlib.Path] = None, timeout: int = 3600, stdin=None) -> subprocess.CompletedProcess:
    """Run subprocess with nice output formatting."""
    here = cwd or pathlib.Path.cwd()
    cmd_str = " ".join(args)
    print(f"\n--- RUN ({here}): {cmd_str}\n")
    try:
        res = subprocess.run(
            args,
            cwd=str(here),
            text=True,
            capture_output=True,
            timeout=timeout,
            check=True,
            stdin=stdin,
        )
        if res.stdout:
            print("---- stdout ----")
            print(res.stdout.strip(), "\n")
        if res.stderr:
            print("---- stderr ----")
            print(res.stderr.strip(), "\n")
        return res
    except subprocess.CalledProcessError as e:
        out = (e.stdout or "").strip()
        err = (e.stderr or "").strip()
        print("❌ Shell command failed.")
        if out:
            print("---- stdout ----")
            print(out)
        if err:
            print("\n---- stderr ----")
            print(err)
        raise


@dataclass
class EpikState:
    """Represents one Epik microstate."""
    mol: Any  # Chem.Mol
    formal_charge: Optional[int]
    penalty: Optional[float]


_PEN_KEYS = [
    "r_epik_State_Penalty",
    "r_epik_Ionization_Penalty",
    "r_epik_Charging_Adjusted_Penalty",
]


def _extract_penalty(props: dict) -> Optional[float]:
    """Extract Epik penalty from mol properties."""
    for k in _PEN_KEYS:
        if k in props:
            try:
                return float(props[k])
            except Exception:
                pass
    for k, v in props.items():
        lk = str(k).lower()
        if "epik" in lk and "penalt" in lk:
            try:
                return float(v)
            except Exception:
                pass
    return None


def _parse_epik_sdf(sdf_path: pathlib.Path) -> List[EpikState]:
    """Parse Epik SDF output to extract states with penalties."""
    from rdkit import Chem
    states: List[EpikState] = []
    try:
        sup = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
    except Exception:
        return states
    for m in sup:
        if m is None:
            continue
        try:
            chg = int(sum(a.GetFormalCharge() for a in m.GetAtoms()))
        except Exception:
            chg = None
        props = {k: str(v) for k, v in m.GetPropsAsDict().items()}
        pen = _extract_penalty(props)
        states.append(EpikState(mol=Chem.Mol(m), formal_charge=chg, penalty=pen))
    return states


def _find_epik_out(workdir: pathlib.Path) -> Optional[pathlib.Path]:
    """Locate Epik's Maestro output after a run."""
    for nm in ("out.mae", "out.maegz"):
        p = workdir / nm
        if p.exists() and p.stat().st_size > 0:
            return p
    jobdirs = [d for d in workdir.iterdir() if d.is_dir() and d.name.startswith(".")]
    for jd in jobdirs:
        cands = sorted(jd.glob("*.mae*"), key=lambda p: p.stat().st_mtime, reverse=True)
        if cands:
            return cands[0]
    cands = sorted((p for p in workdir.glob("*.mae*") if p.name != "in.mae"),
                   key=lambda p: p.stat().st_mtime, reverse=True)
    return cands[0] if cands else None


def epikx_best_state_as_sdf(
    smiles: str,
    name: str,
    ph: float = DEFAULT_PH,
    pht: float = DEFAULT_PHT,
    max_states: int = DEFAULT_MAX_STATES,
    qlo: int = DEFAULT_Q_LO,
    qhi: int = DEFAULT_Q_HI,
) -> Tuple[Any, int]:  # Returns (Chem.Mol, formal_charge)
    """Run Epik on SMILES, return (RDKit Mol with 3D coords, formal_charge) for lowest-penalty state."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    
    with tempfile.TemporaryDirectory(prefix="epikx_") as td:
        tdir = pathlib.Path(td)
        (tdir / "in.smi").write_text(f"{smiles} {name}\n")

        # SMILES -> MAE
        runp(["structconvert", "in.smi", "in.mae"], cwd=tdir)

        # Epik enumeration
        runp([
            "epikx",
            "-WAIT",
            "-ph", str(ph), "-pht", str(pht),
            "-ms", str(max_states), "-q_lo", str(qlo), "-q_hi", str(qhi),
            "in.mae", "out.mae"
        ], cwd=tdir)

        epik_out = _find_epik_out(tdir)
        if not epik_out:
            print("!! Epik temp dir contents:")
            for p in sorted(tdir.iterdir()):
                print("  -", p.name)
            raise RuntimeError("EpikX produced no usable Maestro output")

        # Flatten if compressed
        if epik_out.suffix == ".maegz":
            runp(["structconvert", epik_out.name, "out_flat.mae"], cwd=tdir)
            epik_out = tdir / "out_flat.mae"

        # Convert to SDF and parse states
        runp(["structconvert", epik_out.name, "out.sdf"], cwd=tdir)
        states = _parse_epik_sdf(tdir / "out.sdf")
        if not states:
            raise RuntimeError("No states parsed from Epik SDF")

        # Pick lowest penalty, then smaller |charge|
        states.sort(key=lambda s: (float("inf") if s.penalty is None else s.penalty,
                                   abs(s.formal_charge or 0)))
        best = states[0]
        formal_charge = int(best.formal_charge or 0)

        # Regenerate clean 3D geometry from best state's SMILES
        best_smi = Chem.MolToSmiles(Chem.AddHs(best.mol), isomericSmiles=True)
        m3d = Chem.MolFromSmiles(best_smi)
        m3d = Chem.AddHs(m3d)
        AllChem.EmbedMolecule(m3d, AllChem.ETKDG())
        try:
            AllChem.MMFFOptimizeMolecule(m3d)
        except Exception:
            pass
        m3d.SetProp("_Name", name)
        return m3d, formal_charge


def write_mol_files(mol: Any, out_base: pathlib.Path) -> Tuple[pathlib.Path, pathlib.Path]:
    """Write SDF and MOL2 (via obabel). Returns (sdf_path, mol2_path)."""
    from rdkit import Chem
    
    sdf_path = out_base.with_suffix(".sdf")
    mol2_path = out_base.with_suffix(".mol2")

    w = Chem.SDWriter(str(sdf_path))
    w.write(mol)
    w.close()

    runp(["obabel", str(sdf_path), "-O", str(mol2_path)], cwd=out_base.parent)
    return sdf_path, mol2_path


def prepare_ligand_assets(folder: pathlib.Path, mol_id: str, net_charge: int) -> None:
    """Run antechamber/parmchk2/tleap in folder to create LIG assets."""
    print(f"\n🧪 STEP 1: Preparing Ligand {mol_id}...")
    runp(
        [
            "antechamber",
            "-i", f"{mol_id}.mol2", "-fi", "mol2",
            "-o", "LIG.mol2", "-fo", "mol2",
            "-c", "bcc",
            "-s", "2",
            "-rn", "LIG",
            "-nc", str(net_charge),
        ],
        cwd=folder,
    )
    runp(["parmchk2", "-i", "LIG.mol2", "-f", "mol2", "-o", "LIG.frcmod"], cwd=folder)

    tleap_in = folder / "tleap_ligand.in"
    tleap_in.write_text(
        "source leaprc.gaff2\n"
        "loadamberparams LIG.frcmod\n"
        "LIG = loadmol2 LIG.mol2\n"
        "saveoff LIG LIG.lib\n"
        "savepdb LIG LIG.pdb\n"
        "quit\n"
    )
    runp(["tleap", "-f", "tleap_ligand.in"], cwd=folder)
    print(f"✅ Created {mol_id}.mol2 and prepared ligand assets in {folder}")


# ============================================================================
# Main Builder Function
# ============================================================================

@dataclass(frozen=True)
class BuildOutputs:
    topology: pathlib.Path
    coords: pathlib.Path


def run_builder(
    smiles_norm: str,
    out_dir: pathlib.Path,
    np_files: dict[str, pathlib.Path] | None = None,
    templates_dir: pathlib.Path | None = None,
    dry_run: bool = False,
) -> BuildOutputs:
    """
    Full Epik-based builder with NP complex assembly.
    
    Args:
        smiles_norm: Normalized SMILES string
        out_dir: Output directory for build artifacts
        np_files: Dictionary of NP force field files (lib, frcmod, pdb). If None, builds ligand-only.
        templates_dir: Optional templates directory (currently unused)
        dry_run: If True, skip actual file generation
    
    Workflow:
      1. Epik microstate enumeration
      2. Ligand preparation (antechamber/parmchk2/tleap)
      3. NP complex assembly with packmol
      4. Solvation and final AMBER topology generation
    """
    top = out_dir / "system.prmtop"
    crd = out_dir / "system.inpcrd"

    if dry_run:
        return BuildOutputs(topology=top, coords=crd)

    # Try to use the full Epik + NP pipeline
    try:
        print(f"🧬 Building NP complex with Epik pipeline for SMILES: {smiles_norm}")
        
        # Generate a simple name for this ligand
        ligand_name = ligand_id_from_smiles(smiles_norm)
        
        # STEP 1: Run Epik to get best microstate
        print(f"\n🔬 STEP 1: Epik microstate enumeration...")
        best_mol, formal_charge = epikx_best_state_as_sdf(
            smiles=smiles_norm,
            name=ligand_name,
        )
        
        # STEP 2: Write molecule files and prepare ligand
        print(f"\n📝 STEP 2: Preparing ligand files...")
        out_base = out_dir / ligand_name
        write_mol_files(best_mol, out_base)
        prepare_ligand_assets(out_dir, ligand_name, net_charge=formal_charge)
        
        # STEP 3: Build NP complex (if NP files provided)
        print(f"\n🧩 STEP 3: Building NP complex...")
        if np_files:
            _build_np_complex(out_dir, np_files, ligand_copies=10)
        else:
            print("   ⚠️  No NP files provided, building ligand-only system")
            _build_ligand_only_system(out_dir)
        
        print(f"✅ Full Epik + NP pipeline complete: {top.name}, {crd.name}")
        
    except ImportError as e:
        raise RuntimeError(
            f"Current profile does not contain rdkit,epik. "
            f"Install required dependencies to use Epik pipeline. Error: {e}"
        ) from e
    except Exception as e:
        print(f"⚠️  Epik pipeline failed: {e}")
        import traceback
        traceback.print_exc()
        raise

    return BuildOutputs(topology=top, coords=crd)


def _build_np_complex(
    out_dir: pathlib.Path, 
    np_files: dict[str, pathlib.Path],
    ligand_copies: int = 10
) -> None:
    """
    Build NP complex with ligands using packmol and tleap.
    
    Args:
        out_dir: Output directory for build artifacts
        np_files: Dictionary with keys 'lib', 'charged_frcmod', 'neut_frcmod', 'dry_pdb'
                  Values are absolute paths to NP force field files
        ligand_copies: Number of ligand molecules to pack around NP
    
    Steps:
      1. Load NP structure
      2. Use packmol to add ligands around NP
      3. Solvate and add ions
      4. Generate final AMBER files
    """
    print(f"   Using NP files from: {np_files['lib'].parent}")
    
    # Check if NP_correct.lib has a top-level NP unit
    try:
        lib_txt = np_files['lib'].read_text(errors="ignore")
        np_has_unit = ("!entry.NP.unit" in lib_txt) or ("!entry.NP.unit.atoms" in lib_txt)
    except Exception:
        np_has_unit = False
    
    # Step 3a: Create NP_clean.pdb
    print("   - Creating NP_clean.pdb...")
    tleap_np_clean = out_dir / "tleap_np_clean.in"
    if np_has_unit:
        tleap_np_clean.write_text(
            f"loadoff {np_files['lib']}\n"
            f"NPu = copy NP\n"
            f"savepdb NPu NP_clean.pdb\n"
            f"quit\n"
        )
    else:
        tleap_np_clean.write_text(
            f"loadoff {np_files['lib']}\n"
            f"NPu = loadpdb {np_files['dry_pdb']}\n"
            f"savepdb NPu NP_clean.pdb\n"
            f"quit\n"
        )
    runp(["tleap", "-f", "tleap_np_clean.in"], cwd=out_dir)
    
    # Step 3b: Extract NP bounding box
    print("   - Calculating NP bounding box...")
    runp(["bash", "-c", "grep 'ATOM' NP_clean.pdb | awk '{print $6,$7,$8}' > bbox.txt"], cwd=out_dir)
    
    coords = []
    with open(out_dir / "bbox.txt") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 3:
                x, y, z = map(float, parts)
                coords.append((x, y, z))
    
    if not coords:
        raise RuntimeError("Failed to extract NP coordinates from NP_clean.pdb")
    
    xs, ys, zs = zip(*coords)
    buffer = 20.0
    box = (
        min(xs) - buffer, min(ys) - buffer, min(zs) - buffer,
        max(xs) + buffer, max(ys) + buffer, max(zs) + buffer
    )
    cx = sum(xs) / len(xs)
    cy = sum(ys) / len(ys)
    cz = sum(zs) / len(zs)
    
    # Step 3c: Packmol - place ligands around NP
    print(f"   - Running packmol to add {ligand_copies} ligand copies...")
    packmol_input = out_dir / "packmol.inp"
    with open(packmol_input, "w") as f:
        f.write(
            f"tolerance 2.0\n"
            f"filetype pdb\n"
            f"output combined.pdb\n"
            f"\n"
            f"structure NP_clean.pdb\n"
            f"  number 1\n"
            f"  fixed {cx:.3f} {cy:.3f} {cz:.3f} 0. 0. 0.\n"
            f"end structure\n"
            f"\n"
            f"structure LIG.pdb\n"
            f"  number {ligand_copies}\n"
            f"  inside box {' '.join(map(str, box))}\n"
            f"end structure\n"
        )
    
    with open(packmol_input, "r") as fin:
        runp(["packmol"], cwd=out_dir, timeout=600, stdin=fin)
    
    if not (out_dir / "combined.pdb").exists():
        raise RuntimeError("Packmol failed to produce combined.pdb")
    
    # Step 3d: Extract ligand atoms only
    print("   - Extracting ligand atoms...")
    lig_only = out_dir / "ligands_only.pdb"
    with open(out_dir / "combined.pdb") as fin, open(lig_only, "w") as fout:
        for line in fin:
            if (line.startswith("ATOM") or line.startswith("HETATM")) and line[17:20].strip() == "LIG":
                fout.write(line)
    
    # Step 3e: Final tleap to build solvated system
    print("   - Building final solvated system with tleap...")
    tleap_final = out_dir / "tleap_final.in"
    if np_has_unit:
        np_block = f"loadoff {np_files['lib']}\nNPu = copy NP\n"
    else:
        np_block = f"loadoff {np_files['lib']}\nNPu = loadpdb {np_files['dry_pdb']}\n"
    
    tleap_final.write_text(
        f"source leaprc.gaff2\n"
        f"source leaprc.water.tip3p\n"
        f"loadamberparams {np_files['charged_frcmod']}\n"
        f"loadamberparams {np_files['neut_frcmod']}\n"
        f"{np_block}"
        f"loadamberparams LIG.frcmod\n"
        f"loadoff LIG.lib\n"
        f"LIGS = loadpdb ligands_only.pdb\n"
        f"SYS = combine {{NPu LIGS}}\n"
        f"solvatebox SYS TIP3PBOX 12.0\n"
        f"addionsrand SYS Na+ 0\n"
        f"addionsrand SYS Cl- 0\n"
        f"saveamberparm SYS system.prmtop system.inpcrd\n"
        f"savepdb SYS system_solvated.pdb\n"
        f"quit\n"
    )
    runp(["tleap", "-f", "tleap_final.in"], cwd=out_dir)
    
    print(f"   ✅ NP complex built successfully")


def _build_ligand_only_system(out_dir: pathlib.Path) -> None:
    """Build simple ligand-only system (no NP)."""
    import subprocess
    
    print("   Building ligand-only system...")
    tleap_script = out_dir / "tleap_ligand_only.in"
    tleap_script.write_text(
        "source leaprc.gaff2\n"
        "source leaprc.water.tip3p\n"
        "loadamberparams LIG.frcmod\n"
        "loadoff LIG.lib\n"
        "LIG_mol = loadpdb LIG.pdb\n"
        "solvatebox LIG_mol TIP3PBOX 12.0\n"
        "addionsrand LIG_mol Na+ 0\n"
        "addionsrand LIG_mol Cl- 0\n"
        "saveamberparm LIG_mol system.prmtop system.inpcrd\n"
        "savepdb LIG_mol system_solvated.pdb\n"
        "quit\n"
    )
    subprocess.run(
        ["tleap", "-f", "tleap_ligand_only.in"],
        cwd=str(out_dir),
        capture_output=True,
        check=True,
    )


def _build_simple_ligand_system(smiles_norm: str, out_dir: pathlib.Path) -> None:
    """Fallback: build simple ligand system if Epik fails."""
    # Check if LIG files already exist from partial run
    if (out_dir / "LIG.lib").exists() and (out_dir / "LIG.frcmod").exists():
        _build_ligand_only_system(out_dir)
    else:
        # Complete fallback to placeholder
        top = out_dir / "system.prmtop"
        crd = out_dir / "system.inpcrd"
        _write_placeholder_files(smiles_norm, top, crd)


def _write_placeholder_files(smiles_norm: str, top: pathlib.Path, crd: pathlib.Path) -> None:
    """Write simple placeholder topology and coords files."""
    top.write_text(
        "\n".join([
            "%VERSION  VERSION_STUB",
            f"%COMMENT  StageA placeholder prmtop for {smiles_norm}",
            "%FLAG TITLE",
            "%FORMAT(20a4)",
            "STAGEA_STUB",
            "",
        ]),
        encoding="utf-8",
    )
    crd.write_text(
        "\n".join([
            "StageA placeholder inpcrd",
            f"SMILES {smiles_norm}",
            "   0.0000000   0.0000000   0.0000000",
            "",
        ]),
        encoding="utf-8",
    )


def write_manifest(
    out_dir: pathlib.Path,
    smiles_in: str,
    smiles_norm: str,
    ligand_norm: str,
    outputs: BuildOutputs,
    git_commit: str,
    status: str,
) -> pathlib.Path:
    manifest_path = out_dir / "run_manifest.json"

    data: dict[str, Any] = {
        "project": PROJECT_NAME,
        "ligand_norm": ligand_norm,
        "smiles": smiles_in,
        "stageA": {
            "status": status,
            "builder_version": BUILDER_VERSION,
            "git_commit": git_commit,
            "timestamp": iso_now_utc(),
        },
        "paths": {
            "topology": relpath_str(outputs.topology, out_dir),
            "coords": relpath_str(outputs.coords, out_dir),
        },
    }

    manifest_path.write_text(json.dumps(data, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    return manifest_path


def parse_args(argv: list[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="stageA_builder",
        description="Stage A builder: construct ligand system inputs and emit run manifest.",
    )
    
    # Mutually exclusive group for input mode
    input_group = p.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--smiles", help="Ligand SMILES string (single molecule mode).")
    input_group.add_argument("--csv", help="CSV file with Name,SMILES columns (batch mode).")
    
    p.add_argument("--output", help="Output directory for single molecule build.")
    p.add_argument("--output-dir", help="Output root directory for batch builds (CSV mode).")
    p.add_argument("--templates", default=None, help="Optional templates dir.")
    p.add_argument("--dry-run", action="store_true", help="Skip actual file generation.")
    return p.parse_args(argv)


def process_csv_batch(
    csv_path: pathlib.Path,
    output_root: pathlib.Path,
    templates_dir: pathlib.Path,
    git_commit: str,
    np_files: dict[str, pathlib.Path] | None,
    dry_run: bool = False,
) -> dict[str, Any]:
    """
    Process multiple molecules from CSV file.
    
    Returns batch summary with success/failure counts.
    """
    import csv
    
    # Validate CSV
    if not csv_path.exists():
        raise FileNotFoundError(f"CSV file not found: {csv_path}")
    
    # Read CSV
    with open(csv_path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        if 'Name' not in reader.fieldnames or 'SMILES' not in reader.fieldnames:
            raise ValueError("CSV must contain 'Name' and 'SMILES' columns")
        rows = list(reader)
    
    if not rows:
        raise ValueError("CSV file is empty")
    
    print(f"📊 Processing {len(rows)} molecules from {csv_path.name}\n")
    
    # Track results
    results = {
        'total': len(rows),
        'successful': [],
        'failed': [],
        'skipped': [],
    }
    
    # Process each molecule
    for idx, row in enumerate(rows, 1):
        name = row.get('Name', '').strip()
        smiles = row.get('SMILES', '').strip()
        
        if not name or not smiles:
            print(f"⏭️  [{idx}/{len(rows)}] Skipping row with empty Name or SMILES")
            results['skipped'].append({
                'row': idx,
                'name': name or 'EMPTY',
                'smiles': smiles or 'EMPTY',
                'reason': 'Empty name or SMILES'
            })
            continue
        
        # Create safe directory name from molecule name
        safe_name = name.replace(' ', '_').replace(',', '').replace('/', '_')
        mol_output_dir = output_root / safe_name
        
        print(f"🔬 [{idx}/{len(rows)}] Processing: {name} ({smiles})")
        
        try:
            # Normalize SMILES
            smiles_norm = normalize_smiles(smiles)
            ligand_norm = ligand_id_from_smiles(smiles_norm)
            
            # Create output directory
            safe_mkdir(mol_output_dir)
            
            # Run builder
            outputs = run_builder(
                smiles_norm=smiles_norm,
                out_dir=mol_output_dir,
                np_files=np_files,
                templates_dir=templates_dir,
                dry_run=dry_run,
            )
            
            # Verify outputs
            if not dry_run:
                missing = [p for p in (outputs.topology, outputs.coords) if not p.exists()]
                if missing:
                    raise RuntimeError(f"Builder reported outputs but files are missing: {missing}")
            
            # Write manifest
            write_manifest(
                out_dir=mol_output_dir,
                smiles_in=smiles,
                smiles_norm=smiles_norm,
                ligand_norm=ligand_norm,
                outputs=outputs,
                git_commit=git_commit,
                status="complete",
            )
            
            results['successful'].append({
                'row': idx,
                'name': name,
                'smiles': smiles,
                'output_dir': str(mol_output_dir.relative_to(output_root)),
            })
            
            print(f"   ✅ Success: {mol_output_dir.name}")
            
        except Exception as e:
            error_msg = str(e)
            print(f"   ❌ Failed: {error_msg}")
            
            # Try to write failure manifest
            try:
                safe_mkdir(mol_output_dir)
                smiles_norm_fallback = smiles
                ligand_norm_fallback = ligand_id_from_smiles(smiles_norm_fallback) if smiles else "unknown"
                
                dummy_outputs = BuildOutputs(
                    topology=mol_output_dir / "topology.prmtop",
                    coords=mol_output_dir / "coords.inpcrd",
                )
                write_manifest(
                    out_dir=mol_output_dir,
                    smiles_in=smiles,
                    smiles_norm=smiles_norm_fallback,
                    ligand_norm=ligand_norm_fallback,
                    outputs=dummy_outputs,
                    git_commit=git_commit,
                    status="failed",
                )
            except Exception:
                pass
            
            results['failed'].append({
                'row': idx,
                'name': name,
                'smiles': smiles,
                'error': error_msg,
            })
        
        print()  # Blank line between molecules
    
    # Print summary
    print("=" * 70)
    print("📊 BATCH PROCESSING SUMMARY")
    print("=" * 70)
    print(f"Total molecules:    {results['total']}")
    print(f"✅ Successful:      {len(results['successful'])}")
    print(f"❌ Failed:          {len(results['failed'])}")
    print(f"⏭️  Skipped:         {len(results['skipped'])}")
    print("=" * 70)
    
    # Write batch summary
    summary = {
        'csv_file': str(csv_path),
        'output_directory': str(output_root),
        'timestamp': iso_now_utc(),
        'git_commit': git_commit,
        'statistics': {
            'total': results['total'],
            'successful': len(results['successful']),
            'failed': len(results['failed']),
            'skipped': len(results['skipped']),
        },
        'successful_builds': results['successful'],
        'failed_builds': results['failed'],
        'skipped_rows': results['skipped'],
    }
    
    summary_path = output_root / "batch_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"\n📝 Batch summary written to: {summary_path}\n")
    
    return summary


def main(argv: list[str]) -> int:
    args = parse_args(argv)

    this_file = pathlib.Path(__file__).resolve()
    stageA_root = this_file.parent.parent
    templates_dir = (
        pathlib.Path(args.templates).resolve()
        if args.templates
        else (this_file.parent / "templates").resolve()
    )

    default_np_files = resolve_default_np_files(stageA_root)
    if default_np_files:
        print(f"📦 Using bundled NP assets from {default_np_files['lib'].parent}")
    else:
        print("ℹ️  No bundled NP assets detected; will build ligand-only systems")

    repo_root = find_repo_root(stageA_root)
    git_commit = get_git_commit(repo_root)

    # ========================================================================
    # CSV BATCH MODE
    # ========================================================================
    if args.csv:
        # Validate output directory
        if not args.output_dir:
            sys.stderr.write("ERROR: --output-dir is required when using --csv mode\n")
            return 1
        
        csv_path = pathlib.Path(args.csv).expanduser().resolve()
        output_root = pathlib.Path(args.output_dir).expanduser().resolve()
        safe_mkdir(output_root)
        
        try:
            summary = process_csv_batch(
                csv_path=csv_path,
                output_root=output_root,
                templates_dir=templates_dir,
                git_commit=git_commit,
                np_files=default_np_files,
                dry_run=args.dry_run,
            )
            
            # Return non-zero if any builds failed
            if summary['statistics']['failed'] > 0:
                return 1
            return 0
            
        except Exception as e:
            sys.stderr.write(f"[stageA_builder] CSV BATCH ERROR: {e}\n")
            return 2

    # ========================================================================
    # SINGLE MOLECULE MODE
    # ========================================================================
    else:
        # Validate output directory
        if not args.output:
            sys.stderr.write("ERROR: --output is required when using --smiles mode\n")
            return 1
        
        out_dir = pathlib.Path(args.output).expanduser().resolve()
        safe_mkdir(out_dir)

        smiles_in = args.smiles
        try:
            smiles_norm = normalize_smiles(smiles_in)
            ligand_norm = ligand_id_from_smiles(smiles_norm)

            outputs = run_builder(
                smiles_norm=smiles_norm,
                out_dir=out_dir,
                np_files=default_np_files,
                templates_dir=templates_dir,
                dry_run=args.dry_run,
            )

            if not args.dry_run:
                missing = [p for p in (outputs.topology, outputs.coords) if not p.exists()]
                if missing:
                    raise RuntimeError(f"Builder reported outputs but files are missing: {missing}")

            write_manifest(
                out_dir=out_dir,
                smiles_in=smiles_in,
                smiles_norm=smiles_norm,
                ligand_norm=ligand_norm,
                outputs=outputs,
                git_commit=git_commit,
                status="complete",
            )

            return 0

        except Exception as e:
            try:
                smiles_norm_fallback = smiles_in.strip()
                ligand_norm_fallback = (
                    ligand_id_from_smiles(smiles_norm_fallback)
                    if smiles_norm_fallback
                    else "unknown"
                )

                dummy_outputs = BuildOutputs(
                    topology=out_dir / "topology.prmtop",
                    coords=out_dir / "coords.inpcrd",
                )
                write_manifest(
                    out_dir=out_dir,
                    smiles_in=smiles_in,
                    smiles_norm=smiles_norm_fallback,
                    ligand_norm=ligand_norm_fallback,
                    outputs=dummy_outputs,
                    git_commit=git_commit,
                    status="failed",
                )
            except Exception:
                pass

            sys.stderr.write(f"[stageA_builder] ERROR: {e}\n")
            return 2


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
