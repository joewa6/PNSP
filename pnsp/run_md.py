#!/usr/bin/env python3
"""Prepare or execute the staged AMBER protocol used for PET screening."""

from __future__ import annotations

import argparse
import hashlib
import json
import shlex
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

if __package__:
    from .provenance import execution_provenance
else:
    from provenance import execution_provenance


PET_HEAVY_MASK = ":1-250&!@H="


@dataclass(frozen=True)
class Stage:
    name: str
    coordinate: str
    restart: str
    trajectory: str
    output: str
    reference: str
    mdin: str


def cntrl(title: str, lines: list[str]) -> str:
    return title + "\n &cntrl\n   " + "\n   ".join(lines) + "\n/\n &ewald\n   ew_type=0, use_pme=1,\n/\n"


def protocol(production_ns: float, temperature: float, pressure: float) -> list[Stage]:
    production_steps = round(production_ns * 1_000_000 / 2.0)
    production_stride = round(10.0 / 0.002)
    restraint = f"ntr=1, restraint_wt=300, restraintmask='{PET_HEAVY_MASK}',"
    mild_restraint = f"ntr=1, restraint_wt=20, restraintmask='{PET_HEAVY_MASK}',"
    stages = [
        Stage("min_solv", "system.inpcrd", "min_solv.rst", "min_solv.nc", "min_solv.out", "system.inpcrd", cntrl("Restrained solvent minimization", [
            "imin=1, ntx=1, irest=0, ntxo=2,", restraint,
            "maxcyc=10000, ncyc=5000, dx0=0.001,", "ntpr=1000, ntwr=1000, ntwx=1000, ioutfm=1,",
            "ntc=1, ntf=1, ntb=1, cut=10.0,",
        ])),
        Stage("solvent_nvt", "min_solv.rst", "solvent_nvt.rst", "solvent_nvt.nc", "solvent_nvt.out", "min_solv.rst", cntrl("Restrained solvent NVT equilibration", [
            "imin=0, ntx=1, irest=0, ntxo=2,", restraint,
            "nstlim=2500000, dt=0.002,", f"ntt=2, tempi={temperature}, temp0={temperature}, ig=-1, vrand=5000,",
            "ntpr=10000, ntwr=100000, ntwx=10000, ioutfm=1,", "ntc=2, ntf=2, ntb=1, cut=10.0,",
        ])),
        Stage("solvent_npt", "solvent_nvt.rst", "solvent_npt.rst", "solvent_npt.nc", "solvent_npt.out", "solvent_nvt.rst", cntrl("Restrained solvent NPT equilibration", [
            "imin=0, ntx=5, irest=1, ntxo=2,", restraint,
            "nstlim=2000000, dt=0.002,", f"ntt=2, temp0={temperature}, ig=-1, vrand=5000,",
            f"ntp=1, barostat=2, pres0={pressure}, taup=1.0, mcbarint=5000,",
            "ntpr=1000, ntwr=100000, ntwx=10000, ioutfm=1,", "ntc=2, ntf=2, ntb=2, cut=10.0,",
        ])),
        Stage("min_restrained", "solvent_npt.rst", "min_restrained.rst", "min_restrained.nc", "min_restrained.out", "solvent_npt.rst", cntrl("Restrained full-system minimization", [
            "imin=1, ntx=1, irest=0, ntxo=2,", mild_restraint,
            "maxcyc=20000, ncyc=10000,", "ntpr=1000, ntwr=10000, ntwx=1000, ioutfm=1,",
            "ntc=1, ntf=1, ntb=1, cut=10.0,",
        ])),
        Stage("min_unrestrained", "min_restrained.rst", "min_unrestrained.rst", "min_unrestrained.nc", "min_unrestrained.out", "min_restrained.rst", cntrl("Unrestrained full-system minimization", [
            "imin=1, ntx=1, irest=0, ntxo=2,", "maxcyc=20000, ncyc=10000,",
            "ntpr=1000, ntwr=10000, ntwx=1000, ioutfm=1,", "ntc=1, ntf=1, ntb=1, cut=10.0,",
        ])),
    ]

    heating = [
        ("heat_100_nvt", "min_unrestrained.rst", 100, 100000, False),
        ("heat_100_npt", "heat_100_nvt.rst", 100, 100000, True),
        ("heat_200_npt", "heat_100_npt.rst", 200, 100000, True),
        ("heat_300_npt", "heat_200_npt.rst", temperature, 200000, True),
    ]
    for name, coordinate, target, steps, npt in heating:
        pressure_line = (
            f"ntp=1, barostat=2, pres0={pressure}, taup=1.0, mcbarint=1000,"
            if npt else ""
        )
        stages.append(Stage(name, coordinate, f"{name}.rst", f"{name}.nc", f"{name}.out", coordinate, cntrl(name.replace("_", " "), [
            "imin=0, ntx=5, irest=1, ntxo=2," if coordinate != "min_unrestrained.rst" else "imin=0, ntx=1, irest=0, ntxo=2,",
            mild_restraint, f"nstlim={steps}, dt=0.001,", f"ntt=2, temp0={target}, ig=-1, vrand=1000,",
            pressure_line, "ntpr=10000, ntwr=100000, ntwx=10000, ioutfm=1,",
            f"ntc=2, ntf=2, ntb={2 if npt else 1}, cut=10.0,",
        ])))
    stages.extend([
        Stage("equil", "heat_300_npt.rst", "equil.rst", "equil.nc", "equil.out", "heat_300_npt.rst", cntrl("Unrestrained NPT equilibration", [
            "imin=0, ntx=5, irest=1, ntxo=2,", "nstlim=1500000, dt=0.001,",
            f"ntt=2, temp0={temperature}, ig=-1, vrand=10000,", f"ntp=1, barostat=2, pres0={pressure}, taup=1.0, mcbarint=10000,",
            "ntpr=10000, ntwr=100000, ntwx=10000, ioutfm=1,", "ntc=2, ntf=2, ntb=2, cut=10.0,",
        ])),
        Stage("prod", "equil.rst", "prod.rst", "prod.nc", "prod.out", "equil.rst", cntrl("Production NPT", [
            "imin=0, ntx=5, irest=1, ntxo=2,", f"nstlim={production_steps}, dt=0.002,",
            f"ntt=2, temp0={temperature}, ig=-1, vrand=100000,", f"ntp=1, barostat=2, pres0={pressure}, taup=1.0, mcbarint=100000,",
            f"ntpr=5000, ntwr=100000, ntwx={production_stride}, ioutfm=1,", "ntc=2, ntf=2, ntb=2, cut=10.0,",
        ])),
    ])
    return stages


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def command(engine: str, stage: Stage, mdin_dir: Path) -> list[str]:
    result = [
        engine, "-O", "-p", "system.prmtop", "-i", str(mdin_dir / f"{stage.name}.in"),
        "-c", stage.coordinate, "-r", stage.restart, "-x", stage.trajectory,
        "-o", stage.output,
    ]
    if "restraintmask" in stage.mdin:
        result.extend(["-ref", stage.reference])
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("system", type=Path, help="Directory containing system.prmtop and system.inpcrd")
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--engine", default="pmemd.cuda_SPFP")
    parser.add_argument("--production-ns", type=float, default=15.0)
    parser.add_argument("--temperature", type=float, default=303.0)
    parser.add_argument("--pressure", type=float, default=1.0)
    parser.add_argument("--execute", action="store_true", help="Execute after writing inputs; otherwise prepare only")
    parser.add_argument("--resume", action="store_true", help="Skip stages with nonempty restart and output files")
    args = parser.parse_args()

    source = args.system.resolve()
    output = (args.output or source).resolve()
    source_manifest = source / "system_manifest.json"
    system_provenance = (
        json.loads(source_manifest.read_text()).get("provenance", {})
        if source_manifest.is_file() else {}
    )
    for filename in ("system.prmtop", "system.inpcrd"):
        if not (source / filename).is_file():
            raise FileNotFoundError(source / filename)
    output.mkdir(parents=True, exist_ok=True)
    for filename in ("system.prmtop", "system.inpcrd"):
        destination = output / filename
        if source / filename != destination:
            shutil.copy2(source / filename, destination)

    stages = protocol(args.production_ns, args.temperature, args.pressure)
    mdin_dir = output / "MDIN"
    mdin_dir.mkdir(exist_ok=True)
    for stage in stages:
        (mdin_dir / f"{stage.name}.in").write_text(stage.mdin, encoding="utf-8")
    commands = [command(args.engine, stage, mdin_dir) for stage in stages]
    shell_lines = ["#!/usr/bin/env bash", "set -euo pipefail",
                   'cd -- "$(dirname -- "$0")"', ""]
    shell_lines.extend(" ".join(shlex.quote(part) for part in item) for item in commands)
    run_script = output / "run_md.sh"
    run_script.write_text("\n".join(shell_lines) + "\n", encoding="utf-8")
    run_script.chmod(0o755)

    completed: list[str] = []
    if args.execute:
        engine_path = shutil.which(args.engine)
        if engine_path is None:
            raise FileNotFoundError(f"AMBER engine is not on PATH: {args.engine}")
        for stage, item in zip(stages, commands):
            restart_path = output / stage.restart
            stage_output_path = output / stage.output
            if (
                args.resume
                and restart_path.is_file()
                and stage_output_path.is_file()
                and restart_path.stat().st_size > 0
                and stage_output_path.stat().st_size > 0
            ):
                print(f"[skip] {stage.name}")
                completed.append(stage.name)
                continue
            print(f"[run] {stage.name}", flush=True)
            subprocess.run(item, cwd=output, check=True)
            completed.append(stage.name)

    manifest = {
        "provenance": execution_provenance(
            "run_md", system_provenance.get("contaminant_id"),
            system_provenance.get("replica_id"),
            [source / "system.prmtop", source / "system.inpcrd"], output,
        ),
        "status": "complete" if args.execute and len(completed) == len(stages) else "prepared",
        "engine": args.engine,
        "production_ns": args.production_ns,
        "production_timestep_fs": 2.0,
        "production_output_interval_ps": 10.0,
        "production_expected_frames": round(args.production_ns * 100),
        "temperature_K": args.temperature,
        "pressure_atm": args.pressure,
        "restrained_selection": PET_HEAVY_MASK,
        "stages": [stage.name for stage in stages],
        "mdin_sha256": {
            stage.name: sha256(mdin_dir / f"{stage.name}.in") for stage in stages
        },
    }
    (output / "md_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"[{manifest['status']}] inputs in {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
