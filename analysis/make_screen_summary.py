#!/usr/bin/env python3
"""Combine replica-level copy metrics into the ranked screen summary."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def discover(paths: list[Path]) -> list[Path]:
    files: list[Path] = []
    for path in paths:
        path = path.resolve()
        if path.is_file():
            files.append(path)
        elif path.is_dir():
            files.extend(path.rglob("copy_metrics.csv"))
        else:
            raise FileNotFoundError(path)
    files = sorted(set(files))
    if not files:
        raise FileNotFoundError("No copy_metrics.csv inputs were found")
    return files


def occupancy_bootstrap(values: np.ndarray, n_boot: int, seed: int) -> dict[str, float]:
    rng = np.random.default_rng(seed)
    indices = rng.integers(0, len(values), size=(n_boot, len(values)))
    draws = values[indices].mean(axis=1)
    return {
        "occupancy_bootstrap_mean": float(draws.mean()),
        "occupancy_ci90_low": float(np.percentile(draws, 5)),
        "occupancy_ci90_high": float(np.percentile(draws, 95)),
        "occupancy_lower90_ranking_bound": float(np.percentile(draws, 10)),
    }


def residence_bootstrap(group: pd.DataFrame, n_boot: int, seed: int) -> dict[str, float | int]:
    strata = [
        replica["mean_event_duration_ns"].dropna().to_numpy(float)
        for _, replica in group.groupby("replica", sort=True)
    ]
    if any(len(values) == 0 for values in strata):
        raise RuntimeError(f"{group.system.iloc[0]} has an entire replica with no contact event")
    pooled = np.concatenate(strata)
    rng = np.random.default_rng(seed)
    draws = np.empty(n_boot)
    for index in range(n_boot):
        resampled = [
            values[rng.integers(0, len(values), size=len(values))] for values in strata
        ]
        draws[index] = np.concatenate(resampled).mean()
    return {
        "residence_copy_trials_included": int(len(pooled)),
        "mean_residence_time_ns": float(pooled.mean()),
        "residence_time_sd_ns": float(pooled.std(ddof=1)) if len(pooled) > 1 else np.nan,
        "residence_time_ci90_low_ns": float(np.percentile(draws, 5)),
        "residence_time_ci90_high_ns": float(np.percentile(draws, 95)),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--output", type=Path, default=Path("screen_summary"))
    parser.add_argument("--expected-replicas", type=int, default=3)
    parser.add_argument("--expected-copies", type=int, default=10)
    parser.add_argument("--occupancy-bootstrap", type=int, default=2000)
    parser.add_argument("--residence-bootstrap", type=int, default=10000)
    args = parser.parse_args()

    files = discover(args.inputs)
    copies = pd.concat([pd.read_csv(path) for path in files], ignore_index=True)
    required = {
        "system", "replica", "copy_id", "occupancy", "n_events",
        "mean_event_duration_ns", "maximum_event_duration_ns",
    }
    missing = required - set(copies.columns)
    if missing:
        raise ValueError(f"Copy metrics are missing columns: {sorted(missing)}")
    keys = ["system", "replica", "copy_id"]
    if copies.duplicated(keys).any():
        raise ValueError("Duplicate system/replica/copy rows were found")

    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    for system_index, (system, group) in enumerate(copies.groupby("system", sort=True)):
        replicas = group["replica"].nunique()
        copy_counts = group.groupby("replica")["copy_id"].nunique()
        if replicas != args.expected_replicas or not copy_counts.eq(args.expected_copies).all():
            raise RuntimeError(
                f"{system}: expected {args.expected_replicas} replicas x "
                f"{args.expected_copies} copies; found {replicas} replicas and "
                f"per-replica counts {copy_counts.to_dict()}"
            )
        occupancy = group["occupancy"].to_numpy(float)
        event_positive = group.loc[group["n_events"] > 0]
        row: dict[str, object] = {
            "system": system,
            "replicas": replicas,
            "copy_trials": len(group),
            "occupancy_mean": float(occupancy.mean()),
            "occupancy_sd": float(occupancy.std(ddof=1)),
            "occupancy_median": float(np.median(occupancy)),
            "total_contact_events": int(group["n_events"].sum()),
            "mean_maximum_event_duration_ns": float(
                event_positive["maximum_event_duration_ns"].mean()
            ),
            **occupancy_bootstrap(
                occupancy, args.occupancy_bootstrap, seed=42 + system_index
            ),
            **residence_bootstrap(group, args.residence_bootstrap, seed=42),
        }
        rows.append(row)

    summary = pd.DataFrame(rows).sort_values(
        ["occupancy_mean", "system"], ascending=[False, True]
    ).reset_index(drop=True)
    summary.insert(0, "rank", np.arange(1, len(summary) + 1))
    summary.to_csv(output / "ranked_screen_summary.csv", index=False)
    copies.sort_values(keys).to_csv(output / "all_copy_metrics.csv", index=False)

    figure, axis = plt.subplots(figsize=(7.0, 3.2))
    x = summary["rank"].to_numpy()
    y = summary["occupancy_mean"].to_numpy()
    low = summary["occupancy_ci90_low"].to_numpy()
    high = summary["occupancy_ci90_high"].to_numpy()
    axis.errorbar(
        x, y, yerr=[y - low, high - y], fmt="o", markersize=2.2,
        color="#2a5f8a", ecolor="#8bb9de", elinewidth=0.5, capsize=0,
    )
    axis.plot(
        x, summary["occupancy_lower90_ranking_bound"], color="0.55", linewidth=1.1
    )
    axis.set(
        xlabel="Compound rank (by mean occupancy)", ylabel="Occupancy",
        xlim=(0, len(summary) + 1), ylim=(-0.02, 1.03),
    )
    axis.grid(axis="y", alpha=0.2)
    figure.tight_layout()
    figure.savefig(output / "ranked_screen_occupancy.png", dpi=300)
    figure.savefig(output / "ranked_screen_occupancy.pdf")
    plt.close(figure)

    metadata = {
        "systems": len(summary),
        "copy_records": len(copies),
        "occupancy_definition": "mean of copy occupancies; each copy occupancy is bound frames / total frames",
        "occupancy_uncertainty": f"pooled copy bootstrap; {args.occupancy_bootstrap} resamples; 5th-95th percentile interval",
        "residence_definition": "mean event duration per event-positive copy, then arithmetic mean across all event-positive copies from all replicas",
        "residence_exclusion": "copy trajectories with no contact event",
        "residence_uncertainty": f"copy bootstrap within fixed replica strata; {args.residence_bootstrap} resamples; 5th-95th percentile interval",
        "input_files": [str(path) for path in files],
    }
    (output / "summary_manifest.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"[complete] {len(summary)} systems -> {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
