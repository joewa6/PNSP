#!/usr/bin/env python3
"""Fail closed if released data or regenerated paper outputs drift."""
from __future__ import annotations

import json
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from common import DATA, OUTPUT, PAPER, REPOSITORY


def close(label: str, observed: float, expected: float, atol: float = 5e-6) -> None:
    if not np.isclose(observed, expected, atol=atol, rtol=0):
        raise RuntimeError(f"{label}: observed {observed!r}, expected {expected!r}")


def main() -> None:
    for row in pd.read_csv(PAPER / "manifest.csv").itertuples():
        path = REPOSITORY / row.public_file
        if hashlib.sha256(path.read_bytes()).hexdigest() != row.public_sha256:
            raise RuntimeError(f"Checksum mismatch: {row.public_file}")
        frame = pd.read_csv(path)
        if frame.shape != (row.rows, row.columns):
            raise RuntimeError(f"Dimension mismatch: {row.public_file}")
    long = pd.read_csv(DATA / "convergence_long_summary.csv").set_index("system")
    cumulative = pd.read_csv(DATA / "convergence_long_cumulative.csv")
    if len(long) != 20 or not long.index.is_unique or len(cumulative) != 800:
        raise RuntimeError("Expected twenty systems and 800 cumulative observations")
    wide = cumulative.pivot(index="system", columns="end_ns", values="occupancy")
    for time in (15, 100, 200):
        np.testing.assert_allclose(wide[time], long.loc[wide.index, f"occupancy_{time}ns"])
    close("20-system 15/200 ns Spearman", spearmanr(long.occupancy_15ns, long.occupancy_200ns).statistic, .947368421)
    close("20-system 15/200 ns MAE", (long.occupancy_15ns-long.occupancy_200ns).abs().mean(), .095666167)
    close("20-system 100/200 ns Spearman", spearmanr(long.occupancy_100ns, long.occupancy_200ns).statistic, .992481203)
    close("20-system 100/200 ns MAE", (long.occupancy_100ns-long.occupancy_200ns).abs().mean(), .019575)
    screen = pd.read_csv(DATA/"screen_summary.csv")
    copies = pd.read_csv(DATA/"copy_metrics.csv")
    events = pd.read_csv(DATA/"contact_events.csv")
    if (len(screen), len(copies), len(events)) != (203, 6090, 48334):
        raise RuntimeError("Core dimensions must be 203 systems, 6,090 copies, and 48,334 events")

    close("occupancy minimum", screen.occupancy_mean.min(), .00637778)
    close("occupancy maximum", screen.occupancy_mean.max(), .99377778)
    close("occupancy mean", screen.occupancy_mean.mean(), .69980)
    close("occupancy SD", screen.occupancy_mean.std(ddof=1), .25593)
    close("occupancy median", screen.occupancy_mean.median(), .79580)
    width = screen.occupancy_ci90_hi-screen.occupancy_ci90_lo
    close("median 90% CI width", width.median(), .15453)
    close("minimum mean residence time", screen.residence_time_mean_ns.min(), .02227)
    close("maximum mean residence time", screen.residence_time_mean_ns.max(), 13.58821)

    decomposition = pd.read_csv(DATA/"event_decomposition_system.csv")
    for column, expected in [("total_events", -.38584), ("mean_residence_time_ns", .88198),
                             ("mean_maximum_residence_time_ns", .93306)]:
        close(f"occupancy/{column} Spearman rho",
              spearmanr(decomposition.occupancy, decomposition[column]).statistic, expected)
    close("median hydration loss", pd.read_csv(DATA/"hydration.csv").waters_lost.median(), 3.0163)

    descriptors = pd.read_csv(DATA/"compound_descriptors.csv")
    organochlorines = descriptors.loc[descriptors.chemical_class.eq("Organochlorine pesticides"), "occupancy_mean"]
    if len(organochlorines) != 10:
        raise RuntimeError("Expected ten organochlorine pesticides")
    close("organochlorine mean occupancy", organochlorines.mean(), .93546)
    close("organochlorine occupancy SD", organochlorines.std(ddof=1), .03236)

    importance = pd.read_csv(DATA/"model_all_importance.csv").set_index("feature")
    close("full GAM net-charge importance", importance.loc["NetCharge", "importance"], .67054)
    neutral_importance = pd.read_csv(DATA/"model_neutral_importance.csv").set_index("feature")
    for feature, expected in [("MW", .39887545), ("LogP", .34649315),
                              ("HeavyAtomCount", .31118433),
                              ("NumAromaticRings", .00351937)]:
        close(f"neutral GAM {feature} importance", neutral_importance.loc[feature, "importance"], expected)

    transitions = pd.read_csv(DATA/"mw_transition_summary.csv").set_index("model")
    close("unadjusted MW transition", transitions.loc["unadjusted", "best_transition_Da"], 132)
    close("class-adjusted MW transition", transitions.loc["class-adjusted", "best_transition_Da"], 134)
    multicopy = pd.read_csv(DATA/"multicopy_independence.csv")
    close("multicopy median ICC", multicopy.ICC.median(), -.00706516)
    if int(multicopy.perm_p_bh_reject.sum()) != 0:
        raise RuntimeError("Expected 0/203 BH-adjusted multicopy hits")

    repeated = pd.read_csv(DATA / "neutral_repeated_cv_performance.csv")
    trio = repeated.loc[repeated.feature_set.eq("Size-linked trio")]
    if len(trio) != 20:
        raise RuntimeError("Expected 20 matched CV partitions")
    close("reduced GAM mean CV R2", trio.heldout_r2.mean(), .744701, atol=5e-6)
    close("reduced GAM SD CV R2", trio.heldout_r2.std(ddof=1), .008713, atol=5e-6)

    figures = ["Figure_3_descriptor_distributions", "Figure_4_rank_resolved_occupancy",
               "Figure_5_representative_association_modes", "Figure_6_chemical_classes",
               "Figure_7_model_performance", "Figure_8_neutral_descriptor_model"]
    figures += [f"Figure_S{i}_{name}" for i, name in {
        2:"occupancy_uncertainty", 3:"full_screen_convergence", 4:"long_trajectory_convergence",
        5:"loading_control", 6:"occupancy_event_decomposition", 7:"hydration",
        8:"binary_mixtures", 9:"descriptor_importance", 10:"surface_charge_control",
        11:"MW_class_robustness", 12:"RESP_AM1BCC", 13:"multicopy_independence",
        1:"PET_structural_characterization", 14:"cutoff_sensitivity",
        15:"descriptor_collinearity_and_sensitivity", 16:"MW_response_comparison"}.items()]
    missing = [str(OUTPUT/f"{stem}.{suffix}") for stem in figures for suffix in ("png", "pdf")
               if not (OUTPUT/f"{stem}.{suffix}").is_file()]
    missing += [str(OUTPUT/name) for name in ("Table_S1_GAM_permutation_importance.csv",
                                              "Table_S2_contaminant_ranking.csv")
                if not (OUTPUT/name).is_file()]
    if missing:
        raise RuntimeError("Missing regenerated outputs:\n"+"\n".join(missing))

    report = {"status": "PASS", "systems": 203, "copy_trajectories": 6090,
              "contact_events": 48334, "figures": len(figures), "tables": 2,
              "organochlorine_mean_sd": [float(organochlorines.mean()), float(organochlorines.std(ddof=1))]}
    (OUTPUT/"release_verification.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
