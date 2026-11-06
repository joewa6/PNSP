#!/usr/bin/env python3
"""Render final SI Figures S14–S16 from frozen processed data; no refitting."""
from common import DATA, OUTPUT
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

CUTOFFS = (3.5, 4.0, 4.5, 5.0)
CANONICAL_CUTOFF = 4.0
RNG_SEED = 20261006
CHARGE_GROUP_ORDER = ("Neutral", "Anionic", "Cationic")

def make_figure(
    systems: pd.DataFrame,
    rank_summary: pd.DataFrame,
    charge_summary: pd.DataFrame,
    mode_contrasts: pd.DataFrame,
) -> None:
    plt.rcParams.update({
        "font.size": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.dpi": 150,
    })
    fig, axes = plt.subplots(2, 2, figsize=(10.2, 7.6))
    ax_a, ax_b, ax_c, ax_d = axes.ravel()

    rng = np.random.default_rng(RNG_SEED)
    values_by_cutoff = [systems.loc[systems.cutoff_A.eq(c), "occupancy"].to_numpy() for c in CUTOFFS]
    box = ax_a.boxplot(
        values_by_cutoff, positions=CUTOFFS, widths=0.20, patch_artist=True,
        showfliers=False, medianprops={"color": "black", "linewidth": 1.4},
        whiskerprops={"color": "#555555"}, capprops={"color": "#555555"},
        boxprops={"facecolor": "#d9d9d9", "edgecolor": "#555555"},
    )
    for cutoff, values in zip(CUTOFFS, values_by_cutoff):
        jitter = rng.uniform(-0.055, 0.055, len(values))
        ax_a.scatter(np.full(len(values), cutoff) + jitter, values, s=6, alpha=0.22,
                     color="#555555", edgecolors="none", rasterized=True)
    ax_a.set(xlabel="Contact cutoff (Å)", ylabel="Mean contaminant occupancy",
             xticks=CUTOFFS, ylim=(-0.03, 1.03))

    nonreference = rank_summary.loc[~rank_summary.cutoff_A.eq(CANONICAL_CUTOFF)]
    ax_b.plot(nonreference.cutoff_A, nonreference.spearman_vs_4A, "o-", color="#333333")
    for row in nonreference.itertuples(index=False):
        if row.cutoff_A == 3.5:
            offset, horizontal, vertical = (12, -2), "left", "center"
        else:
            offset, horizontal, vertical = (0, 8), "center", "bottom"
        ax_b.annotate(
            f"ρ={row.spearman_vs_4A:.3f}\nmedian |Δrank|={row.median_absolute_rank_shift:.0f}",
            (row.cutoff_A, row.spearman_vs_4A), xytext=offset,
            textcoords="offset points", ha=horizontal, va=vertical, fontsize=8,
        )
    ax_b.set(xlabel="Alternative cutoff (Å)", ylabel="Spearman correlation with 4.0 Å ranking",
             xticks=nonreference.cutoff_A, ylim=(max(0.8, nonreference.spearman_vs_4A.min() - 0.04), 1.01))

    styles = {
        "Neutral": ("#222222", "o", "-"),
        "Anionic": ("#777777", "s", "--"),
        "Cationic": ("#aaaaaa", "^", ":"),
    }
    for group in CHARGE_GROUP_ORDER:
        color, marker, linestyle = styles[group]
        data = charge_summary.loc[charge_summary.charge_group.eq(group)].sort_values("cutoff_A")
        n_systems = int(data.systems.iloc[0])
        ax_c.errorbar(
            data.cutoff_A, data.mean_occupancy,
            yerr=[data.mean_occupancy - data.ci90_lo, data.ci90_hi - data.mean_occupancy],
            color=color, marker=marker, linestyle=linestyle, linewidth=1.2,
            markersize=4, capsize=2, label=f"{group} (n={n_systems})",
        )
    ax_c.set(xlabel="Contact cutoff (Å)", ylabel="Mean occupancy",
             xticks=CUTOFFS, ylim=(-0.03, 1.03))
    ax_c.legend(frameon=False, fontsize=8, loc="lower right")

    order = ["Rigid polyaromatic", "Flexible neutral", "Cationic", "Anionic"]
    data = mode_contrasts.set_index("operational_group").loc[order].reset_index()
    y = np.arange(len(data))
    ax_d.errorbar(
        data.mean_change_5A_minus_3p5A, y,
        xerr=[data.mean_change_5A_minus_3p5A - data.ci90_lo,
              data.ci90_hi - data.mean_change_5A_minus_3p5A],
        fmt="o", color="#333333", ecolor="#777777", capsize=3,
    )
    labels = [f"{g} (n={n})" for g, n in zip(data.operational_group, data.systems)]
    ax_d.set(xlabel="Mean occupancy change (5.0 Å minus 3.5 Å)",
             yticks=y, yticklabels=labels)
    ax_d.invert_yaxis()

    for label, ax in zip("ABCD", axes.ravel()):
        ax.text(-0.14, 1.05, label, transform=ax.transAxes, fontsize=13,
                fontweight="bold", va="top")
        ax.grid(axis="y", color="#e6e6e6", linewidth=0.7)
    fig.tight_layout()
    fig.savefig(OUTPUT / "Figure_S14_cutoff_sensitivity.png", dpi=400, bbox_inches="tight")
    fig.savefig(OUTPUT / "Figure_S14_cutoff_sensitivity.pdf", bbox_inches="tight")
    plt.close(fig)


FEATURES = [
    "MW", "LogP", "TPSA", "HBD", "HBA", "RotatableBonds",
    "NumAromaticRings", "HeavyAtomCount", "NumHalogens",
]
LABELS = {
    "MW": "MW", "LogP": "logP", "TPSA": "TPSA", "HBD": "HBD",
    "HBA": "HBA", "RotatableBonds": "Rotatable bonds",
    "NumAromaticRings": "Aromatic rings", "HeavyAtomCount": "Heavy atoms",
    "NumHalogens": "Halogens",
}
ABLATION_ORDER = [
    ("Full descriptors", "Full model"),
    ("MW only", "MW only"),
    ("Size-linked trio", "MW + heavy atoms\n+ logP only"),
    ("Full minus MW", "Full − MW"),
    ("Full minus heavy atoms", "Full − heavy atoms"),
    ("Full minus logP", "Full − logP"),
    ("Full minus size-linked trio", "Full − MW − heavy atoms\n− logP"),
]
PERMUTATION_ORDER = [
    "MW", "Heavy-atom count", "logP", "MW + heavy-atom count",
    "Size-linked trio", "Polar block", "All descriptors",
]


def s15() -> None:
    correlations = pd.read_csv(DATA / "neutral_descriptor_spearman_matrix.csv", index_col=0)
    vifs = pd.read_csv(DATA / "descriptor_vif.csv")
    repeated = pd.read_csv(DATA / "neutral_repeated_cv_performance.csv")
    grouped = pd.read_csv(DATA / "neutral_grouped_permutation_importance.csv")

    if correlations.shape != (len(FEATURES), len(FEATURES)):
        raise RuntimeError("Unexpected neutral correlation-matrix dimensions")
    if set(repeated.feature_set) != {name for name, _ in ABLATION_ORDER}:
        raise RuntimeError("Repeated-CV feature sets do not match the expected ablation sets")
    if set(grouped.group) != set(PERMUTATION_ORDER):
        raise RuntimeError("Grouped-permutation sets do not match the expected groups")

    plt.rcParams.update({
        "font.size": 8.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.dpi": 150,
    })
    fig, axes = plt.subplots(2, 2, figsize=(13.4, 10.4), layout="constrained")
    ax_a, ax_b, ax_c, ax_d = axes.ravel()

    values = correlations.loc[FEATURES, FEATURES].to_numpy(float)
    image = ax_a.imshow(values, vmin=-1, vmax=1, cmap="coolwarm")
    ticks = np.arange(len(FEATURES))
    axis_labels = [LABELS[feature] for feature in FEATURES]
    ax_a.set_xticks(ticks, axis_labels, rotation=48, ha="right")
    ax_a.set_yticks(ticks, axis_labels)
    for row in range(len(FEATURES)):
        for column in range(len(FEATURES)):
            if row == column:
                continue
            value = values[row, column]
            ax_a.text(
                column, row, f"{value:.2f}", ha="center", va="center", fontsize=5.5,
                color="white" if abs(value) > 0.65 else "black",
            )
    colorbar = fig.colorbar(image, ax=ax_a, fraction=0.046, pad=0.03)
    colorbar.set_label("Spearman ρ")
    ax_a.set_title("A. Neutral-contaminant correlations (n = 168)",
                   loc="left", fontweight="bold")

    full_vif = vifs.loc[vifs.dataset.eq("full")].set_index("descriptor")
    neutral_vif = vifs.loc[vifs.dataset.eq("neutral")].set_index("descriptor")
    descriptor_order = sorted(FEATURES, key=lambda feature: neutral_vif.loc[feature, "vif"])
    y = np.arange(len(descriptor_order))
    ax_b.barh(
        y - 0.18, [full_vif.loc[feature, "vif"] for feature in descriptor_order],
        height=0.34, color="#a0a0a0", label="All contaminants (n = 203)",
    )
    ax_b.barh(
        y + 0.18, [neutral_vif.loc[feature, "vif"] for feature in descriptor_order],
        height=0.34, color="#4c78a8", label="Neutral subset (n = 168)",
    )
    ax_b.axvline(5, color="#666666", linestyle="--", linewidth=0.9)
    ax_b.axvline(10, color="#222222", linestyle=":", linewidth=0.9)
    ax_b.set_yticks(y, [LABELS[feature] for feature in descriptor_order])
    ax_b.set_xlabel("Variance inflation factor")
    ax_b.set_title("B. Variance inflation", loc="left", fontweight="bold")
    ax_b.legend(frameon=False, fontsize=8, loc="lower right")
    ax_b.grid(axis="x", color="#e5e5e5")

    pivot = repeated.pivot(index="seed", columns="feature_set", values="heldout_r2")
    rng = np.random.default_rng(42)
    x = np.arange(len(ABLATION_ORDER))
    for position, (feature_set, _) in enumerate(ABLATION_ORDER):
        scores = pivot[feature_set].to_numpy(float)
        if len(scores) != 20 or not np.isfinite(scores).all():
            raise RuntimeError(f"Invalid repeated-CV scores for {feature_set}")
        mean = scores.mean()
        sd = scores.std(ddof=1)
        ax_c.bar(
            position, mean, width=0.65,
            color="#59a14f" if feature_set == "Size-linked trio" else "#a6bdd5",
            zorder=2,
        )
        ax_c.scatter(
            np.full(len(scores), position) + rng.uniform(-0.10, 0.10, len(scores)),
            scores, s=10, color="#555555", alpha=0.5, linewidth=0, zorder=3,
        )
        ax_c.errorbar(position, mean, yerr=sd, fmt="none", ecolor="black",
                      capsize=3, zorder=4)
        ax_c.text(
            position, max(scores.max(), mean + sd) + 0.012,
            f"{mean:.3f}\n± {sd:.3f}", ha="center", va="bottom", fontsize=7,
        )
    ax_c.set_xticks(x, [label for _, label in ABLATION_ORDER], rotation=36, ha="right")
    ax_c.set_ylim(0, 0.88)
    ax_c.set_ylabel("Mean out-of-fold GAM R²")
    ax_c.set_title(
        "C. Predictive sensitivity to descriptor selection\n"
        "Mean ± SD across 20 five-fold partitions",
        loc="left", fontweight="bold",
    )
    ax_c.grid(axis="y", color="#e5e5e5")

    permutation = grouped.set_index("group").loc[PERMUTATION_ORDER]
    y = np.arange(len(permutation))
    ax_d.barh(
        y, permutation.mean_delta_r2, xerr=permutation.sd_delta_r2,
        color="#b185ad", error_kw={"ecolor": "#333333", "capsize": 2},
    )
    ax_d.set_yticks(y, permutation.index)
    ax_d.invert_yaxis()
    ax_d.set_xlabel("Decrease in fitted-model R² after joint permutation")
    ax_d.set_title(
        "D. Grouped permutation sensitivity\n"
        "Mean ± SD across 200 permutations",
        loc="left", fontweight="bold",
    )
    ax_d.grid(axis="x", color="#e5e5e5")
    ax_d.text(
        0.98, 0.03, "Baseline fitted R² = 0.831",
        transform=ax_d.transAxes, ha="right", va="bottom", fontsize=8,
    )

    for axis in axes.ravel():
        axis.set_axisbelow(True)

    fig.savefig(OUTPUT / "Figure_S15_descriptor_collinearity_and_sensitivity.png",
                dpi=400, bbox_inches="tight")
    fig.savefig(OUTPUT / "Figure_S15_descriptor_collinearity_and_sensitivity.pdf",
                bbox_inches="tight")
    plt.close(fig)


def s16():
    curves = pd.read_csv(DATA / "MW_response_curves.csv")
    data = pd.read_csv(DATA / "neutral_support.csv")
    low, high = curves.MW.min(), curves.MW.max()
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig = plt.figure(figsize=(12, 8), layout="constrained")
    gs = fig.add_gridspec(2, 3, height_ratios=[1.4, 1])
    ax = fig.add_subplot(gs[0, :])
    for label, group in curves.groupby("model", sort=False):
        ax.plot(group.MW, group.centered_MW_term, label=label, lw=2)
    ax.axvspan(130, 150, color="grey", alpha=.14)
    ax.axhline(0, color="grey", lw=.7)
    ax.set(xlim=(low, high), xlabel="Molecular weight (Da)", ylabel="Centered MW contribution to occupancy",
           title="A. MW response under four descriptor specifications")
    ax.legend(frameon=False, fontsize=9)
    for i, (column, ylabel) in enumerate([("occupancy_mean", "Mean occupancy"), ("LogP", "logP"), ("HeavyAtomCount", "Heavy-atom count")]):
        panel = fig.add_subplot(gs[1, i])
        panel.scatter(data.MW, data[column], s=17, alpha=.55, color="#4c78a8")
        panel.axvspan(130, 150, color="grey", alpha=.14)
        panel.axvline(low, color="grey", ls=":", lw=.8)
        panel.axvline(high, color="grey", ls=":", lw=.8)
        panel.set(xlabel="Molecular weight (Da)", ylabel=ylabel,
                  title=f"{'BCD'[i]}. Observed compounds (n = 168)")
    for suffix in ["png", "pdf"]:
        fig.savefig(OUTPUT / f"Figure_S16_MW_response_comparison.{suffix}", dpi=300)
    plt.close(fig)

def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    make_figure(*(pd.read_csv(DATA / name) for name in (
        "cutoff_sensitivity_by_system.csv", "cutoff_sensitivity_rank_stability.csv",
        "cutoff_sensitivity_charge_groups.csv", "cutoff_sensitivity_mode_contrasts.csv")))
    s15()
    s16()

if __name__ == "__main__":
    main()
