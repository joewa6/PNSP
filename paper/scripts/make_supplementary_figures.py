#!/usr/bin/env python3
"""Regenerate submitted numerical Supporting-Information Figures S2-S13.

This is only a rendering layer: calculations are frozen in processed_data.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from common import DATA, configure_style, save

BLUE, GREEN, RED = "#4C78A8", "#59A14F", "#E15759"
ORANGE, PURPLE = "#F28E2B", "#B279A2"


def clean(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(color="#DDDDDD", linewidth=.45, alpha=.6)


def compare(ax, x, y, xl, yl, color):
    rho = spearmanr(x, y).statistic
    mae = np.mean(np.abs(np.asarray(x)-np.asarray(y)))
    ax.scatter(x, y, s=11, alpha=.38, color=color, edgecolors="none")
    ax.plot([0, 1], [0, 1], "--", color="#555555", lw=.9)
    ax.text(.04, .95, f"Spearman ρ = {rho:.3f}\nMAE = {mae:.3f}",
            transform=ax.transAxes, va="top", fontsize=8)
    ax.set(xlim=(-.02, 1.02), ylim=(-.02, 1.02), xlabel=xl, ylabel=yl)
    clean(ax)


def s2():
    summary = pd.read_csv(DATA/"screen_summary.csv")
    copies = pd.read_csv(DATA/"copy_metrics.csv")
    summary["width"] = summary.occupancy_ci90_hi-summary.occupancy_ci90_lo
    selected = ["Bis_2-ethylhexyl_hexanedioate", "Carboxin", "Cyanazine",
                "Perfluorooctanoic_acid", "Benzene", "Formaldehyde"]
    labels = ["Bis(2-ethylhexyl)\nhexanedioate", "Carboxin", "Cyanazine",
              "Perfluorooctanoic acid", "Benzene", "Formaldehyde"]
    fig = plt.figure(figsize=(11.2, 8.2))
    grid = fig.add_gridspec(2, 1, height_ratios=(.82, 1.35), hspace=.42)
    ax = fig.add_subplot(grid[0])
    points = ax.scatter(summary["rank"], summary.width, c=summary.occupancy_mean,
                        cmap="viridis", vmin=0, vmax=1, s=24, alpha=.78,
                        edgecolors="none")
    ax.axhline(summary.width.median(), color="black", ls="--", lw=1)
    ax.text(.02, .94, f"Median = {summary.width.median():.3f}",
            transform=ax.transAxes, va="top")
    ax.set(xlabel="Contaminant rank by mean occupancy", ylabel="90% bootstrap CI width",
           xlim=(0, 204), ylim=(0, max(.30, summary.width.max()*1.06)))
    bar = fig.colorbar(points, ax=ax, pad=.015, fraction=.025)
    bar.set_label("Mean occupancy")
    ax.text(-.055, 1.03, "A", transform=ax.transAxes, fontweight="bold", fontsize=12)
    clean(ax)
    ax = fig.add_subplot(grid[1]); rng = np.random.default_rng(204)
    replica_colors = {1: BLUE, 2: "#F58518", 3: "#54A24B"}
    for i, system in enumerate(selected):
        values = copies.loc[copies.system.eq(system), "occupancy"].to_numpy()
        violin = ax.violinplot(values, positions=[i], widths=.72, showextrema=False)
        for body in violin["bodies"]:
            body.set_facecolor("#C6D7E8"); body.set_edgecolor("#6889AA"); body.set_alpha(.42)
        ax.boxplot(values, positions=[i], widths=.23, showfliers=False, patch_artist=True,
                   boxprops={"facecolor":"white", "edgecolor":"#333333"},
                   medianprops={"color":"#333333", "linewidth":1.4})
        group = copies.loc[copies.system.eq(system)].sort_values(["replica", "copy"])
        for replica in (1, 2, 3):
            rv = group.loc[group.replica.eq(f"rep{replica}"), "occupancy"].to_numpy()
            ax.scatter(i+rng.normal(0, .055, len(rv)), rv, s=20, alpha=.72,
                       color=replica_colors[replica], edgecolor="white", linewidth=.3,
                       label=f"Replica {replica}" if i == 0 else None, zorder=3)
        row = summary.loc[summary.system.eq(system)].iloc[0]
        ax.errorbar(i+.30, row.occupancy_mean,
                    yerr=[[row.occupancy_mean-row.occupancy_ci90_lo],
                          [row.occupancy_ci90_hi-row.occupancy_mean]],
                    fmt="D", color="black", markersize=4, capsize=2.5, zorder=5)
    ax.set(xticks=range(6), xticklabels=labels, ylabel="Copy-level occupancy",
           ylim=(-.03, 1.07), xlim=(-.55, 5.55))
    ax.legend(frameon=False, ncol=3, loc="lower center", bbox_to_anchor=(.5, -.23))
    ax.text(-.055, 1.03, "B", transform=ax.transAxes, fontweight="bold", fontsize=12)
    clean(ax)
    fig.subplots_adjust(left=.09, right=.96, bottom=.12, top=.98)
    save(fig, "Figure_S2_occupancy_uncertainty")


def s3():
    d = pd.read_csv(DATA/"convergence_full.csv")
    w = d.pivot(index="system", columns="window", values="occupancy")
    fig, axes = plt.subplots(1, 3, figsize=(9, 3.05))
    compare(axes[0], w.block1, w.block2, "Occupancy, 0–5 ns block", "Occupancy, 5–10 ns block", BLUE)
    compare(axes[1], w.block2, w.block3, "Occupancy, 5–10 ns block", "Occupancy, 10–15 ns block", GREEN)
    axes[2].scatter(w.cumulative5, w.cumulative15, s=11, alpha=.30, color=RED, label="5 ns")
    axes[2].scatter(w.cumulative10, w.cumulative15, s=11, alpha=.30, color=BLUE, label="10 ns")
    axes[2].plot([0, 1], [0, 1], "--", color="#555555", lw=.9)
    r5 = spearmanr(w.cumulative5, w.cumulative15).statistic
    r10 = spearmanr(w.cumulative10, w.cumulative15).statistic
    axes[2].text(.04, .95, f"ρ(5,15) = {r5:.3f}\nρ(10,15) = {r10:.3f}",
                 transform=axes[2].transAxes, va="top", fontsize=8)
    axes[2].set(xlim=(-.02, 1.02), ylim=(-.02, 1.02), xlabel="Cumulative occupancy",
                ylabel="Final 15 ns occupancy")
    axes[2].legend(frameon=False, loc="lower right", fontsize=7)
    clean(axes[2])
    for label, ax in zip("ABC", axes):
        ax.text(.02, .98, label, transform=ax.transAxes, va="top", fontweight="bold", fontsize=11)
    fig.subplots_adjust(left=.075, right=.995, bottom=.18, top=.97, wspace=.34)
    save(fig, "Figure_S3_full_screen_convergence")


def s4():
    d = pd.read_csv(DATA/"convergence_long_cumulative.csv")
    order = ["Ammonia", "Bisphenol_A", "Chloroform", "Ethylene_glycol", "Phenol", "Pyrene", "Toluene"]
    colors = dict(zip(order, plt.cm.tab10(np.linspace(0, .9, 7))))
    fig = plt.figure(figsize=(9, 5.8)); grid = fig.add_gridspec(2, 2, height_ratios=(1.25, 1))
    axes = [fig.add_subplot(grid[0, :]), fig.add_subplot(grid[1, 0]), fig.add_subplot(grid[1, 1])]
    for system in order:
        g = d.loc[d.system.eq(system)]
        axes[0].plot(g.end_ns, g.occupancy, lw=1.55, color=colors[system], label=system.replace("_", " "))
    axes[0].axvline(15, color="#333333", ls="--", lw=.9); axes[0].axvline(100, color="#888888", ls=":", lw=.9)
    axes[0].set(xlim=(5, 200), ylim=(-.02, 1.02), xlabel="Cumulative production time (ns)", ylabel="Mean occupancy")
    axes[0].legend(frameon=False, fontsize=6.8, ncol=2, loc="lower right"); clean(axes[0])
    w = d.pivot(index="system", columns="end_ns", values="occupancy")
    compare(axes[1], w[15], w[200], "Occupancy at 15 ns", "Occupancy at 200 ns", BLUE)
    compare(axes[2], w[100], w[200], "Occupancy at 100 ns", "Occupancy at 200 ns", GREEN)
    for label, ax in zip("ABC", axes):
        ax.text(.02, .98, label, transform=ax.transAxes, va="top", fontweight="bold", fontsize=11)
    fig.subplots_adjust(left=.085, right=.99, bottom=.09, top=.98, wspace=.25, hspace=.34)
    save(fig, "Figure_S4_long_trajectory_convergence")


def stacked(file, xcol, levels, order, display, colors, xlabel, stem, negative=False):
    d = pd.read_csv(DATA/file)
    fig, axes = plt.subplots(len(order), 1, figsize=(6.2, 8), sharex=True, sharey=True)
    x = np.arange(len(levels), dtype=float)
    for ax, compound in zip(axes, order):
        g = d.loc[d.compound.eq(compound)].set_index(xcol).loc[levels].reset_index(); color = colors[compound]
        ax.errorbar(x, g.mean_occupancy,
                    yerr=np.vstack((g.mean_occupancy-g.ci90_low, g.ci90_high-g.mean_occupancy)),
                    color=color, marker="o", markeredgecolor="black", markeredgewidth=.5,
                    markersize=6.2, linewidth=1.8, capsize=3)
        ax.set(ylim=(-.04, 1.04), xlim=(-.34, len(levels)-.66), yticks=(0, .5, 1))
        low = compound in {"Pyrene", "Toluene", "Perfluorooctanoic_acid"}
        ax.text(.02, .12 if low else .83, display.get(compound, compound.replace("_", " ")),
                transform=ax.transAxes, color=color, fontweight="bold", va="bottom" if low else "top")
        clean(ax)
    axes[len(axes)//2].set_ylabel("Mean per-copy occupancy")
    axes[-1].set_xticks(x, [f"−{v}" if negative else str(v) for v in levels]); axes[-1].set_xlabel(xlabel)
    fig.subplots_adjust(left=.16, right=.98, bottom=.08, top=.99, hspace=.10); save(fig, stem)


def s5():
    order = ["Pyrene", "Toluene", "Perfluorooctanoic_acid", "Ethylene_glycol", "Nitrate"]
    display = {"Perfluorooctanoic_acid": "Perfluorooctanoate", "Ethylene_glycol": "Ethylene glycol"}
    stacked("loading_control_summary.csv", "loading_copies", [5, 10, 20], order, display,
            dict(zip(order, [BLUE, ORANGE, RED, GREEN, PURPLE])), "Contaminant copies", "Figure_S5_loading_control")


def s6():
    d = pd.read_csv(DATA/"event_decomposition_system.csv")
    specs = [("total_events", "Total number of copy events"), ("mean_residence_time_ns", "Mean residence time (ns)"),
             ("mean_maximum_residence_time_ns", "Mean maximum residence time (ns)")]
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.6))
    for label, ax, (column, ylabel) in zip("ABC", axes, specs):
        g = d[["occupancy", column]].dropna(); rho = spearmanr(g.occupancy, g[column]).statistic
        ax.scatter(g.occupancy, g[column], s=22, alpha=.64, color="#7189A6", edgecolors="none")
        ax.text(.04, .95, f"Spearman ρ = {rho:.3f}\n$n$ = {len(g)}", transform=ax.transAxes, va="top")
        ax.set(xlabel="Mean contaminant occupancy", ylabel=ylabel)
        ax.text(-.13, 1.02, label, transform=ax.transAxes, fontweight="bold", fontsize=12); clean(ax)
    fig.tight_layout(); save(fig, "Figure_S6_occupancy_event_decomposition")


def s7():
    d = pd.read_csv(DATA/"hydration.csv"); fig, axes = plt.subplots(1, 2, figsize=(7.8, 3.8)); rng = np.random.default_rng(42)
    for row in d.itertuples(index=False):
        axes[0].plot([0, 1], [row.hydration_unbound, row.hydration_bound], color="#C9C9C9", lw=.45, alpha=.55)
    axes[0].scatter(rng.normal(0, .035, len(d)), d.hydration_unbound, s=9, alpha=.24, color="#6BAED6", edgecolors="none")
    axes[0].scatter(rng.normal(1, .035, len(d)), d.hydration_bound, s=9, alpha=.24, color="#FB6A4A", edgecolors="none")
    axes[0].plot([0, 1], [d.hydration_unbound.median(), d.hydration_bound.median()], "o-", color="black", lw=1.5)
    axes[0].set(xticks=[0, 1], xticklabels=["Unbound", "Bound"], ylabel="First-shell waters")
    axes[1].hist(d.waters_lost, bins=24, color="#74C476", edgecolor="white")
    axes[1].axvline(d.waters_lost.median(), color="black", ls="--"); axes[1].axvline(0, color="#777777", ls=":")
    axes[1].text(.98, .95, f"Median = {d.waters_lost.median():.2f}", transform=axes[1].transAxes, ha="right", va="top")
    axes[1].set(xlabel="First-shell waters lost on binding", ylabel="Systems")
    for label, ax in zip("AB", axes):
        ax.text(-.12, 1.02, label, transform=ax.transAxes, fontweight="bold", fontsize=12); clean(ax)
    fig.tight_layout(); save(fig, "Figure_S7_hydration")


def s8():
    copies = pd.read_csv(DATA/"mixture_control_copies.csv"); summary = pd.read_csv(DATA/"mixture_control_summary.csv")
    pairs = [("similar", "Bisphenol A + Pyrene"), ("different", "Ethylene glycol + Pyrene")]
    colors = {"Bisphenol_A": BLUE, "Ethylene_glycol": GREEN, "Pyrene": RED}
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.6), sharey=True); rng = np.random.default_rng(51)
    for ax, (role, title) in zip(axes, pairs):
        block = summary.loc[summary.comparison_role.eq(role)]; species = sorted(block.species.unique())
        for li, loading in enumerate([5, 10]):
            for si, name in enumerate(species):
                row = block.loc[(block.copies_per_species.eq(loading)) & block.species.eq(name)].iloc[0]
                vals = copies.loc[(copies.condition.eq(row.condition)) & copies.species.eq(name), "occupancy"].to_numpy()
                x = li+(si-(len(species)-1)/2)*.20
                ax.scatter(x+rng.uniform(-.035, .035, len(vals)), vals, s=13, alpha=.35, color=colors[name], edgecolors="none")
                ax.errorbar(x, row.mean_occupancy, yerr=[[row.mean_occupancy-row.ci90_low], [row.ci90_high-row.mean_occupancy]],
                            fmt="o", color=colors[name], markeredgecolor="black", markeredgewidth=.4, capsize=3,
                            label=name.replace("_", " ") if li == 0 else None)
        ax.set(xticks=[0, 1], xticklabels=["5 + 5", "10 + 10"], xlabel="Copies of each species", title=title, ylim=(-.03, 1.03))
        ax.legend(frameon=False, fontsize=8); clean(ax)
    axes[0].set_ylabel("Per-copy occupancy")
    for label, ax in zip("AB", axes):
        ax.text(-.10, 1.02, label, transform=ax.transAxes, fontweight="bold", fontsize=12)
    fig.tight_layout(); save(fig, "Figure_S8_binary_mixtures")


def s9():
    ridge = pd.read_csv(DATA/"model_all_ridge_coefficients.csv").sort_values("standardized_coefficient")
    gam = pd.read_csv(DATA/"model_all_importance.csv").sort_values("importance")
    fig, axes = plt.subplots(1, 2, figsize=(7.7, 4))
    axes[0].barh(ridge.feature, ridge.standardized_coefficient, color=np.where(ridge.standardized_coefficient >= 0, BLUE, "#C44E52"))
    axes[0].axvline(0, color="black", lw=.8); axes[0].set(xlabel="Standardized Ridge coefficient", title="A. Ridge coefficients")
    axes[1].barh(gam.feature, gam.importance, xerr=gam.importance_sd, color="#B07AA1", error_kw={"elinewidth": .8, "capsize": 2})
    axes[1].set(xlabel=r"GAM permutation $\Delta R^2$", title="B. GAM permutation importance")
    for ax in axes: clean(ax)
    fig.tight_layout(); save(fig, "Figure_S9_descriptor_importance")


def s10():
    order = ["Pyrene", "Ethylene_glycol", "Perfluorooctanoic_acid", "2_4-Dichlorophenoxyacetic_acid", "Nitrate"]
    display = {"Pyrene": "Pyrene (neutral)", "Ethylene_glycol": "Ethylene glycol (neutral)",
               "Perfluorooctanoic_acid": "Perfluorooctanoate (anionic)", "2_4-Dichlorophenoxyacetic_acid": "2,4-D (anionic)",
               "Nitrate": "Nitrate (anionic)"}
    stacked("surface_charge_control_summary.csv", "charge_magnitude_e", [10, 25, 35, 50], order, display,
            dict(zip(order, [BLUE, GREEN, RED, ORANGE, PURPLE])), "PET-NPP net charge (e)", "Figure_S10_surface_charge_control", True)


def s11():
    corr = pd.read_csv(DATA/"mw_descriptor_correlations.csv").sort_values("spearman_rho")
    profiles = pd.read_csv(DATA/"mw_breakpoint_profiles.csv"); summary = pd.read_csv(DATA/"mw_transition_summary.csv")
    compounds = pd.read_csv(DATA/"compound_descriptors.csv"); neutral = compounds.loc[compounds.NetCharge.eq(0), ["MW", "occupancy_mean"]]
    fig, axes = plt.subplots(1, 3, figsize=(12.2, 3.7))
    bar_colors = ["#C44E52" if abs(value) >= .5 else BLUE for value in corr.spearman_rho]
    axes[0].barh(corr.descriptor, corr.spearman_rho, color=bar_colors, alpha=.9)
    axes[0].axvline(0, color="black", lw=.8); axes[0].set(xlabel=r"Spearman $\rho$ with MW", xlim=(-.1, .82))
    x = np.linspace(neutral.MW.min(), neutral.MW.max(), 400)
    colors = {"unadjusted": "#222222", "class-adjusted": "#C44E52"}
    for row in summary.itertuples(index=False):
        below, above, knot = row.slope_below_per_100_Da/100, row.slope_above_per_100_Da/100, row.best_transition_Da
        basis = below*neutral.MW+(above-below)*np.maximum(0, neutral.MW-knot); intercept = neutral.occupancy_mean.mean()-basis.mean()
        y = intercept+below*x+(above-below)*np.maximum(0, x-knot)
        axes[1].plot(x, y, color=colors[row.model], lw=2.2,
                     ls="--" if row.model == "unadjusted" else "-",
                     label=f"{row.model}: {knot:.0f} Da")
    axes[1].scatter(neutral.MW, neutral.occupancy_mean, s=18, alpha=.35, color="#6E87A8", edgecolors="none")
    axes[1].set(xlabel="Molecular weight (Da)", ylabel="Occupancy", ylim=(-.03, 1.04)); axes[1].legend(frameon=False, fontsize=9)
    for model, group in profiles.groupby("model", sort=False):
        axes[2].plot(group.knot_Da, group.rss-group.rss.min(), lw=2, color=colors[model], label=model)
        knot = summary.loc[summary.model.eq(model), "best_transition_Da"].iloc[0]
        axes[2].axvline(knot, color=colors[model], ls=":" if model == "class-adjusted" else "--", alpha=.8)
    axes[2].set(xlabel="Candidate transition (Da)", ylabel="RSS above model minimum"); axes[2].legend(frameon=False, fontsize=9)
    for ax, title in zip(axes, ["A. Descriptor collinearity", "B. Piecewise MW response", "C. Transition-profile sensitivity"]):
        ax.set_title(title, loc="left", fontweight="bold"); clean(ax)
    fig.tight_layout(); save(fig, "Figure_S11_MW_class_robustness")


def s12():
    grid = pd.read_csv(DATA/"resp_esp_grid.csv"); stats = pd.read_csv(DATA/"resp_esp_summary.csv").set_index("molecule")
    order = ["Pyrene", "Phenol", "Ethylene_glycol"]; fig, axes = plt.subplots(2, 3, figsize=(10.2, 6.2)); factor = 627.509474
    for col, molecule in enumerate(order):
        d = grid.loc[grid.molecule.eq(molecule)]; x = d.resp_esp_hartree_per_e*factor; y = d.am1bcc_esp_hartree_per_e*factor
        low, high = min(x.min(), y.min()), max(x.max(), y.max()); axes[0, col].scatter(x, y, s=5, alpha=.22, color=BLUE, edgecolors="none")
        axes[0, col].plot([low, high], [low, high], "--", color="black", lw=.8); row = stats.loc[molecule]
        axes[0, col].text(.04, .95, f"r = {row.pearson_r:.3f}\nRMSE = {row.esp_rmse_kcal_mol_per_e:.2f}", transform=axes[0, col].transAxes, va="top", fontsize=8)
        axes[0, col].set(title=molecule.replace("_", " "), xlabel="RESP ESP (kcal mol$^{-1}$ e$^{-1}$)",
                         ylabel="AM1-BCC ESP (kcal mol$^{-1}$ e$^{-1}$)")
        delta = d.delta_am1bcc_minus_resp_hartree_per_e*factor; axes[1, col].hist(delta, bins=35, color=GREEN, edgecolor="white")
        axes[1, col].axvline(0, color="black", lw=.8); axes[1, col].set(xlabel="AM1-BCC − RESP ESP\n(kcal mol$^{-1}$ e$^{-1}$)", ylabel="Grid points")
        clean(axes[0, col]); clean(axes[1, col])
    for label, ax in zip("ABC", axes[0]): ax.text(-.13, 1.10, label, transform=ax.transAxes, fontweight="bold", fontsize=12)
    for label, ax in zip("DEF", axes[1]): ax.text(-.13, 1.08, label, transform=ax.transAxes, fontweight="bold", fontsize=12)
    fig.tight_layout(); save(fig, "Figure_S12_RESP_AM1BCC")


def s13():
    d = pd.read_csv(DATA/"multicopy_independence.csv"); fig, axes = plt.subplots(1, 3, figsize=(8.6, 3.35))
    limit = max(d.sigma_predicted.max(), d.sigma_between.max())*1.06
    axes[0].scatter(d.sigma_predicted, d.sigma_between, s=18, alpha=.62, color=BLUE, edgecolors="none"); axes[0].plot([0, limit], [0, limit], "--", color="black", lw=1.1)
    ratio = (d.sigma_predicted/d.sigma_between).replace([np.inf, -np.inf], np.nan).median()
    axes[0].text(.03, .97, f"median ratio = {ratio:.2f}", transform=axes[0].transAxes, va="top", bbox={"boxstyle":"round,pad=.2", "facecolor":"white", "edgecolor":"#BBBBBB"})
    axes[0].set(xlim=(-.005, limit), ylim=(-.005, limit), xlabel=r"Predicted $\sigma_B$ ($\sigma_W/\sqrt{10}$)", ylabel=r"Observed $\sigma_B$")
    axes[1].hist(d.ICC.dropna(), bins=np.linspace(-.12, .52, 25), color=GREEN, edgecolor="white", linewidth=.5); axes[1].axvline(0, color="#777777", ls="--", lw=1); axes[1].axvline(d.ICC.median(), color="black", lw=1.5)
    axes[1].text(.03, .97, f"median = {d.ICC.median():.3f}", transform=axes[1].transAxes, va="top", bbox={"boxstyle":"round,pad=.2", "facecolor":"white", "edgecolor":"#BBBBBB"}); axes[1].set(xlabel="Intraclass correlation (ICC)", ylabel="System count")
    axes[2].hist(d.perm_p.dropna(), bins=np.linspace(0, 1, 21), color="#C44E52", edgecolor="white", linewidth=.5); axes[2].axhline(len(d)/20, color="#777777", ls="--", lw=1)
    hits = int(d.perm_p_bh_reject.sum()); axes[2].text(.97, .97, f"BH hits = {hits}/203", transform=axes[2].transAxes, ha="right", va="top", bbox={"boxstyle":"round,pad=.2", "facecolor":"white", "edgecolor":"#BBBBBB"}); axes[2].set(xlim=(0, 1), xlabel="Permutation p-value", ylabel="System count")
    for ax, title in zip(axes, ["A  Variance scaling", "B  Replica clustering", "C  Replica-label test"]):
        ax.set_title(title, loc="left", fontweight="bold"); ax.spines[["top", "right"]].set_visible(False); ax.grid(axis="y", alpha=.18)
    fig.tight_layout(); save(fig, "Figure_S13_multicopy_independence")


def main():
    configure_style()
    for function in (s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13): function()
    print("Generated submitted numerical Supporting-Information Figures S2-S13")


if __name__ == "__main__":
    main()
