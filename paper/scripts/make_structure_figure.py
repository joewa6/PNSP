#!/usr/bin/env python3
"""Render Figure S1 from the released PET structural measurements."""
from common import DATA, OUTPUT
import matplotlib.pyplot as plt
import pandas as pd

CLASS_LABELS = ["Aromatic", "Aliphatic", "Carboxylate", "Ester"]
COLORS = ["#4C78A8", "#F58518", "#9E2F3B", "#72B7B2"]

def draw_figure(atoms, groups, composition, summary, output):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"})
    fig = plt.figure(figsize=(9.4, 7.5))
    grid = fig.add_gridspec(2, 2, left=.08, right=.98, bottom=.08, top=.94,
                           hspace=.35, wspace=.29)
    ax_a, ax_b, ax_c, ax_d = [fig.add_subplot(grid[i, j]) for i, j in [(0,0),(0,1),(1,0),(1,1)]]
    # All reference atoms are drawn; coordinates and masks also exported to CSV.
    def scatter(ax, frame, color, size=3.0, alpha=.45, **kwargs):
        return ax.scatter(frame.pc1, frame.pc2, c=color, s=size, alpha=alpha,
                          linewidths=0, rasterized=True, **kwargs)
    xlim = (float(min(-44, atoms.pc1.min()-3)), float(max(44, atoms.pc1.max()+3)))
    ylim = (float(min(-27, atoms.pc2.min()-3)), 40)
    for ax in (ax_a, ax_b, ax_d):
        ax.set(xlim=xlim, ylim=ylim, xlabel="PC1 (Å)", ylabel="PC2 (Å)")
        ax.set_aspect("equal", adjustable="box")
        ax.set_xticks([-40, -20, 0, 20, 40])
        ax.set_yticks([-20, 0, 20, 40])
        ax.tick_params(labelsize=8)

    scatter(ax_a, atoms, "#C7CED6", size=3.0, alpha=.48)
    a, b = summary["major_axis_A"]/2, summary["intermediate_axis_A"]/2
    ax_a.plot([-a, a], [0, 0], color="#0072B2", lw=2.1)
    ax_a.plot([0, 0], [-b, b], color="#56B4E9", lw=1.8)
    ax_a.set_title("Global morphology", fontsize=10, weight="bold", pad=8)
    ax_a.text(.02, .97, f"$R_g$ = {summary['radius_of_gyration_A']:.1f} Å\n"
              f"Effective axes = {summary['major_axis_A']:.1f}, {summary['intermediate_axis_A']:.1f}, "
              f"{summary['minor_axis_A']:.1f} Å\nMajor/minor ratio = {summary['axis_ratio_major_minor']:.1f}",
              transform=ax_a.transAxes, fontsize=8, va="top")

    scatter(ax_b, atoms.loc[~atoms.is_exposed], "#D0D0D0", size=2.5, alpha=.18)
    scatter(ax_b, atoms.loc[atoms.is_exposed], "#E74C3C", size=4, alpha=.72)
    ax_b.set_title("Geometric exposure proxy", fontsize=10, weight="bold", pad=8)
    # Explicit legend handles keep the unselected marker visible without changing the data layer.
    handles = [plt.Line2D([], [], linestyle="", marker="o", color=c, markersize=4, label=label)
               for c, label in [("#B0B0B0", "Not selected"), ("#E74C3C", "Proxy-selected")]]
    ax_b.legend(handles=handles, loc="upper left", fontsize=7.5, frameon=False)
    ax_b.text(.98, .97, f"{summary['total_exposed_atoms']:,}/{summary['n_atoms']:,} atoms\n"
              f"({summary['frac_exposed_atoms']:.1%}) selected", transform=ax_b.transAxes,
              ha="right", va="top", fontsize=7.7)

    values = composition.fraction.to_numpy()
    wedges, _ = ax_c.pie(values, colors=COLORS, startangle=90, counterclock=False,
                         radius=.78, center=(0, .12),
                         wedgeprops=dict(width=.48, edgecolor="white", linewidth=1.2))
    ax_c.set_title("Composition of proxy-selected atoms", fontsize=10, weight="bold", pad=8)
    ax_c.legend(wedges, [f"{label} ({value:.1%})" for label, value in zip(CLASS_LABELS, values)],
                loc="lower center", bbox_to_anchor=(.5, -.02), ncol=2, fontsize=8,
                frameon=False, columnspacing=1.2, handlelength=1.5)
    ax_c.set(xlim=(-1.08, 1.08), ylim=(-1.05, 1.05))
    ax_c.set_aspect("equal")

    scatter(ax_d, atoms, "#BFC5CC", size=2.2, alpha=.20)
    selected = groups.loc[groups.is_exposed]
    ax_d.scatter(selected.pc1, selected.pc2, s=42, c="#C44E52", edgecolors="#222222",
                 linewidths=.6, alpha=.94, zorder=5)
    ax_d.set_title("Proxy-selected terminal carboxylates", fontsize=10, weight="bold", pad=8)
    ax_d.text(.02, .97, f"{len(selected)}/{len(groups)} groups selected\n"
              f"Mean 3D NN spacing\n= {summary['mean_terminal_carboxylate_nn_distance_A']:.2f} Å",
              transform=ax_d.transAxes, fontsize=7.7, va="top")
    # The inset occupies the blank header band, not the carboxylate projection.
    inset = ax_d.inset_axes([.62, .81, .35, .17])
    nn = selected.nearest_selected_neighbor_A.to_numpy()
    inset.hist(nn, bins=12, color="#C44E52", edgecolor="white", linewidth=.4, alpha=.85)
    inset.axvline(nn.mean(), color="#333333", linestyle="--", lw=1)
    inset.set_xlabel("3D NN distance (Å)", fontsize=6, labelpad=1)
    inset.set_ylabel("Groups", fontsize=6, labelpad=1)
    inset.tick_params(labelsize=6, length=2, pad=1)
    inset.set_xticks([5, 10, 15])

    # Figure-coordinate panel letters avoid collisions with ticks or titles.
    fig.canvas.draw()
    for ax, label in zip((ax_a, ax_b, ax_c, ax_d), "ABCD"):
        box = ax.get_position()
        fig.text(box.x0-.047, box.y1+.016, label, fontsize=14, weight="bold", va="bottom")
    stem = output / "Figure_S1_PET_structural_characterization"
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(stem.with_suffix("."+suffix), dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(fig)



def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    summary = pd.read_csv(DATA / "pet_model_summary.csv").iloc[0].to_dict()
    for key in ("n_atoms", "n_residues", "total_exposed_atoms", "total_terminal_carboxylates",
                "exposed_terminal_carboxylates"):
        summary[key] = int(summary[key])
    draw_figure(pd.read_csv(DATA / "pet_atom_geometry.csv"),
                pd.read_csv(DATA / "pet_terminal_carboxylates.csv"),
                pd.read_csv(DATA / "pet_selected_composition.csv"),
                summary, OUTPUT)

if __name__ == "__main__":
    main()
