#!/usr/bin/env python3
"""Regenerate manuscript Figures 3--8 from the released processed data."""

from __future__ import annotations

from io import BytesIO
import sys
import textwrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from sklearn.impute import SimpleImputer
from sklearn.metrics import r2_score
from sklearn.preprocessing import StandardScaler

from common import DATA, REPOSITORY, bootstrap_mean, configure_style, panel_label, save


CLASS_ORDER = [
    "Phthalates", "PAHs", "Organochlorine pesticides",
    "Triazine herbicides", "Organophosphates", "Trihalomethanes", "PFAS",
]
FEATURES_NEUTRAL = [
    "MW", "LogP", "TPSA", "HBD", "HBA", "RotatableBonds",
    "NumAromaticRings", "HeavyAtomCount", "NumHalogens",
]


def figure_3(descriptors: pd.DataFrame) -> None:
    metadata = [
        ("MW", "Mol Weight (Da)"), ("LogP", "LogP"),
        ("TPSA", r"TPSA ($\AA^2$)"), ("HBD", "H-bond donors"),
        ("HBA", "H-Bond Acceptors"), ("RotatableBonds", "Rotatable Bonds"),
        ("NumAromaticRings", "Aromatic Ring Count"), ("NumHalogens", "Halogen Count"),
    ]
    rng = np.random.default_rng(42)
    figure, axes = plt.subplots(2, 4, figsize=(8.0, 4.4), constrained_layout=True)
    for axis, (column, label) in zip(axes.flat, metadata):
        values = descriptors[column].to_numpy(float)
        axis.boxplot(values, widths=0.40, patch_artist=True, showfliers=False,
                     boxprops={"facecolor": "#4c72b0", "alpha": 0.35},
                     medianprops={"color": "#2a4f7a"})
        axis.scatter(rng.normal(1, 0.035, len(values)), values, s=6,
                     alpha=0.16, color="#4c72b0", linewidth=0)
        axis.set(xticks=[], xlabel=label, ylabel="Value")
    save(figure, "Figure_3_descriptor_distributions")


def figure_4() -> None:
    frame = pd.read_csv(DATA / "figure4_bootstrap.csv").sort_values("rank")
    figure, axis = plt.subplots(figsize=(7.0, 3.3))
    axis.errorbar(frame["rank"], frame["mean_occ"],
                  yerr=[frame.mean_occ-frame.ci90_lo, frame.ci90_hi-frame.mean_occ],
                  fmt="o", ms=2.3, color="#2a5f8a", ecolor="#9fc0e5",
                  elinewidth=0.55, capsize=0)
    axis.plot(frame["rank"], frame["lower90_bound"], color="0.55", lw=1.1)
    axis.set(xlabel="Compound rank (by mean occupancy)", ylabel="Occupancy", ylim=(-0.02, 1.03),
             xlim=(0, len(frame)+1))
    save(figure, "Figure_4_rank_resolved_occupancy")


def molecule_image(smiles: str, size=(700, 430)):
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    AllChem.Compute2DCoords(molecule)
    drawer = rdMolDraw2D.MolDraw2DCairo(*size)
    options = drawer.drawOptions()
    options.useBWAtomPalette()
    options.bondLineWidth = 1.35
    drawer.DrawMolecule(molecule)
    drawer.FinishDrawing()
    return np.asarray(Image.open(BytesIO(drawer.GetDrawingText())).convert("RGB"))


def figure_5() -> None:
    examples = pd.read_csv(DATA / "figure5_examples.csv")
    figure = plt.figure(figsize=(6.2, 6.1), facecolor="white")
    grid = figure.add_gridspec(3, 4, height_ratios=(1, 1, 0.52), hspace=0.58, wspace=0.22)
    axes = [figure.add_subplot(grid[0,:2]), figure.add_subplot(grid[0,2:]),
            figure.add_subplot(grid[1,:2]), figure.add_subplot(grid[1,2:]),
            figure.add_subplot(grid[2,1:3])]
    grey = "#555555"
    for index, (axis, row) in enumerate(zip(axes, examples.itertuples(index=False))):
        if row.system == "Ammonia":
            axis.text(0.5, 0.60, r"$\mathrm{NH_4^+}$", transform=axis.transAxes,
                      ha="center", va="center", fontsize=13, color=grey)
            name_y, occupancy_y = 0.27, 0.02
        else:
            axis.imshow(molecule_image(row.smiles))
            name_y, occupancy_y = -0.01, -0.22
        axis.axis("off")
        axis.set_title(chr(65 + index), loc="left", fontsize=10, fontweight="bold")
        axis.text(0.5, name_y, "\n".join(textwrap.wrap(row.display_name, 24)),
                  transform=axis.transAxes, ha="center", va="top", fontsize=10.2,
                  fontstyle="italic", color=grey)
        axis.text(0.5, occupancy_y, f"Occupancy = {row.occupancy_mean:.3f}",
                  transform=axis.transAxes, ha="center", va="top", fontsize=8.8, color=grey)
    figure.subplots_adjust(left=.035, right=.965, top=.985, bottom=.055)
    save(figure, "Figure_5_representative_association_modes")


def figure_6(descriptors: pd.DataFrame, copies: pd.DataFrame) -> None:
    system_class = descriptors.dropna(subset=["chemical_class"])[["system","chemical_class","occupancy_mean"]]
    class_map = dict(zip(system_class.system, system_class.chemical_class))
    selected = copies.loc[copies.system.isin(class_map)].copy()
    selected["chemical_class"] = selected.system.map(class_map)
    rng = np.random.default_rng(42)
    figure, axis = plt.subplots(figsize=(7.0, 3.3))
    for index, group in enumerate(CLASS_ORDER):
        raw = selected.loc[selected.chemical_class.eq(group), "occupancy"].to_numpy()
        axis.scatter(rng.normal(index, .06, len(raw)), raw, s=9, alpha=.27,
                     color="#4c72b0", linewidth=0)
        values = system_class.loc[system_class.chemical_class.eq(group), "occupancy_mean"]
        mean, low, high = bootstrap_mean(values, 5000, 42 + index)
        axis.errorbar(index, mean, yerr=[[mean-low],[high-mean]], fmt="o", color="black",
                      capsize=4, ms=5)
        axis.text(index, 1.025, f"n={len(values)}", ha="center", va="bottom", fontsize=9)
    axis.set(xticks=range(len(CLASS_ORDER)), xticklabels=CLASS_ORDER,
             xlabel="Chemical class", ylabel="Occupancy", ylim=(-.02,1.07))
    axis.tick_params(axis="x", rotation=22)
    for label in axis.get_xticklabels(): label.set_ha("right")
    save(figure, "Figure_6_chemical_classes")


def parity(axis, actual, predicted, title, label, color="#4c78a8"):
    score = r2_score(actual, predicted)
    axis.scatter(actual, predicted, s=22, alpha=.75, color=color, edgecolor="white", lw=.3)
    low=min(np.min(actual),np.min(predicted)); high=max(np.max(actual),np.max(predicted))
    axis.plot([low,high],[low,high],"--",color="black",lw=1)
    axis.set(xlabel="Actual occupancy", ylabel="Predicted occupancy", title=f"{label}. {title}")
    axis.text(.04,.94,f"$R^2$ = {score:.3f}\nn = {len(actual)}",transform=axis.transAxes,va="top")


def figure_7() -> None:
    data = pd.read_csv(DATA / "model_all_predictions.csv")
    figure, axes = plt.subplots(2,2,figsize=(7.2,6.8))
    specs=[("ridge_logp_oof","Ridge (LogP only)"),("gam_logp_oof","GAM (LogP only)"),
           ("ridge_all_oof","Ridge (all descriptors)"),("gam_all_oof","GAM (all descriptors)")]
    for label, axis, (column,title) in zip("ABCD",axes.flat,specs):
        parity(axis,data.actual,data[column],title,label)
    figure.tight_layout()
    save(figure,"Figure_7_model_performance")


def figure_8(descriptors: pd.DataFrame) -> None:
    predictions=pd.read_csv(DATA/"model_neutral_predictions.csv")
    importance=pd.read_csv(DATA/"model_neutral_importance.csv").sort_values("importance")
    neutral=descriptors.loc[descriptors.NetCharge.eq(0)].merge(predictions[["system"]],on="system")
    y=neutral.occupancy_mean.to_numpy(float)
    sys.path.insert(0,str(REPOSITORY))
    from analysis.make_descriptor_models import fit_gam
    imputer=SimpleImputer(strategy="median"); scaler=StandardScaler()
    scaled=scaler.fit_transform(imputer.fit_transform(neutral[FEATURES_NEUTRAL]))
    model=fit_gam(scaled,y)
    mw_index=FEATURES_NEUTRAL.index("MW")
    grid=model.generate_X_grid(term=mw_index,n=200)
    partial,ci=model.partial_dependence(term=mw_index,X=grid,width=.90)
    mw=scaler.inverse_transform(grid)[:,mw_index]; order=np.argsort(mw)

    figure=plt.figure(figsize=(9.4,14.5)); layout=figure.add_gridspec(
        3,2,height_ratios=[1.55,1.35,1.55],hspace=.50,wspace=.48
    )
    a=figure.add_subplot(layout[0,0]); b=figure.add_subplot(layout[0,1]); c=figure.add_subplot(layout[1,:]); d=figure.add_subplot(layout[2,0]); e=figure.add_subplot(layout[2,1])
    parity(a,predictions.actual,predictions.ridge_all_oof,"Ridge (all descriptors)","A","#55a868")
    parity(b,predictions.actual,predictions.gam_all_oof,"GAM (all descriptors)","B","#55a868")
    c.barh(importance.feature,importance.importance,xerr=importance.importance_sd,color="#c44e52",alpha=.9)
    c.set(xlabel=r"Permutation $\Delta R^2$",title="C. Feature importance")
    per_compound_partial,_=model.partial_dependence(term=mw_index,X=scaled,width=.90)
    residuals=y-model.predict(scaled)+per_compound_partial
    residual_frame=pd.DataFrame({"MW":neutral.MW.to_numpy(float),"residual":residuals})
    residual_frame["bin"]=pd.qcut(
        residual_frame.MW,min(12,residual_frame.MW.nunique()),duplicates="drop"
    )
    residual_summary=(
        residual_frame.groupby("bin",observed=True)
        .agg(MW=("MW","median"),residual=("residual","mean"))
    )
    d.fill_between(mw[order],ci[order,0],ci[order,1],color="#1f4f8b",alpha=.2)
    d.plot(mw[order],partial[order],color="#1f4f8b",lw=2)
    d.scatter(residual_summary.MW,residual_summary.residual,s=36,color="black",
              edgecolor="white",linewidth=.4,zorder=5)
    d.axhline(0,color="grey",ls=":",lw=.8); d.set(xlabel="Molecular weight (Da)",ylabel="Partial effect on occupancy",title="D. GAM smooth: MW")
    e.scatter(neutral.MW,y,s=22,alpha=.7,color="#55a868")
    bins=pd.qcut(neutral.MW,min(12,neutral.MW.nunique()),duplicates="drop")
    trend=neutral.assign(_bin=bins).groupby("_bin",observed=True).agg(MW=("MW","median"),occupancy=("occupancy_mean","mean"))
    e.plot(trend.MW,trend.occupancy,color="black",lw=1.7)
    e.set(xlabel="Molecular weight (Da)",ylabel="Observed occupancy",title="E. MW vs occupancy",ylim=(0,1.05))
    save(figure,"Figure_8_neutral_descriptor_model")


def main() -> None:
    configure_style()
    descriptors=pd.read_csv(DATA/"compound_descriptors.csv")
    copies=pd.read_csv(DATA/"copy_metrics.csv")
    if len(descriptors)!=203 or len(copies)!=6090:
        raise RuntimeError("Expected 203 compounds and 6,090 copy trajectories")
    figure_3(descriptors); figure_4(); figure_5(); figure_6(descriptors,copies)
    figure_7(); figure_8(descriptors)
    print("Generated manuscript Figures 3--8")


if __name__ == "__main__":
    main()
