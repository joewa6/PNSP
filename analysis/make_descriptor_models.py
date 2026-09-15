#!/usr/bin/env python3
"""Compute RDKit descriptors and reproduce the Ridge/GAM screen models."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pygam import LinearGAM, s
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors
from sklearn.impute import SimpleImputer
from sklearn.linear_model import RidgeCV
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import KFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


FEATURES = [
    "MW", "LogP", "TPSA", "HBD", "HBA", "RotatableBonds",
    "NumAromaticRings", "HeavyAtomCount", "NumHalogens", "NetCharge",
]
NEUTRAL_FEATURES = [feature for feature in FEATURES if feature != "NetCharge"]
SEED = 42


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip()).strip("_.")


def descriptor_row(row: pd.Series) -> dict[str, object]:
    selected = row.get("selected_smiles")
    smiles = str(selected if pd.notna(selected) and str(selected).strip() else row["smiles"])
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"RDKit could not parse {row['name']}: {smiles}")
    calculated_charge = int(sum(atom.GetFormalCharge() for atom in molecule.GetAtoms()))
    charge = int(row["formal_charge"]) if "formal_charge" in row and pd.notna(row["formal_charge"]) else calculated_charge
    if charge != calculated_charge:
        raise ValueError(
            f"{row['name']}: formal_charge={charge}, but SMILES encodes {calculated_charge}"
        )
    return {
        "system": safe_name(str(row["name"])),
        "smiles": Chem.MolToSmiles(molecule, isomericSmiles=True),
        "MW": Descriptors.MolWt(molecule),
        "LogP": Crippen.MolLogP(molecule),
        "TPSA": rdMolDescriptors.CalcTPSA(molecule),
        "HBD": Lipinski.NumHDonors(molecule),
        "HBA": Lipinski.NumHAcceptors(molecule),
        "RotatableBonds": Lipinski.NumRotatableBonds(molecule),
        "NumAromaticRings": Lipinski.NumAromaticRings(molecule),
        "HeavyAtomCount": molecule.GetNumHeavyAtoms(),
        "NumHalogens": sum(
            atom.GetAtomicNum() in {9, 17, 35, 53} for atom in molecule.GetAtoms()
        ),
        "NetCharge": charge,
    }


def gam_terms(features: int, splines: int):
    terms = s(0, n_splines=splines)
    for index in range(1, features):
        terms += s(index, n_splines=splines)
    return terms


def fit_gam(X: np.ndarray, y: np.ndarray) -> LinearGAM:
    best: LinearGAM | None = None
    best_gcv = np.inf
    for splines in (5, 6):
        model = LinearGAM(gam_terms(X.shape[1], splines))
        model.gridsearch(X, y, lam=np.logspace(-2, 2, 5), progress=False)
        gcv = float(model.statistics_.get("GCV", np.inf))
        if gcv < best_gcv:
            best, best_gcv = model, gcv
    if best is None:
        raise RuntimeError("All GAM fits failed")
    return best


def gam_oof(X: pd.DataFrame, y: np.ndarray, folds: int) -> np.ndarray:
    cv = KFold(n_splits=folds, shuffle=True, random_state=SEED)
    predictions = np.zeros(len(y), dtype=float)
    X = X.reset_index(drop=True)
    for train, test in cv.split(X):
        imputer = SimpleImputer(strategy="median")
        scaler = StandardScaler()
        train_X = scaler.fit_transform(imputer.fit_transform(X.iloc[train]))
        test_X = scaler.transform(imputer.transform(X.iloc[test]))
        predictions[test] = fit_gam(train_X, y[train]).predict(test_X)
    return predictions


def ridge_oof(X: pd.DataFrame, y: np.ndarray, folds: int) -> np.ndarray:
    cv = KFold(n_splits=folds, shuffle=True, random_state=SEED)
    model = Pipeline(
        [
            ("impute", SimpleImputer(strategy="median")),
            ("scale", StandardScaler()),
            ("ridge", RidgeCV(alphas=np.logspace(-3, 3, 40))),
        ]
    )
    return cross_val_predict(model, X, y, cv=cv, n_jobs=1)


def fitted_gam_importance(
    X: pd.DataFrame, y: np.ndarray, features: list[str], repeats: int = 40
) -> pd.DataFrame:
    imputer = SimpleImputer(strategy="median")
    scaler = StandardScaler()
    scaled = scaler.fit_transform(imputer.fit_transform(X[features]))
    model = fit_gam(scaled, y)
    baseline = r2_score(y, model.predict(scaled))
    rng = np.random.default_rng(SEED)
    rows = []
    clean = X[features].reset_index(drop=True)
    for feature in features:
        deltas = []
        for _ in range(repeats):
            permuted = clean.copy()
            permuted[feature] = rng.permutation(permuted[feature].to_numpy())
            transformed = scaler.transform(imputer.transform(permuted))
            deltas.append(baseline - r2_score(y, model.predict(transformed)))
        rows.append(
            {
                "feature": feature,
                "importance_delta_r2": float(np.mean(deltas)),
                "importance_sd": float(np.std(deltas, ddof=1)),
                "permutation_repeats": repeats,
            }
        )
    return pd.DataFrame(rows).sort_values("importance_delta_r2", ascending=False)


def metric_row(scope: str, model: str, y: np.ndarray, prediction: np.ndarray) -> dict[str, object]:
    return {
        "scope": scope,
        "model": model,
        "n": len(y),
        "heldout_r2": float(r2_score(y, prediction)),
        "heldout_rmse": float(np.sqrt(mean_squared_error(y, prediction))),
        "heldout_mae": float(mean_absolute_error(y, prediction)),
    }


def analyze_scope(
    frame: pd.DataFrame,
    features: list[str],
    scope: str,
    folds: int,
    output: Path,
) -> tuple[list[dict[str, object]], pd.DataFrame]:
    y = frame["occupancy_mean"].to_numpy(float)
    predictions = {
        "ridge_logp": ridge_oof(frame[["LogP"]], y, folds),
        "gam_logp": gam_oof(frame[["LogP"]], y, folds),
        "ridge_all": ridge_oof(frame[features], y, folds),
        "gam_all": gam_oof(frame[features], y, folds),
    }
    prediction_table = pd.DataFrame(
        {"system": frame["system"], "actual_occupancy": y, **predictions}
    )
    prediction_table.to_csv(output / f"{scope}_oof_predictions.csv", index=False)
    importance = fitted_gam_importance(frame[features], y, features)
    importance.to_csv(output / f"{scope}_gam_permutation_importance.csv", index=False)
    metrics = [metric_row(scope, name, y, values) for name, values in predictions.items()]
    return metrics, importance


def plot_full_parity(frame: pd.DataFrame, output: Path) -> None:
    predictions = pd.read_csv(output / "full_oof_predictions.csv")
    specifications = [
        ("ridge_logp", "A. Ridge (LogP only)"),
        ("gam_logp", "B. GAM (LogP only)"),
        ("ridge_all", "C. Ridge (all descriptors)"),
        ("gam_all", "D. GAM (all descriptors)"),
    ]
    figure, axes = plt.subplots(2, 2, figsize=(7.2, 6.8))
    actual = predictions["actual_occupancy"].to_numpy()
    for axis, (column, title) in zip(axes.flat, specifications):
        predicted = predictions[column].to_numpy()
        low = min(actual.min(), predicted.min())
        high = max(actual.max(), predicted.max())
        axis.scatter(actual, predicted, s=22, alpha=0.75, color="#4C78A8")
        axis.plot([low, high], [low, high], "--", color="black", linewidth=1)
        axis.set(title=title, xlabel="Actual occupancy", ylabel="Predicted occupancy")
        axis.text(0.03, 0.95, f"$R^2$ = {r2_score(actual, predicted):.3f}", transform=axis.transAxes, va="top")
        axis.grid(alpha=0.2)
    figure.tight_layout()
    figure.savefig(output / "full_model_parity.png", dpi=300)
    figure.savefig(output / "full_model_parity.pdf")
    plt.close(figure)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--compounds", required=True, type=Path)
    parser.add_argument("--output", type=Path, default=Path("descriptor_models"))
    parser.add_argument("--folds", type=int, default=5)
    args = parser.parse_args()

    summary = pd.read_csv(args.summary)
    compounds = pd.read_csv(args.compounds)
    for column in ("name", "smiles"):
        if column not in compounds:
            raise ValueError(f"Compound table is missing {column!r}")
    descriptors = pd.DataFrame([descriptor_row(row) for _, row in compounds.iterrows()])
    if descriptors["system"].duplicated().any():
        raise ValueError("Compound identifiers are not unique")
    frame = summary[["system", "occupancy_mean"]].merge(
        descriptors, on="system", how="inner", validate="one_to_one"
    )
    if len(frame) != len(summary) or len(frame) != len(descriptors):
        missing_summary = sorted(set(descriptors.system) - set(summary.system))
        missing_descriptors = sorted(set(summary.system) - set(descriptors.system))
        raise ValueError(
            f"Summary/descriptor mismatch; absent from summary={missing_summary}; "
            f"absent from compounds={missing_descriptors}"
        )
    if len(frame) < max(12, args.folds * 2):
        raise ValueError("At least 12 matched compounds are required for descriptor modelling")

    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output / "descriptor_model_input.csv", index=False)
    metrics, _ = analyze_scope(frame, FEATURES, "full", args.folds, output)
    neutral = frame.loc[frame["NetCharge"].eq(0)].reset_index(drop=True)
    if len(neutral) >= max(12, args.folds * 2):
        neutral_metrics, _ = analyze_scope(
            neutral, NEUTRAL_FEATURES, "neutral", args.folds, output
        )
        metrics.extend(neutral_metrics)
    pd.DataFrame(metrics).to_csv(output / "model_performance.csv", index=False)
    plot_full_parity(frame, output)
    provenance = {
        "systems": len(frame),
        "neutral_systems": len(neutral),
        "descriptors": FEATURES,
        "cv": f"KFold(n_splits={args.folds}, shuffle=True, random_state=42)",
        "ridge": "median imputation and standardization inside each outer fold; RidgeCV alpha grid logspace(-3,3,40)",
        "gam": "median imputation and standardization inside each outer fold; additive LinearGAM; 5 or 6 splines and lambda logspace(-2,2,5) selected by training-fold GCV",
        "permutation_importance": "separate full-data fitted GAM; baseline minus permuted in-sample R2; 40 repeats per descriptor; seed 42",
    }
    (output / "model_provenance.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(pd.DataFrame(metrics).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
