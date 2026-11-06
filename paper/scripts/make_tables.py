#!/usr/bin/env python3
"""Regenerate submitted Supporting-Information Tables S1 and S2."""
from __future__ import annotations

import pandas as pd
from common import DATA, OUTPUT

FEATURE_NAMES = {
    "NetCharge": "Net charge",
    "LogP": "logP",
    "TPSA": "TPSA",
    "MW": "Molecular weight",
    "RotatableBonds": "Rotatable bonds",
    "NumAromaticRings": "Aromatic rings",
    "HeavyAtomCount": "Heavy-atom count",
    "HBA": "H-bond acceptors",
    "NumHalogens": "Halogen count",
    "HBD": "H-bond donors",
}


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)

    importance = pd.read_csv(DATA/"model_all_importance.csv")
    table_s1 = pd.DataFrame({
        "Descriptor": importance.feature.map(FEATURE_NAMES),
        "Mean decrease in R2 (delta R2)": importance.importance.round(3),
        "SD": importance.importance_sd.round(3),
        "Permutation repeats": 40,
    })
    if table_s1.Descriptor.isna().any() or len(table_s1) != 10:
        raise RuntimeError("Table S1 requires the ten canonical descriptors")
    table_s1.to_csv(OUTPUT/"Table_S1_GAM_permutation_importance.csv", index=False,
                    float_format="%.3f")

    screen = pd.read_csv(DATA/"screen_summary.csv").sort_values("rank")
    if len(screen) != 203 or screen["rank"].tolist() != list(range(1, 204)):
        raise RuntimeError("Table S2 requires the complete canonical 203-system ranking")
    table_s2 = screen[["rank", "contaminant"]].rename(
        columns={"rank": "Rank", "contaminant": "Contaminant"}
    )
    table_s2["Contaminant"] = table_s2.Contaminant.str.replace("_", " ", regex=False).str.strip()
    table_s2.to_csv(OUTPUT/"Table_S2_contaminant_ranking.csv", index=False)
    print("Generated submitted Tables S1 and S2")


if __name__ == "__main__":
    main()
