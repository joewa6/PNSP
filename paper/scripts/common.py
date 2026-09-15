"""Shared paths and plotting helpers for the paper-reproduction scripts."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


PAPER = Path(__file__).resolve().parents[1]
REPOSITORY = PAPER.parent
DATA = PAPER / "processed_data"
OUTPUT = PAPER / "generated"


def configure_style() -> None:
    plt.rcParams.update(
        {
            "figure.dpi": 120,
            "savefig.dpi": 300,
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.alpha": 0.2,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def save(figure: plt.Figure, stem: str) -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT / f"{stem}.png", dpi=300, bbox_inches="tight")
    figure.savefig(OUTPUT / f"{stem}.pdf", bbox_inches="tight")
    plt.close(figure)


def bootstrap_mean(values, repeats: int = 5000, seed: int = 42):
    values = np.asarray(values, dtype=float)
    rng = np.random.default_rng(seed)
    draws = values[rng.integers(0, len(values), size=(repeats, len(values)))].mean(1)
    return float(values.mean()), float(np.percentile(draws, 5)), float(np.percentile(draws, 95))


def panel_label(axis, label: str) -> None:
    axis.text(-0.12, 1.05, label, transform=axis.transAxes, fontsize=12, fontweight="bold")
