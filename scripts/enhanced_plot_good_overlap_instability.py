#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Enhanced plots for Extension 03:
good-overlap denominator instability.

This cleaned version keeps only the two diagnostic figures that carry the main
message of the extension:
1. the observed-weight variance denominator line plot;
2. the epsilon sensitivity curve for the stabilized diagnostic.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


FIGSIZE = (6.2, 3.8)
DPI = 400
PAD = 0.04


def read_csv(input_dir: Path, name: str) -> pd.DataFrame:
    path = input_dir / name
    if not path.exists():
        raise FileNotFoundError(f"Missing required file: {path}")
    return pd.read_csv(path)


def save(fig: plt.Figure, output_dir: Path, name: str) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_dir / name, dpi=DPI, bbox_inches="tight", pad_inches=PAD)
    plt.close(fig)


def fmt_cell(value: float) -> str:
    """Readable labels for small denominators; avoids many misleading 0.00 labels."""
    if not np.isfinite(value):
        return "NA"
    if abs(value) < 0.005:
        return "<0.005"
    return f"{value:.3f}" if abs(value) < 0.1 else f"{value:.2f}"


def plot_denominator_lineplot(alpha_hidden: pd.DataFrame, output_dir: Path) -> None:
    """Plot one representative denominator curve for a clearer diagnostic figure."""
    fig, ax = plt.subplots(figsize=FIGSIZE)

    # The hidden-strength lines almost perfectly overlap in this diagnostic.
    # Keep one representative curve to avoid visual clutter.
    target_hidden = 0.12
    available = np.array(sorted(alpha_hidden["hidden_strength"].dropna().unique()), dtype=float)
    chosen_hidden = available[np.argmin(np.abs(available - target_hidden))]

    d = (
        alpha_hidden.loc[np.isclose(alpha_hidden["hidden_strength"], chosen_hidden)]
        .sort_values("observed_strength")
        .copy()
    )

    ax.plot(
        d["observed_strength"],
        d["var_observed_control_median"],
        marker="o",
        linewidth=2.4,
        markersize=5.5,
    )

    ax.set_xlabel("Observed propensity signal (alpha)")
    ax.set_ylabel("Median Var(observed control weights)")
    ax.set_title("VBM denominator is near zero under good overlap")
    ax.grid(True, linewidth=0.5, alpha=0.35)

    ax.text(
        0.03,
        0.92,
        f"Representative hidden-confounder strength = {chosen_hidden:g}",
        transform=ax.transAxes,
        fontsize=8.5,
        ha="left",
        va="top",
        bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "0.8", "alpha": 0.9},
    )

    save(fig, output_dir, "enhanced_good_overlap_denominator_lineplot_single_curve.png")


def plot_epsilon_sensitivity(epsilon: pd.DataFrame, output_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=FIGSIZE)

    for eps, d in epsilon.groupby("epsilon_fraction"):
        d = d.sort_values("observed_strength")
        ax.plot(
            d["observed_strength"],
            d["stabilized_variance_ratio_median"],
            marker="o",
            linewidth=1.8,
            markersize=4,
            label=f"epsilon={eps:g}",
        )

    ax.axhline(1, linestyle="--", linewidth=1, color="0.45")
    ax.set_yscale("log")
    ax.set_xlabel("Observed propensity signal (alpha)")
    ax.set_ylabel("Median stabilized variance ratio, log scale")
    ax.set_title("Stabilization is dominated by epsilon when alpha is weak")
    ax.legend(title="Stabilizer", fontsize=8, title_fontsize=8, frameon=True)
    ax.text(
        0.01,
        1.5,
        "ratio = 1 means no artificial inflation",
        fontsize=8,
        ha="left",
        va="bottom",
    )

    save(fig, output_dir, "enhanced_good_overlap_epsilon_sensitivity.png")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    input_dir = Path(args.input)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    alpha_hidden = read_csv(input_dir, "good_overlap_alpha_hidden_summary.csv")
    epsilon = read_csv(input_dir, "good_overlap_epsilon_summary.csv")

    plot_denominator_lineplot(alpha_hidden, output_dir)
    plot_epsilon_sensitivity(epsilon, output_dir)

    print(f"Saved cleaned Extension 03 plots to {output_dir}")


if __name__ == "__main__":
    main()
