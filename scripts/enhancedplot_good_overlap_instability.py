#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Enhanced plots for Version10 Extension 03:
good-overlap denominator instability.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def read_csv(input_dir: Path, name: str) -> pd.DataFrame:
    path = input_dir / name
    if not path.exists():
        raise FileNotFoundError(f"Missing required file: {path}")
    return pd.read_csv(path)


def save(fig: plt.Figure, output_dir: Path, name: str) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_dir / name, dpi=300, bbox_inches="tight")
    plt.close(fig)


def heatmap(df: pd.DataFrame, index: str, columns: str, values: str, title: str, output_dir: Path, name: str) -> None:
    pivot = df.pivot_table(index=index, columns=columns, values=values, aggfunc="mean").sort_index()
    fig, ax = plt.subplots(figsize=(9.2, 5.5))
    im = ax.imshow(pivot.values, aspect="auto")
    ax.set_xticks(np.arange(pivot.shape[1]))
    ax.set_xticklabels([str(x) for x in pivot.columns], rotation=45, ha="right")
    ax.set_yticks(np.arange(pivot.shape[0]))
    ax.set_yticklabels([str(x) for x in pivot.index])
    ax.set_title(title)
    ax.set_xlabel(columns)
    ax.set_ylabel(index)
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            val = pivot.values[i, j]
            txt = f"{val:.2f}" if np.isfinite(val) else "NA"
            ax.text(j, i, txt, ha="center", va="center", fontsize=8)
    fig.colorbar(im, ax=ax, label=values)
    save(fig, output_dir, name)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    input_dir = Path(args.input)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    alpha_hidden = read_csv(input_dir, "good_overlap_alpha_hidden_summary.csv")
    sample_size = read_csv(input_dir, "good_overlap_sample_size_summary.csv")
    epsilon = read_csv(input_dir, "good_overlap_epsilon_summary.csv")
    protocol = read_csv(input_dir, "good_overlap_diagnostic_protocol.csv")

    heatmap(
        alpha_hidden,
        index="hidden_strength",
        columns="observed_strength",
        values="flag_rate",
        title="Good-overlap instability flag rate across alpha and hidden strength",
        output_dir=output_dir,
        name="enhanced_good_overlap_flag_heatmap.png",
    )

    heatmap(
        alpha_hidden,
        index="hidden_strength",
        columns="observed_strength",
        values="var_observed_control_median",
        title="Median observed-weight variance denominator",
        output_dir=output_dir,
        name="enhanced_good_overlap_denominator_heatmap.png",
    )

    fig, ax = plt.subplots(figsize=(9.5, 5.5))
    for n, d in sample_size.groupby("sample_size"):
        d = d.sort_values("observed_strength")
        ax.plot(d["observed_strength"], d["flag_rate"], marker="o", label=f"n={n}")
    ax.set_xlabel("Observed propensity strength alpha")
    ax.set_ylabel("Instability flag rate")
    ax.set_title("Sample-size robustness of good-overlap denominator risk")
    ax.legend()
    save(fig, output_dir, "enhanced_good_overlap_sample_size_flag_rate.png")

    fig, ax = plt.subplots(figsize=(9.5, 5.5))
    for eps, d in epsilon.groupby("epsilon_fraction"):
        d = d.sort_values("observed_strength")
        ax.plot(d["observed_strength"], d["stabilized_variance_ratio_median"], marker="o", label=f"epsilon={eps:g}")
    ax.set_yscale("log")
    ax.set_xlabel("Observed propensity strength alpha")
    ax.set_ylabel("Median stabilized variance ratio, log scale")
    ax.set_title("Epsilon sensitivity for the stabilized diagnostic")
    ax.legend()
    save(fig, output_dir, "enhanced_good_overlap_epsilon_sensitivity.png")

    with (output_dir / "enhanced_good_overlap_protocol.txt").open("w", encoding="utf-8") as f:
        f.write("Version10 good-overlap diagnostic protocol\n")
        f.write("=" * 48 + "\n")
        for _, row in protocol.iterrows():
            f.write(f"{int(row['step'])}. {row['diagnostic']} [{row['version10_rule']}]\n")

    print(f"Saved enhanced Extension 03 plots to {output_dir}")


if __name__ == "__main__":
    main()
