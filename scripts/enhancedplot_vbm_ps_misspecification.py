#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Enhanced plots for Version10 Extension 02:
VBM under observed propensity-score misspecification.
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    input_dir = Path(args.input)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    summary = read_csv(input_dir, "vbm_ps_misspec_factorial_summary.csv")
    decomp = read_csv(input_dir, "vbm_ps_misspec_decomposition.csv")
    curves = read_csv(input_dir, "vbm_ps_misspec_bootstrap_curves.csv")

    order = ["no_hidden_correct", "no_hidden_misspecified", "hidden_correct", "hidden_misspecified"]
    labels = ["A1\nNeither", "A2\nMisspec only", "B1\nHidden only", "B2\nBoth"]
    summary["scenario"] = pd.Categorical(summary["scenario"], order, ordered=True)
    summary = summary.sort_values("scenario")

    fig, ax = plt.subplots(figsize=(9.5, 5.5))
    x = np.arange(len(summary))
    ax.bar(x, summary["combined_error_R2"].to_numpy(), label="Combined error R²")
    ax.scatter(x, summary["VBM_R2_star"].to_numpy(), marker="D", label="VBM R²*")
    ax.scatter(x, summary["hidden_only_R2"].to_numpy(), marker="o", label="Hidden-only R²")
    ax.scatter(x, summary["misspec_only_R2"].to_numpy(), marker="s", label="Misspec-only R²")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("R² scale")
    ax.set_title("2 x 2 decomposition: hidden confounding versus PS misspecification")
    ax.legend()
    save(fig, output_dir, "enhanced_vbm_ps_misspec_decomposition.png")

    fig, ax = plt.subplots(figsize=(9.5, 5.5))
    ax.axhline(0, linestyle="--", linewidth=1)
    ax.bar(x, summary["observed_minus_oracle_ATT"].to_numpy())
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Observed analysis ATT - oracle ATT")
    ax.set_title("ATT distortion induced by misspecification and hidden U")
    save(fig, output_dir, "enhanced_vbm_ps_misspec_att_gap.png")

    pivot = summary.pivot(index="hidden_present", columns="ps_setting", values="combined_error_R2")
    fig, ax = plt.subplots(figsize=(6.8, 4.8))
    im = ax.imshow(pivot.values, aspect="auto")
    ax.set_xticks(np.arange(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns)
    ax.set_yticks(np.arange(len(pivot.index)))
    ax.set_yticklabels(["Hidden U absent" if not v else "Hidden U present" for v in pivot.index])
    ax.set_title("Combined error R² matrix")
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            val = pivot.values[i, j]
            ax.text(j, i, f"{val:.3f}" if np.isfinite(val) else "NA", ha="center", va="center")
    fig.colorbar(im, ax=ax, label="Combined error R²")
    save(fig, output_dir, "enhanced_vbm_ps_misspec_matrix.png")

    if {"stage", "R2", "lower", "upper", "scenario"}.issubset(curves.columns):
        fig, ax = plt.subplots(figsize=(10, 5.8))
        for scenario, d in curves.groupby("scenario"):
            d = d.sort_values("R2")
            fine = d[d["stage"].astype(str).str.contains("fine", case=False, na=False)]
            if fine.empty:
                fine = d
            ax.plot(fine["R2"], fine["lower"], label=f"{scenario} lower")
            ax.plot(fine["R2"], fine["upper"], linestyle="--", label=f"{scenario} upper")
        ax.axhline(0, linestyle=":", linewidth=1)
        ax.set_title("VBM bootstrap lower/upper curves by factorial cell")
        ax.set_xlabel("R²")
        ax.set_ylabel("Bootstrap ATT interval")
        ax.legend(fontsize=8, ncol=2)
        save(fig, output_dir, "enhanced_vbm_ps_misspec_bootstrap_curves.png")

    decomp.to_csv(output_dir / "enhanced_vbm_ps_misspec_readme_table.csv", index=False)
    print(f"Saved enhanced Extension 02 plots to {output_dir}")


if __name__ == "__main__":
    main()
