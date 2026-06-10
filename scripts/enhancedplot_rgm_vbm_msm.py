#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Enhanced plots for Version10 Extension 04:
RGM model exploration and VBM/MSM/RGM comparison.
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

    thresholds = read_csv(input_dir, "rgm_threshold_summary.csv")
    method = read_csv(input_dir, "rgm_method_comparison.csv")
    rgm_curve = read_csv(input_dir, "plot_rgm_sharp_bootstrap_curve.csv")
    benchmark = read_csv(input_dir, "plot_covariate_benchmark_vbm_msm_rgm.csv")

    fig, ax = plt.subplots(figsize=(8.5, 5))
    finite_thresholds = thresholds[np.isfinite(thresholds["value"])]
    ax.bar(finite_thresholds["quantity"], finite_thresholds["value"])
    ax.set_title("Extension 04 threshold summary")
    ax.set_ylabel("Value")
    ax.tick_params(axis="x", rotation=20)
    save(fig, output_dir, "enhanced_rgm_threshold_summary.png")

    if {"T", "lower", "upper"}.issubset(rgm_curve.columns):
        d = rgm_curve.sort_values("T")
        fig, ax = plt.subplots(figsize=(8.5, 5.5))
        ax.fill_between(d["T"].to_numpy(), d["lower"].to_numpy(), d["upper"].to_numpy(), alpha=0.35)
        ax.plot(d["T"], d["lower"], label="Lower")
        ax.plot(d["T"], d["upper"], label="Upper")
        t_star = thresholds.loc[thresholds["quantity"] == "RGM_sharp_T_star", "value"]
        if len(t_star) and np.isfinite(float(t_star.iloc[0])):
            ax.axvline(float(t_star.iloc[0]), linestyle=":", label="T*")
        ax.axhline(0, linestyle="--", linewidth=1)
        ax.set_xlabel("RGM total-variation radius T")
        ax.set_ylabel("Bootstrap ATT interval")
        ax.set_title("Sharp RGM bootstrap curve")
        ax.legend()
        save(fig, output_dir, "enhanced_rgm_sharp_curve.png")

    usable = benchmark.copy()
    if {"variable", "method", "lower", "upper"}.issubset(usable.columns):
        method_order = ["MSM", "VBM", "RGM-sharp", "RGM-conservative"]
        usable["method_simple"] = usable["method"].astype(str)
        for prefix in ["MSM", "VBM", "RGM-sharp", "RGM-conservative"]:
            usable.loc[usable["method"].astype(str).str.contains(prefix, regex=False), "method_simple"] = prefix
        variables = list(dict.fromkeys(usable["variable"].astype(str)))
        fig, ax = plt.subplots(figsize=(11, 6))
        offsets = {m: (i - 1.5) * 0.16 for i, m in enumerate(method_order)}
        xbase = np.arange(len(variables))
        for m in method_order:
            d = usable[usable["method_simple"] == m]
            if d.empty:
                continue
            xpos = np.array([variables.index(v) for v in d["variable"].astype(str)]) + offsets[m]
            ymid = (d["lower"].to_numpy() + d["upper"].to_numpy()) / 2
            yerr = np.vstack([ymid - d["lower"].to_numpy(), d["upper"].to_numpy() - ymid])
            ax.errorbar(xpos, ymid, yerr=yerr, fmt="o", label=m, capsize=2)
        ax.axhline(0, linestyle="--", linewidth=1)
        ax.set_xticks(xbase)
        ax.set_xticklabels(variables, rotation=35, ha="right")
        ax.set_ylabel("ATT interval")
        ax.set_title("Covariate benchmark comparison across VBM, MSM, and RGM")
        ax.legend()
        save(fig, output_dir, "enhanced_rgm_vbm_msm_benchmark_comparison.png")

    with (output_dir / "enhanced_rgm_method_comparison.md").open("w", encoding="utf-8") as f:
        f.write("# RGM / MSM / VBM comparison\n\n")
        f.write("| method | sensitivity_parameter | geometry | main_assumption | strength | weakness |\n")
        f.write("|---|---|---|---|---|---|\n")
        for _, row in method.iterrows():
            vals = [str(row.get(col, "")).replace("|", "/") for col in ["method", "sensitivity_parameter", "geometry", "main_assumption", "strength", "weakness"]]
            f.write("| " + " | ".join(vals) + " |\n")

    print(f"Saved enhanced Extension 04 plots to {output_dir}")


if __name__ == "__main__":
    main()
