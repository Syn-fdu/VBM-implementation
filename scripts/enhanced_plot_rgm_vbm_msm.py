#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Enhanced plots for Extension 04:
RGM model exploration and VBM/MSM/RGM comparison.

This script produces report-ready figures:
1. Flatter layout.
2. Tighter within-covariate spacing.
3. Calmer custom colors.
4. Thicker interval lines and caps.
5. Cleaner legends, grids, and labels.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# Global style
# ============================================================
plt.rcParams.update({
    "font.size": 10,
    "axes.titlesize": 15,
    "axes.labelsize": 11,
    "xtick.labelsize": 9.5,
    "ytick.labelsize": 9.5,
    "legend.fontsize": 9.5,
    "axes.linewidth": 0.9,
    "savefig.dpi": 300,
})

PALETTE = {
    # Muted, paper-like palette aligned with the main benchmark figures.
    "MSM": "#D9992B",              # muted ochre
    "VBM": "#4C78A8",              # muted blue
    "RGM-sharp": "#3F3B7A",        # muted indigo
    "RGM-conservative": "#B65C5A", # muted brick red
    "lower": "#6F7FAF",
    "upper": "#6F7FAF",
    "fill": "#B7C9D9",
    "zero": "#4D4D4D",
    "tstar": "#4C78A8",
}

LABEL_MAP = {
    "gender": "Gender",
    "age": "Age",
    "income": "Income",
    "income_missing": "Income missing",
    "education": "Education",
    "cig_smoked": "Cig. smoked",
    "smoking_history": "Smoking history",
    "race": "Race",
}

METHOD_LABEL_MAP = {
    "MSM": "MSM",
    "VBM": "VBM",
    "RGM-sharp": "RGM-sharp",
    "RGM-conservative": "RGM-conservative",
}


# ============================================================
# Utilities
# ============================================================
def read_csv(input_dir: Path, name: str) -> pd.DataFrame:
    path = input_dir / name
    if not path.exists():
        raise FileNotFoundError(f"Missing required file: {path}")
    return pd.read_csv(path)


def save(fig: plt.Figure, output_dir: Path, name: str) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / name, dpi=300, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)


def clean_axes(ax: plt.Axes) -> None:
    """Light report-style axes."""
    ax.grid(True, axis="y", linewidth=0.7, alpha=0.25)
    ax.grid(False, axis="x")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def normalize_method_name(x: str) -> str:
    x = str(x)
    for prefix in ["RGM-conservative", "RGM-sharp", "MSM", "VBM"]:
        if prefix in x:
            return prefix
    return x


# ============================================================
# Plot 1: threshold summary
# ============================================================
def plot_threshold_summary(thresholds: pd.DataFrame, output_dir: Path) -> None:
    finite_thresholds = thresholds[np.isfinite(thresholds["value"])].copy()
    if finite_thresholds.empty:
        return

    fig, ax = plt.subplots(figsize=(8.2, 3.8))

    x = np.arange(len(finite_thresholds))
    bars = ax.bar(
        x,
        finite_thresholds["value"].to_numpy(dtype=float),
        width=0.58,
        color="#4C78A8",
        alpha=0.88,
    )

    ax.set_xticks(x)
    ax.set_xticklabels(
        finite_thresholds["quantity"].astype(str),
        rotation=18,
        ha="right",
    )

    ax.set_ylabel("Value")
    ax.set_title("Extension 04 threshold summary", pad=8)
    clean_axes(ax)

    for bar in bars:
        h = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            h,
            f"{h:.3g}",
            ha="center",
            va="bottom",
            fontsize=8.5,
        )

    fig.subplots_adjust(left=0.10, right=0.98, top=0.85, bottom=0.28)
    save(fig, output_dir, "enhanced_rgm_threshold_summary.png")


# ============================================================
# Plot 2: RGM sharp bootstrap curve
# ============================================================
def plot_rgm_curve(
    rgm_curve: pd.DataFrame,
    thresholds: pd.DataFrame,
    output_dir: Path,
) -> None:
    if not {"T", "lower", "upper"}.issubset(rgm_curve.columns):
        return

    d = rgm_curve.sort_values("T").copy()

    x = d["T"].to_numpy(dtype=float)
    lower = d["lower"].to_numpy(dtype=float)
    upper = d["upper"].to_numpy(dtype=float)

    fig, ax = plt.subplots(figsize=(10.8, 4.3))

    ax.fill_between(
        x,
        lower,
        upper,
        color=PALETTE["fill"],
        alpha=0.55,
        linewidth=0,
        label="Bootstrap interval",
        zorder=1,
    )

    ax.plot(
        x,
        lower,
        color=PALETTE["lower"],
        linewidth=2.4,
        label="Lower bound",
        zorder=3,
    )

    ax.plot(
        x,
        upper,
        color=PALETTE["upper"],
        linewidth=2.4,
        label="Upper bound",
        zorder=3,
    )

    ax.axhline(
        0,
        color=PALETTE["zero"],
        linestyle="--",
        linewidth=1.3,
        alpha=0.8,
        label="Null effect",
        zorder=2,
    )

    t_star = thresholds.loc[thresholds["quantity"] == "RGM_sharp_T_star", "value"]
    if len(t_star) and np.isfinite(float(t_star.iloc[0])):
        t_star_value = float(t_star.iloc[0])
        ax.axvline(
            t_star_value,
            color=PALETTE["tstar"],
            linestyle=":",
            linewidth=2.2,
            label=f"T* = {t_star_value:.3f}",
            zorder=4,
        )
        ax.text(
            t_star_value + 0.012,
            ax.get_ylim()[0] + 0.08 * (ax.get_ylim()[1] - ax.get_ylim()[0]),
            f"T* = {t_star_value:.3f}",
            rotation=90,
            va="bottom",
            ha="left",
            fontsize=9,
            color=PALETTE["tstar"],
        )

    ax.set_xlabel("RGM total-variation radius T")
    ax.set_ylabel("Bootstrap ATT interval")
    ax.set_title("Sharp RGM sensitivity curve", pad=8)

    clean_axes(ax)

    # Slightly more compact and less intrusive legend.
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.16),
        ncol=5,
        frameon=False,
        handlelength=2.2,
        columnspacing=1.5,
    )

    fig.subplots_adjust(left=0.08, right=0.985, top=0.88, bottom=0.27)
    save(fig, output_dir, "enhanced_rgm_sharp_curve.png")


# ============================================================
# Plot 3: benchmark comparison
# ============================================================
def plot_benchmark_comparison(benchmark: pd.DataFrame, output_dir: Path) -> None:
    usable = benchmark.copy()

    needed = {"variable", "method", "lower", "upper"}
    if not needed.issubset(usable.columns):
        return

    method_order = ["MSM", "VBM", "RGM-sharp", "RGM-conservative"]

    usable["method_simple"] = usable["method"].astype(str).map(normalize_method_name)
    usable = usable[usable["method_simple"].isin(method_order)].copy()

    if usable.empty:
        return

    # Prefer paper-like order when present.
    paper_order = [
        "gender",
        "age",
        "income",
        "income_missing",
        "education",
        "cig_smoked",
        "smoking_history",
        "race",
    ]
    existing = list(dict.fromkeys(usable["variable"].astype(str)))
    variables = [v for v in paper_order if v in existing] + [v for v in existing if v not in paper_order]

    # Tighter within-covariate spacing than the coarser T-grid step.
    # Old offsets were approximately [-0.24, -0.08, 0.08, 0.24].
    # New offsets make the four intervals read as one compact cluster.
    offsets = {
        "MSM": -0.135,
        "VBM": -0.045,
        "RGM-sharp": 0.045,
        "RGM-conservative": 0.135,
    }

    xbase = np.arange(len(variables))
    fig, ax = plt.subplots(figsize=(11.6, 4.9))

    for m in method_order:
        d = usable[usable["method_simple"] == m].copy()
        if d.empty:
            continue

        d["xbase"] = d["variable"].astype(str).map({v: i for i, v in enumerate(variables)})
        d = d.sort_values("xbase")

        xpos = d["xbase"].to_numpy(dtype=float) + offsets[m]

        lower = pd.to_numeric(d["lower"], errors="coerce").to_numpy(dtype=float)
        upper = pd.to_numeric(d["upper"], errors="coerce").to_numpy(dtype=float)

        # Use tau_hat when available; otherwise keep midpoint.
        if "tau_hat" in d.columns:
            y = pd.to_numeric(d["tau_hat"], errors="coerce").to_numpy(dtype=float)
            if np.any(~np.isfinite(y)):
                y = (lower + upper) / 2
        else:
            y = (lower + upper) / 2

        yerr = np.vstack([y - lower, upper - y])

        ax.errorbar(
            xpos,
            y,
            yerr=yerr,
            fmt="o",
            color=PALETTE[m],
            ecolor=PALETTE[m],
            elinewidth=2.25,      # thicker vertical intervals
            capsize=4.2,
            capthick=2.0,
            markersize=6.4,
            markeredgewidth=0.7,
            markeredgecolor="white",
            label=METHOD_LABEL_MAP[m],
            zorder=3,
        )

    ax.axhline(
        0,
        color=PALETTE["zero"],
        linestyle="--",
        linewidth=1.3,
        alpha=0.8,
        label="Null effect",
        zorder=1,
    )

    ax.set_xticks(xbase)
    ax.set_xticklabels(
        [LABEL_MAP.get(v, v.replace("_", " ").title()) for v in variables],
        rotation=20,
        ha="right",
    )

    ax.set_xlim(-0.55, len(variables) - 0.45)

    ax.set_ylabel("ATT interval")
    ax.set_title("Covariate benchmark comparison across VBM, MSM, and RGM", pad=8)

    clean_axes(ax)

    # Keep the legend away from data and from x tick labels.
    handles, labels = ax.get_legend_handles_labels()
    order_labels = ["MSM", "VBM", "RGM-sharp", "RGM-conservative", "Null effect"]
    order_idx = [labels.index(lbl) for lbl in order_labels if lbl in labels]
    handles = [handles[i] for i in order_idx]
    labels = [labels[i] for i in order_idx]

    ax.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.19),
        ncol=5,
        frameon=False,
        handlelength=2.0,
        columnspacing=1.4,
    )

    fig.subplots_adjust(left=0.075, right=0.985, top=0.86, bottom=0.31)
    save(fig, output_dir, "enhanced_rgm_vbm_msm_benchmark_comparison.png")


# ============================================================
# Method comparison markdown
# ============================================================
def write_method_table(method: pd.DataFrame, output_dir: Path) -> None:
    with (output_dir / "enhanced_rgm_method_comparison.md").open("w", encoding="utf-8") as f:
        f.write("# RGM / MSM / VBM comparison\n\n")
        f.write("| method | sensitivity_parameter | geometry | main_assumption | strength | weakness |\n")
        f.write("|---|---|---|---|---|---|\n")
        for _, row in method.iterrows():
            vals = [
                str(row.get(col, "")).replace("|", "/")
                for col in [
                    "method",
                    "sensitivity_parameter",
                    "geometry",
                    "main_assumption",
                    "strength",
                    "weakness",
                ]
            ]
            f.write("| " + " | ".join(vals) + " |\n")


# ============================================================
# Main
# ============================================================
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

    plot_threshold_summary(thresholds, output_dir)
    plot_rgm_curve(rgm_curve, thresholds, output_dir)
    plot_benchmark_comparison(benchmark, output_dir)
    write_method_table(method, output_dir)

    print(f"Saved enhanced Extension 04 plots to {output_dir}")


if __name__ == "__main__":
    main()
