#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Enhanced plotting script for the main NHANES VBM/MSM workflow.

Layout rules:
1. Figure 1 benchmark labels are centered above their vertical reference lines.
2. Closely spaced benchmark labels are separated vertically.
3. Legends are placed below the plotting area.
4. Output PNG files only.
5. Existing image files in the enhanced output folder are cleared before saving.

Usage:
    python scripts/enhanced_plots.py \
        --input output/tables \
        --output output/figures_enhanced

In version10.2 this script is launched automatically by main.R when
config$run_main_enhanced_plots is TRUE.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# Keep SVG text editable
plt.rcParams["font.size"] = 9


# ============================================================
# Utilities
# ============================================================
def required_files_present(input_dir: Path, names: list[str], figure_label: str) -> bool:
    missing = [name for name in names if not (input_dir / name).exists()]
    if missing:
        print(f"Skipping {figure_label}; missing required file(s): {', '.join(missing)}")
        return False
    return True

def clear_output_dir(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    exts = {
        ".png",
        ".jpg",
        ".jpeg",
        ".svg",
        ".pdf",
        ".html",
        ".webp",
    }

    for p in output_dir.iterdir():
        if p.is_file() and p.suffix.lower() in exts:
            p.unlink()


def read_csv(input_dir: Path, name: str) -> pd.DataFrame:
    path = input_dir / name
    if not path.exists():
        raise FileNotFoundError(f"Missing required file: {path}")
    return pd.read_csv(path)



def metric_value(final_results: pd.DataFrame, metric: str, default=np.nan) -> float:
    if "metric" not in final_results.columns or "current_version" not in final_results.columns:
        return default

    s = final_results.loc[final_results["metric"].eq(metric), "current_version"]

    if s.empty:
        return default

    try:
        return float(s.iloc[0])
    except Exception:
        return default


def interpolate_curve(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    xs: np.ndarray,
) -> np.ndarray:
    d = df[[x_col, y_col]].dropna().sort_values(x_col)
    return np.interp(xs, d[x_col].to_numpy(), d[y_col].to_numpy())


def save_figures(fig: plt.Figure, output_dir: Path, stem: str) -> None:
    """Save a figure as PNG only."""
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / f"{stem}.png", format="png", dpi=300, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)

def assign_benchmark_label_heights(bmk: pd.DataFrame) -> pd.DataFrame:
    """
    Keep each label directly above its benchmark vertical line.

    Only y_text is changed to avoid overlap.
    x_text is always equal to sensitivity_value.
    """
    bmk = bmk.sort_values("sensitivity_value").reset_index(drop=True).copy()

    if bmk.empty:
        bmk["x_text"] = []
        bmk["y_text"] = []
        return bmk

    x = bmk["sensitivity_value"].to_numpy(dtype=float)

    bmk["x_text"] = x
    bmk["y_text"] = 0.50

    # Treat nearby x-values as a cluster.
    # Because labels are not horizontally shifted, close labels are placed on different rows.
    cluster_threshold = 0.045

    clusters = []
    current = [0]

    for i in range(1, len(x)):
        if abs(x[i] - x[i - 1]) < cluster_threshold:
            current.append(i)
        else:
            clusters.append(current)
            current = [i]

    clusters.append(current)

    # More levels = less vertical collision.
    # These values are in the top annotation band coordinate system, from 0 to 1.
    y_levels_for_cluster = [
        0.18,
        0.30,
        0.42,
        0.54,
        0.66,
    ]

    single_label_y = 0.42

    for cluster in clusters:
        if len(cluster) == 1:
            idx = cluster[0]
            bmk.loc[idx, "x_text"] = x[idx]
            bmk.loc[idx, "y_text"] = single_label_y
        else:
            for j, idx in enumerate(cluster):
                bmk.loc[idx, "x_text"] = x[idx]
                bmk.loc[idx, "y_text"] = y_levels_for_cluster[j % len(y_levels_for_cluster)]

    return bmk


# ============================================================
# Figure 1
# ============================================================
def make_figure1(input_dir: Path, output_dir: Path) -> None:
    if not required_files_present(input_dir, [
        "plot_vbm_bootstrap_curve.csv",
        "final_results_vs_paper.csv",
        "plot_covariate_benchmark_vbm_msm.csv",
    ], "main enhanced Figure 1"):
        return

    bootstrap = read_csv(input_dir, "plot_vbm_bootstrap_curve.csv")
    final_results = read_csv(input_dir, "final_results_vs_paper.csv")
    benchmark = read_csv(input_dir, "plot_covariate_benchmark_vbm_msm.csv")

    xs = np.round(np.arange(0.0, 0.71, 0.1), 2)

    boot = bootstrap.copy()
    if "stage" in boot.columns:
        coarse = boot.loc[boot["stage"].astype(str).eq("coarse")]
        if not coarse.empty:
            boot = coarse

    ci_lower = interpolate_curve(boot, "R2", "lower", xs)
    ci_upper = interpolate_curve(boot, "R2", "upper", xs)

    # Deterministic point-bound curves are intentionally absent from the
    # current pipeline. Figure 1 therefore reports only bootstrap intervals
    # and their midpoints, with no duplicated interval-bound overlay.
    midpoints = (ci_lower + ci_upper) / 2

    r2_star = metric_value(
        final_results,
        "VBM_R2_star_bootstrap",
        default=np.nan,
    )

    label_map = {
        "gender": "Gender",
        "age": "Age",
        "income": "Income",
        "income_missing": "Income (Missing)",
        "education": "Education",
        "cig_smoked": "Cig. Smoked",
        "smoking_history": "Smoking history",
        "race": "Race",
    }

    bmk = benchmark.loc[
        benchmark["method"].astype(str).eq("VBM benchmark (bootstrap)"),
        ["variable", "sensitivity_value"],
    ].dropna().drop_duplicates()

    bmk["label"] = bmk["variable"].map(label_map).fillna(bmk["variable"].astype(str))

    # Wider x range leaves room at both sides.
    xlim = (-0.08, 0.78)

    bmk = bmk.loc[
        (bmk["sensitivity_value"] >= xlim[0])
        & (bmk["sensitivity_value"] <= xlim[1])
    ].copy()

    bmk = assign_benchmark_label_heights(bmk)

    # Large figure:
    # top band for benchmark labels,
    # main panel for curve,
    # bottom room for legend.
    fig = plt.figure(figsize=(16, 7.0))

    gs = fig.add_gridspec(
        nrows=2,
        ncols=1,
        height_ratios=[2.6, 7.8],
        hspace=0.14,
    )

    ax_top = fig.add_subplot(gs[0, 0])
    ax = fig.add_subplot(gs[1, 0], sharex=ax_top)

    # ---------------- Main plot ----------------
    ax.vlines(
        xs,
        ci_lower,
        ci_upper,
        linewidth=3,
        color="#1f77b4",
        alpha=0.95,
        label="95% bootstrap CI",
        zorder=1,
    )

    ax.scatter(
        xs,
        midpoints,
        s=58,
        label="Bootstrap interval midpoint",
        zorder=3,
    )

    ax.axhline(
        0,
        linewidth=2.0,
        linestyle="--",
        label="Null effect",
        zorder=0,
    )

    if np.isfinite(r2_star):
        ax.axvline(
            r2_star,
            linewidth=2.2,
            linestyle=":",
            zorder=0,
        )

    for _, row in bmk.iterrows():
        ax.axvline(
            float(row["sensitivity_value"]),
            linewidth=1.1,
            alpha=0.25,
            zorder=0,
        )

    ymin = min(ci_lower.min(), -0.15)
    ymax = max(ci_upper.max(), 4.7)

    ax.set_xlim(*xlim)
    ax.set_ylim(ymin - 0.12, ymax + 0.45)

    # The title is drawn in the top annotation band to avoid overlap
    # with benchmark labels when the figure is flattened.

    ax.set_xlabel("R²", fontsize=11)
    ax.set_ylabel("Estimated ATT", fontsize=11)
    ax.tick_params(axis="both", labelsize=9)

    # Put legend below the plot, away from the title and top labels.
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.12),
        ncol=3,
        frameon=False,
        fontsize=8.5,
        handlelength=2.0,
        columnspacing=2.1,
    )

    # ---------------- Top annotation band ----------------
    ax_top.set_xlim(*xlim)
    ax_top.set_ylim(0, 1)
    ax_top.axis("off")

    ax_top.text(
        0.5,
        0.98,
        "Figure 1-style VBM sensitivity curve",
        transform=ax_top.transAxes,
        ha="center",
        va="top",
        fontsize=13,
        fontweight="bold",
    )

    ax_top.text(
        xlim[0] + 0.005,
        0.76,
        "Benchmarked covariates",
        ha="left",
        va="top",
        fontsize=8.8,
        fontweight="bold",
    )

    # R² star label
    if np.isfinite(r2_star):
        ax_top.axvline(
            r2_star,
            ymin=0.05,
            ymax=0.95,
            linewidth=2.2,
            linestyle=":",
        )

        ax_top.text(
            r2_star + 0.006,
            0.55,
            f"R²* = {r2_star:.2f}",
            rotation=90,
            va="center",
            ha="left",
            fontsize=9.5,
        )

    # Benchmark stems and labels.
    # Label x-position equals vertical-line x-position.
    for _, row in bmk.iterrows():
        x = float(row["sensitivity_value"])
        y_text = float(row["y_text"])
        label = str(row["label"])

        # short stem from bottom of annotation band
        ax_top.plot(
            [x, x],
            [0.02, 0.16],
            linewidth=1.1,
            alpha=0.45,
        )

        # centered exactly above the vertical line
        ax_top.annotate(
            label,
            xy=(x, 0.16),
            xytext=(x, y_text),
            textcoords="data",
            ha="center",
            va="center",
            fontsize=8.2,
            arrowprops=dict(
                arrowstyle="-",
                lw=0.9,
                alpha=0.6,
                shrinkA=4,
                shrinkB=2,
            ),
            clip_on=False,
        )

    # Leave enough bottom space for legend.
    fig.subplots_adjust(bottom=0.20)

    save_figures(fig, output_dir, "figure1_vbm_paper_style")


# ============================================================
# Figure 3
# ============================================================
def make_figure3(input_dir: Path, output_dir: Path) -> None:
    if not required_files_present(input_dir, ["plot_covariate_benchmark_vbm_msm.csv"], "main enhanced Figure 3"):
        return

    benchmark = read_csv(input_dir, "plot_covariate_benchmark_vbm_msm.csv")

    label_map = {
        "gender": "Gender",
        "age": "Age",
        "income": "Income",
        "income_missing": "Income (Missing)",
        "education": "Education",
        "cig_smoked": "Cig. Smoked",
        "smoking_history": "Smoking history",
        "race": "Race",
    }

    order = list(label_map.keys())

    df = benchmark.copy()

    df = df[
        df["method"].astype(str).str.contains(
            "VBM|MSM|Qbal|Corr",
            case=False,
            regex=True,
            na=False,
        )
    ].copy()

    def normalize_method(x: str) -> str:
        x = str(x)
        xu = x.upper()

        if "QBAL" in xu:
            return "MSM (Qbal)"

        if "CORR" in xu:
            return "VBM, w/ Corr."

        if "MSM" in xu:
            return "MSM"

        if "VBM" in xu:
            return "VBM"

        return x

    df["method_label"] = df["method"].map(normalize_method)

    df["x_base"] = df["variable"].apply(
        lambda v: order.index(v) if v in order else len(order)
    )

    method_order = ["MSM", "MSM (Qbal)", "VBM", "VBM, w/ Corr."]
    methods = [m for m in method_order if m in df["method_label"].unique()]

    offsets = {
        "MSM": -0.30,
        "MSM (Qbal)": -0.10,
        "VBM": 0.10,
        "VBM, w/ Corr.": 0.30,
    }

    fig, ax = plt.subplots(figsize=(16, 5.5))

    for method in methods:
        sub = df.loc[df["method_label"].eq(method)].sort_values("x_base")

        x = sub["x_base"].to_numpy(dtype=float) + offsets.get(method, 0.0)
        y = sub["tau_hat"].to_numpy(dtype=float)
        lower = sub["lower"].to_numpy(dtype=float)
        upper = sub["upper"].to_numpy(dtype=float)

        yerr = np.vstack(
            [
                y - lower,
                upper - y,
            ]
        )

        ax.errorbar(
            x,
            y,
            yerr=yerr,
            fmt="o",
            capsize=3,
            linewidth=1.6,
            markersize=4.8,
            label=method,
        )

    tau = pd.to_numeric(df["tau_hat"], errors="coerce").dropna()

    if not tau.empty:
        ax.axhline(
            float(tau.iloc[0]),
            linewidth=2.0,
            linestyle=":",
            label="Original ATT",
        )

    ax.axhline(
        0,
        linewidth=2.0,
        linestyle="--",
        label="Null effect",
    )

    ax.set_xticks(range(len(order)))

    ax.set_xticklabels(
        [label_map[v] for v in order],
        rotation=12,
        ha="right",
        fontsize=8.8,
    )

    ax.tick_params(axis="y", labelsize=9)

    ax.set_ylabel("Estimated ATT", fontsize=11)
    ax.set_xlabel("")

    ax.set_title(
        "Figure 3-style benchmark intervals: four model comparison",
        fontsize=15,
        pad=10,
    )

    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.13),
        ncol=6,
        frameon=False,
        fontsize=8.5,
        handlelength=1.8,
        columnspacing=2.1,
    )

    fig.subplots_adjust(bottom=0.25)

    save_figures(fig, output_dir, "figure3_benchmark_paper_style")



# ============================================================
# Optional RGM benchmark extension figure
# ============================================================
def make_figure4_rgm_extension(input_dir: Path, output_dir: Path) -> None:
    path = input_dir / "plot_covariate_benchmark_vbm_msm_rgm.csv"
    if not path.exists():
        return

    benchmark = pd.read_csv(path)

    label_map = {
        "gender": "Gender",
        "age": "Age",
        "income": "Income",
        "income_missing": "Income (Missing)",
        "education": "Education",
        "cig_smoked": "Cig. Smoked",
        "smoking_history": "Smoking history",
        "race": "Race",
    }

    order = [v for v in label_map if v in set(benchmark["variable"].astype(str))]
    remaining = [v for v in benchmark["variable"].astype(str).unique() if v not in order]
    order = order + remaining

    def normalize_method(x: str) -> str:
        x = str(x)
        xu = x.upper()
        if "RGM-SHARP" in xu or "RGM (SHARP" in xu:
            return "RGM (sharp)"
        if "RGM-CONSERVATIVE" in xu or xu == "RGM":
            return "RGM"
        if xu == "MSM" or "MSM BENCHMARK" in xu:
            return "MSM"
        if xu == "VBM" or "VBM BENCHMARK" in xu:
            return "VBM"
        return x

    df = benchmark.copy()
    df["method_label"] = df["method"].map(normalize_method)
    df = df[df["method_label"].isin(["MSM", "VBM", "RGM", "RGM (sharp)"])].copy()

    if df.empty:
        return

    df["x_base"] = df["variable"].astype(str).apply(
        lambda v: order.index(v) if v in order else len(order)
    )

    method_order = ["MSM", "VBM", "RGM", "RGM (sharp)"]
    methods = [m for m in method_order if m in df["method_label"].unique()]
    offsets = {
        "MSM": -0.27,
        "VBM": -0.09,
        "RGM": 0.09,
        "RGM (sharp)": 0.27,
    }

    fig, ax = plt.subplots(figsize=(16, 5.5))

    for method in methods:
        sub = df.loc[df["method_label"].eq(method)].sort_values("x_base")
        x = sub["x_base"].to_numpy(dtype=float) + offsets.get(method, 0.0)
        y = pd.to_numeric(sub["tau_hat"], errors="coerce").to_numpy(dtype=float)
        lower = pd.to_numeric(sub["lower"], errors="coerce").to_numpy(dtype=float)
        upper = pd.to_numeric(sub["upper"], errors="coerce").to_numpy(dtype=float)
        yerr = np.vstack([y - lower, upper - y])

        ax.errorbar(
            x,
            y,
            yerr=yerr,
            fmt="o",
            capsize=5,
            linewidth=2.3,
            markersize=7.5,
            label=method,
        )

    tau = pd.to_numeric(df["tau_hat"], errors="coerce").dropna()
    if not tau.empty:
        ax.axhline(float(tau.iloc[0]), linewidth=2.0, linestyle=":", label="Original ATT")

    ax.axhline(0, linewidth=2.0, linestyle="--", label="Null effect")

    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(
        [label_map.get(v, v) for v in order],
        rotation=10,
        ha="right",
        fontsize=15,
    )

    ax.tick_params(axis="y", labelsize=15)
    ax.set_ylabel("Estimated ATT", fontsize=11)
    ax.set_xlabel("Benchmark covariate", fontsize=17)
    ax.set_title(
        "Benchmark intervals with RGM extension",
        fontsize=24,
        pad=22,
    )

    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.10),
        ncol=6,
        frameon=False,
        fontsize=14,
        handlelength=2.3,
        columnspacing=2.1,
    )

    fig.subplots_adjust(bottom=0.25)
    save_figures(fig, output_dir, "figure4_benchmark_rgm_extension_paper_style")

# ============================================================
# Main
# ============================================================
def main() -> None:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input",
        type=Path,
        default=Path("output/tables"),
        help="Directory containing CSV files.",
    )

    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/figures_enhanced"),
        help="Directory to save figures.",
    )

    args = parser.parse_args()

    clear_output_dir(args.output)

    make_figure1(args.input, args.output)
    make_figure3(args.input, args.output)
    make_figure4_rgm_extension(args.input, args.output)

    print(f"Done. PNG figures saved to: {args.output.resolve()}")


if __name__ == "__main__":
    main()
