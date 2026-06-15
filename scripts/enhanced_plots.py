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
from matplotlib.lines import Line2D


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
    point_bounds_path = input_dir / "plot_vbm_point_bounds_curve.csv"
    point_bounds = read_csv(input_dir, "plot_vbm_point_bounds_curve.csv") if point_bounds_path.exists() else None

    xs = np.round(np.arange(0.0, 0.71, 0.1), 2)

    boot = bootstrap.copy()
    if "stage" in boot.columns:
        coarse = boot.loc[boot["stage"].astype(str).eq("coarse")]
        if not coarse.empty:
            boot = coarse

    ci_lower = interpolate_curve(boot, "R2", "lower", xs)
    ci_upper = interpolate_curve(boot, "R2", "upper", xs)
    midpoints = (ci_lower + ci_upper) / 2

    if point_bounds is not None and {"R2", "lower", "upper"}.issubset(point_bounds.columns):
        point_lower = interpolate_curve(point_bounds, "R2", "lower", xs)
        point_upper = interpolate_curve(point_bounds, "R2", "upper", xs)
        point_midpoints = (point_lower + point_upper) / 2
    else:
        point_lower = ci_lower
        point_upper = ci_upper
        point_midpoints = midpoints

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

    xlim = (-0.08, 0.78)
    bmk = bmk.loc[
        (bmk["sensitivity_value"] >= xlim[0])
        & (bmk["sensitivity_value"] <= xlim[1])
    ].copy()

    # The left-cluster benchmark lines are visually indistinguishable at this scale.
    # Plot them as a single grouped reference line at R2 = 0.
    other_vars = {
        "gender",
        "age",
        "income_missing",
        "cig_smoked",
        "smoking_history",
    }
    left_cluster_benchmark = bmk["variable"].isin(other_vars)
    bmk["plot_x"] = bmk["sensitivity_value"].where(~left_cluster_benchmark, 0.0)

    fig, ax = plt.subplots(figsize=(14.5, 4.7))

    # ---------------- Main plot ----------------
    ci_color = "#8CB6D8"
    point_color = "#1F77B4"
    point_interval_color = "#3C8DC7"
    null_color = "#7C8794"
    benchmark_color = "#A8C6DA"
    label_color = "#536879"

    yerr = np.vstack(
        [
            midpoints - ci_lower,
            ci_upper - midpoints,
        ]
    )

    ax.errorbar(
        xs,
        midpoints,
        yerr=yerr,
        fmt="none",
        ecolor=ci_color,
        elinewidth=1.8,
        capsize=5.0,
        capthick=1.45,
        alpha=0.9,
        label="95% bootstrap CI",
        zorder=2,
    )

    point_interval_label_used = False
    for x, lo, hi in zip(xs, point_lower, point_upper):
        ax.plot(
            [x, x],
            [lo, hi],
            color=point_interval_color,
            linewidth=7.2,
            solid_capstyle="butt",
            alpha=0.92,
            label="ATT point-estimate bounds" if not point_interval_label_used else "_nolegend_",
            zorder=2.7,
        )
        point_interval_label_used = True

    ax.plot(
        xs,
        point_midpoints,
        color=point_color,
        linewidth=1.15,
        linestyle="--",
        alpha=0.38,
        zorder=2.5,
    )

    ax.scatter(
        xs,
        point_midpoints,
        s=60,
        color=point_color,
        edgecolor="white",
        linewidth=1.0,
        label="ATT point estimate",
        zorder=3,
    )

    ax.axhline(
        0,
        linewidth=1.45,
        linestyle="--",
        color=null_color,
        label="Null effect",
        zorder=0,
    )

    if np.isfinite(r2_star):
        ax.axvline(
            r2_star,
            linewidth=1.9,
            linestyle=":",
            color=point_color,
            zorder=0,
        )

    for x in sorted(bmk["plot_x"].dropna().unique()):
        ax.axvline(
            float(x),
            linewidth=1.0,
            color=benchmark_color,
            alpha=0.28,
            linestyle=(0, (1.2, 2.4)),
            zorder=0,
        )

    ymin = min(ci_lower.min(), point_lower.min(), -0.15)
    ymax = max(ci_upper.max(), point_upper.max(), 4.7)

    ax.set_xlim(*xlim)
    ax.set_ylim(ymin - 0.12, ymax + 0.45)
    y_top = ax.get_ylim()[1]

    ax.grid(True, axis="y", color="0.88", linewidth=0.75)
    ax.grid(False, axis="x")
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.9)
    ax.spines["bottom"].set_linewidth(0.9)

    ax.set_xlabel(r"$R^2$", fontsize=11)
    ax.set_ylabel("Estimated ATT", fontsize=11)
    ax.tick_params(axis="both", labelsize=9)

    handles, labels = ax.get_legend_handles_labels()
    legend_order = [
        "95% bootstrap CI",
        "ATT point-estimate bounds",
        "ATT point estimate",
        "Null effect",
    ]
    by_label = dict(zip(labels, handles))

    ax.legend(
        [by_label[label] for label in legend_order if label in by_label],
        [label for label in legend_order if label in by_label],
        loc="upper center",
        bbox_to_anchor=(0.5, -0.14),
        ncol=4,
        frameon=False,
        fontsize=8.5,
        handlelength=2.0,
        columnspacing=1.7,
    )

    # R2* annotation inside the main panel.
    if np.isfinite(r2_star):
        ax.annotate(
            rf"$R^{{2*}}$ = {r2_star:.2f}",
            xy=(r2_star, y_top - 0.85),
            xytext=(r2_star + 0.025, y_top - 0.55),
            ha="left",
            va="top",
            fontsize=9.2,
            bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="0.75", alpha=0.92),
            arrowprops=dict(
                arrowstyle="-",
                lw=0.9,
                alpha=0.7,
                color="0.35",
            ),
        )

    label_box = dict(
        boxstyle="round,pad=0.23",
        fc="white",
        ec=label_color,
        lw=0.8,
        alpha=0.94,
    )

    label_arrow = dict(
        arrowstyle="-",
        lw=0.9,
        alpha=0.62,
        color=label_color,
        shrinkA=2,
        shrinkB=2,
        connectionstyle="angle3,angleA=0,angleB=90",
    )

    def annotate_benchmark_label(
        text: str,
        x_anchor: float,
        x_text: float,
        y_anchor: float,
        y_text: float,
        ha: str,
        fontsize: float,
    ) -> None:
        ax.scatter(
            [x_anchor],
            [y_anchor],
            marker="D",
            s=32,
            color=label_color,
            edgecolor="white",
            linewidth=0.7,
            label="_nolegend_",
            zorder=5,
        )
        ax.annotate(
            text,
            xy=(x_anchor, y_anchor),
            xytext=(x_text, y_text),
            ha=ha,
            va="center",
            fontsize=fontsize,
            color=label_color,
            arrowprops=label_arrow,
            bbox=label_box,
            zorder=6,
        )

    # Group left-cluster covariates into a single central annotation.
    other = bmk.loc[bmk["variable"].isin(other_vars)].copy()
    mid = bmk.loc[bmk["variable"].isin(["income", "race", "education"])].copy()

    # Single label for the left cluster.
    if not other.empty:
        x_anchor = float(other["plot_x"].median())
        x_text = x_anchor + 0.035
        y_anchor = y_top - 0.62
        y_text = y_top - 1.02

        annotate_benchmark_label(
            "Other covariates",
            x_anchor=x_anchor,
            x_text=x_text,
            y_anchor=y_anchor,
            y_text=y_text,
            ha="left",
            fontsize=8.9,
        )

    # Key benchmark labels inside the main panel.
    preferred_positions = {
        "income": {"dx": -0.012, "dy": 1.42, "anchor_dy": 1.18, "ha": "right"},
        "race": {"dx": 0.005, "dy": 1.26, "anchor_dy": 0.94, "ha": "left"},
        "education": {"dx": -0.02, "dy": 0.42, "anchor_dy": 0.76, "ha": "right"},
    }

    for _, row in mid.iterrows():
        var = str(row["variable"])
        x = float(row["plot_x"])
        pos = preferred_positions.get(
            var,
            {"dx": 0.01, "dy": 0.8, "anchor_dy": 0.62, "ha": "left"},
        )
        y_anchor = y_top - pos["anchor_dy"]
        y_text = y_top - pos["dy"]

        annotate_benchmark_label(
            str(row["label"]),
            x_anchor=x,
            x_text=x + pos["dx"],
            y_anchor=y_anchor,
            y_text=y_text,
            ha=pos["ha"],
            fontsize=8.8,
        )

    fig.subplots_adjust(bottom=0.24, top=0.97)

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

    group_spacing = 0.68
    df["x_base"] = df["variable"].apply(
        lambda v: (order.index(v) if v in order else len(order)) * group_spacing
    )

    method_order = ["MSM", "MSM (Qbal)", "VBM", "VBM, w/ Corr."]
    methods = [m for m in method_order if m in df["method_label"].unique()]

    # Keep the four method intervals readable within each covariate group.
    offsets = {
        "MSM": -0.108,
        "MSM (Qbal)": -0.037,
        "VBM": 0.037,
        "VBM, w/ Corr.": 0.108,
    }

    # Muted, paper-like palette with enough contrast across methods.
    method_colors = {
        "MSM": "#F2A13A",
        "MSM (Qbal)": "#A01820",
        "VBM": "#5B9BD5",
        "VBM, w/ Corr.": "#33308C",
    }

    fig, ax = plt.subplots(figsize=(12.9, 5.35))

    tau = pd.to_numeric(df["tau_hat"], errors="coerce").dropna()

    if not tau.empty:
        ax.axhline(
            float(tau.iloc[0]),
            linewidth=0.9,
            linestyle=(0, (1.5, 2.8)),
            color="#2F3437",
            alpha=0.62,
            label="Original ATT",
            zorder=1,
        )

    ax.axhline(
        0,
        linewidth=1.35,
        linestyle="-",
        color="#111111",
        alpha=1.0,
        label="Null effect",
        zorder=2,
    )

    for method in methods:
        sub = df.loc[df["method_label"].eq(method)].sort_values("x_base")

        x = sub["x_base"].to_numpy(dtype=float) + offsets.get(method, 0.0)
        y = sub["tau_hat"].to_numpy(dtype=float)
        boot_lower = pd.to_numeric(sub.get("bootstrap_lower", sub["lower"]), errors="coerce").to_numpy(dtype=float)
        boot_upper = pd.to_numeric(sub.get("bootstrap_upper", sub["upper"]), errors="coerce").to_numpy(dtype=float)
        det_lower = pd.to_numeric(sub.get("deterministic_lower", sub["lower"]), errors="coerce").to_numpy(dtype=float)
        det_upper = pd.to_numeric(sub.get("deterministic_upper", sub["upper"]), errors="coerce").to_numpy(dtype=float)

        boot_lower = np.where(np.isfinite(boot_lower), boot_lower, det_lower)
        boot_upper = np.where(np.isfinite(boot_upper), boot_upper, det_upper)
        det_lower = np.where(np.isfinite(det_lower), det_lower, boot_lower)
        det_upper = np.where(np.isfinite(det_upper), det_upper, boot_upper)

        color = method_colors.get(method, "#4C78A8")

        cap_width = 0.052
        for xi, blo, bhi, dlo, dhi in zip(x, boot_lower, boot_upper, det_lower, det_upper):
            ax.plot(
                [xi, xi],
                [blo, bhi],
                color=color,
                linewidth=2.25,
                alpha=0.38,
                solid_capstyle="butt",
                zorder=2,
            )
            ax.plot(
                [xi - cap_width, xi + cap_width],
                [blo, blo],
                color=color,
                linewidth=1.8,
                alpha=0.38,
                solid_capstyle="butt",
                zorder=2,
            )
            ax.plot(
                [xi - cap_width, xi + cap_width],
                [bhi, bhi],
                color=color,
                linewidth=1.8,
                alpha=0.38,
                solid_capstyle="butt",
                zorder=2,
            )
            ax.plot(
                [xi, xi],
                [dlo, dhi],
                color=color,
                linewidth=7.0,
                alpha=0.94,
                solid_capstyle="butt",
                zorder=3,
            )

        ax.scatter(
            x,
            y,
            s=42,
            color=color,
            edgecolor="white",
            linewidth=0.9,
            label=method,
            zorder=4,
        )

    # Paper-style gridlines: visible enough to structure the panel, but still secondary.
    ax.grid(True, axis="y", color="#C7C9CC", linewidth=1.05)
    ax.grid(True, axis="x", color="#C7C9CC", linewidth=1.05)
    ax.set_axisbelow(True)

    # Reduce side padding so the covariate groups occupy the plot area more evenly.
    ax.set_xlim(-0.5 * group_spacing, (len(order) - 0.5) * group_spacing)

    ax.set_xticks([i * group_spacing for i in range(len(order))])

    ax.set_xticklabels(
        [label_map[v] for v in order],
        rotation=0,
        ha="center",
        fontsize=9.3,
        fontfamily="serif",
    )

    ax.set_yticks([0, 1, 2, 3])
    ax.tick_params(axis="y", labelsize=9.5)
    ax.tick_params(axis="x", length=0)

    ax.set_ylabel("Estimated ATT", fontsize=11.5, fontfamily="serif")
    ax.set_xlabel("")

    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)

    handles = [
        Line2D([0], [0], color="#2F3437", linewidth=0.9, alpha=0.62, linestyle=(0, (1.5, 2.8)), label="Original ATT"),
        Line2D([0], [0], color="#111111", linewidth=1.35, linestyle="-", label="Null effect"),
        Line2D([0], [0], color="#8FAFC8", linewidth=1.65, marker="_", markersize=10, label="95% bootstrap CI"),
        Line2D([0], [0], color="#4C8FC2", linewidth=7.0, solid_capstyle="butt", label="ATT bounds"),
    ]

    handles.extend(
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markerfacecolor=method_colors.get(method, "#4C78A8"),
            markeredgecolor="white",
            markeredgewidth=0.9,
            markersize=6.4,
            label=method,
        )
        for method in methods
    )

    ax.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.14),
        ncol=len(handles),
        frameon=False,
        fontsize=7.8,
        handlelength=1.35,
        columnspacing=0.95,
    )

    fig.subplots_adjust(bottom=0.26, top=0.97, left=0.07, right=0.99)

    save_figures(fig, output_dir, "figure3_benchmark_paper_style")




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

    print(f"Done. PNG figures saved to: {args.output.resolve()}")


if __name__ == "__main__":
    main()
