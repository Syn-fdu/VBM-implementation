#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Enhanced synthetic VBM/MSM plots.

This script is separate from scripts/enhanced_plots.py so the original NHANES
figures are not affected. It reads the current synthetic-extension CSV files
produced by:

    source("extension_01_hidden_strength_vbm_msm.R")

and writes PNG-only Figure-1-style sensitivity plots for every synthetic
dataset, VBM/MSM benchmark interval comparisons, and a cross-dataset detection
power plot along the omitted-confounding strength axis.

Usage:
    python scripts/enhanced_synthetic_plots.py \
        --input output/synthetic_extension/tables \
        --output output/synthetic_extension/figures_enhanced
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


plt.rcParams["font.size"] = 12


LABEL_MAP = {
    "gender": "Gender",
    "age": "Age",
    "income": "Income",
    "income_missing": "Income (Missing)",
    "education": "Education",
    "cig_smoked": "Cig. Smoked",
    "smoking_history": "Smoking history",
    "race": "Race",
    "location_region": "Location (oracle)",
    "hidden_confounder": "Hidden confounder (oracle)",
}

MODEL_ORDER = ["VBM", "MSM"]


def clear_output_dir(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    for p in output_dir.iterdir():
        if p.is_file() and p.suffix.lower() in {".png", ".svg", ".pdf", ".jpg", ".jpeg", ".webp"}:
            p.unlink()


def read_csv(input_dir: Path, name: str) -> pd.DataFrame:
    path = input_dir / name
    if not path.exists():
        raise FileNotFoundError(f"Missing required file: {path}")
    return pd.read_csv(path)


def safe_stem(s: str) -> str:
    out = re.sub(r"[^A-Za-z0-9_]+", "_", str(s)).strip("_").lower()
    return out or "dataset"


def save_figures(fig: plt.Figure, output_dir: Path, stem: str) -> None:
    """Save a figure as PNG only."""
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / f"{stem}.png", format="png", dpi=400, bbox_inches="tight")
    plt.close(fig)

def metric_from_summary(summary: pd.DataFrame, dataset_id: str, column: str) -> float:
    if column not in summary.columns:
        return np.nan
    sub = summary.loc[summary["dataset_id"].astype(str).eq(str(dataset_id)), column]
    if sub.empty:
        return np.nan
    return pd.to_numeric(sub, errors="coerce").iloc[0]


def ordered_summary(summary: pd.DataFrame) -> pd.DataFrame:
    df = summary.copy()
    if "scenario_order" in df.columns:
        df["scenario_order"] = pd.to_numeric(df["scenario_order"], errors="coerce")
        df = df.sort_values(["scenario_order", "dataset_id"], na_position="last")
    else:
        df = df.sort_values("dataset_id")
    return df.reset_index(drop=True)


def as_bool_series(s: pd.Series) -> pd.Series:
    if s.dtype == bool:
        return s
    return s.astype(str).str.lower().isin({"true", "t", "1", "yes", "y"})


def interpolate_curve(df: pd.DataFrame, xs: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    d = df[["R2", "lower", "upper"]].copy()
    d["R2"] = pd.to_numeric(d["R2"], errors="coerce")
    d["lower"] = pd.to_numeric(d["lower"], errors="coerce")
    d["upper"] = pd.to_numeric(d["upper"], errors="coerce")
    d = d.replace([np.inf, -np.inf], np.nan).dropna().sort_values("R2")

    # If coarse and fine rows contain the same R², keep the last one after sorting.
    d = d.groupby("R2", as_index=False).last()

    return (
        np.interp(xs, d["R2"].to_numpy(), d["lower"].to_numpy()),
        np.interp(xs, d["R2"].to_numpy(), d["upper"].to_numpy()),
    )


def assign_label_heights(bmk: pd.DataFrame) -> pd.DataFrame:
    bmk = bmk.sort_values("sensitivity_value").reset_index(drop=True).copy()
    if bmk.empty:
        bmk["x_text"] = []
        bmk["y_text"] = []
        return bmk

    x = bmk["sensitivity_value"].to_numpy(dtype=float)
    bmk["x_text"] = x
    bmk["y_text"] = 0.50

    cluster_threshold = 0.045
    clusters: list[list[int]] = []
    current = [0]

    for i in range(1, len(x)):
        if abs(x[i] - x[i - 1]) < cluster_threshold:
            current.append(i)
        else:
            clusters.append(current)
            current = [i]
    clusters.append(current)

    levels = [0.18, 0.34, 0.50, 0.66, 0.82]

    for cluster in clusters:
        if len(cluster) == 1:
            bmk.loc[cluster[0], "y_text"] = 0.50
        else:
            for j, idx in enumerate(cluster):
                bmk.loc[idx, "y_text"] = levels[j % len(levels)]

    return bmk


def make_synthetic_figure1(
    input_dir: Path,
    output_dir: Path,
    summary: pd.DataFrame,
    dataset_id: str,
    dataset_label: str,
) -> None:
    curve_name = f"{dataset_id}_vbm_bootstrap_curve.csv"
    benchmark_name = f"{dataset_id}_benchmark_vbm_msm.csv"

    curve = read_csv(input_dir, curve_name)
    benchmark = read_csv(input_dir, benchmark_name)

    curve["R2"] = pd.to_numeric(curve["R2"], errors="coerce")
    max_r2_curve = float(np.nanmax(curve["R2"].to_numpy()))

    r2_star = metric_from_summary(summary, dataset_id, "VBM_R2_star_observed_only")
    oracle_bmk_r2 = metric_from_summary(summary, dataset_id, "oracle_omitted_VBM_benchmark_R2")

    bmk = benchmark.loc[
        benchmark["method"].astype(str).eq("VBM benchmark (bootstrap)"),
        ["variable", "sensitivity_value", "benchmark_role"],
    ].copy()

    bmk["sensitivity_value"] = pd.to_numeric(bmk["sensitivity_value"], errors="coerce")
    bmk = bmk.replace([np.inf, -np.inf], np.nan).dropna(subset=["sensitivity_value"])
    bmk = bmk.drop_duplicates(subset=["variable", "sensitivity_value"])

    bmk["label"] = bmk["variable"].astype(str).map(LABEL_MAP).fillna(bmk["variable"].astype(str))

    x_max_candidates = [
        max_r2_curve,
        float(np.nanmax(bmk["sensitivity_value"].to_numpy())) if not bmk.empty else np.nan,
        r2_star,
        oracle_bmk_r2,
        0.70,
    ]

    x_max = np.nanmax(x_max_candidates)
    xlim = (-0.06, min(0.95, x_max + 0.08))

    bmk = bmk.loc[
        (bmk["sensitivity_value"] >= xlim[0]) &
        (bmk["sensitivity_value"] <= xlim[1])
    ].copy()

    bmk = assign_label_heights(bmk)

    # Use the actual grid values, but thin dense fine-grid labels by plotting
    # at a regular 0.1 grid plus the searched R²* and oracle benchmark.
    xs = np.round(np.arange(0.0, min(0.90, xlim[1]) + 1e-9, 0.1), 2)
    extra = [x for x in [r2_star, oracle_bmk_r2] if np.isfinite(x) and xlim[0] <= x <= xlim[1]]
    xs = np.array(sorted(set(np.round(np.concatenate([xs, np.array(extra)]), 3))))

    ci_lower, ci_upper = interpolate_curve(curve, xs)
    midpoints = (ci_lower + ci_upper) / 2

    fig = plt.figure(figsize=(22, 14))
    gs = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[3.2, 9.8], hspace=0.08)

    ax_top = fig.add_subplot(gs[0, 0])
    ax = fig.add_subplot(gs[1, 0], sharex=ax_top)

    fig.suptitle(
        f"Figure 1-style VBM sensitivity curve: {dataset_label}",
        fontsize=24,
        y=0.985,
    )

    ax.vlines(xs, ci_lower, ci_upper, linewidth=3, alpha=0.95, label="95% bootstrap CI", zorder=1)
    ax.scatter(xs, midpoints, s=58, label="Bootstrap interval midpoint", zorder=3)

    ax.axhline(0, linewidth=2.0, linestyle="--", label="Null effect", zorder=0)

    if np.isfinite(r2_star):
        ax.axvline(r2_star, linewidth=2.2, linestyle=":", zorder=0)

    for _, row in bmk.iterrows():
        lw = 1.8 if str(row.get("benchmark_role", "")) == "oracle_omitted_confounder" else 1.1
        alpha = 0.60 if str(row.get("benchmark_role", "")) == "oracle_omitted_confounder" else 0.23
        ax.axvline(float(row["sensitivity_value"]), linewidth=lw, alpha=alpha, zorder=0)

    ymin = min(float(np.nanmin(ci_lower)), -0.15)
    ymax = max(float(np.nanmax(ci_upper)), float(np.nanmax(midpoints)) + 0.5)

    ax.set_xlim(*xlim)
    ax.set_ylim(ymin - 0.12, ymax + 0.45)
    ax.set_xlabel("R²", fontsize=17)
    ax.set_ylabel("Estimated ATT", fontsize=17)
    ax.tick_params(axis="both", labelsize=15)

    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.12),
        ncol=3,
        frameon=False,
        fontsize=14,
        handlelength=2.3,
        columnspacing=2.1,
    )

    ax_top.set_xlim(*xlim)
    ax_top.set_ylim(0, 1)
    ax_top.axis("off")

    ax_top.text(
        xlim[0] + 0.005,
        0.88,
        "Benchmarked covariates",
        ha="left",
        va="top",
        fontsize=15,
        fontweight="bold",
    )

    if np.isfinite(r2_star):
        ax_top.axvline(r2_star, ymin=0.05, ymax=0.78, linewidth=2.2, linestyle=":")
        ax_top.text(
            r2_star + 0.006,
            0.76,
            f"R²* = {r2_star:.2f}",
            rotation=90,
            va="top",
            ha="left",
            fontsize=15,
        )

    for _, row in bmk.iterrows():
        x = float(row["sensitivity_value"])
        y_text = float(row["y_text"])
        label = str(row["label"])

        is_oracle = str(row.get("benchmark_role", "")) == "oracle_omitted_confounder"
        lw = 1.8 if is_oracle else 1.1
        stem_alpha = 0.75 if is_oracle else 0.45

        ax_top.plot([x, x], [0.02, 0.14], linewidth=lw, alpha=stem_alpha)

        ax_top.annotate(
            label,
            xy=(x, 0.14),
            xytext=(x, y_text),
            textcoords="data",
            ha="center",
            va="center",
            fontsize=13,
            fontweight="bold" if is_oracle else "normal",
            bbox=dict(
                boxstyle="round,pad=0.25",
                fc="white",
                ec="0.55" if is_oracle else "0.75",
                alpha=0.97,
            ),
            arrowprops=dict(
                arrowstyle="-",
                lw=1.1 if is_oracle else 0.9,
                alpha=0.70 if is_oracle else 0.55,
                shrinkA=4,
                shrinkB=2,
            ),
            clip_on=False,
        )

    fig.subplots_adjust(top=0.93, bottom=0.13)
    save_figures(fig, output_dir, f"synthetic_figure1_vbm_paper_style_{safe_stem(dataset_id)}")



def make_synthetic_benchmark_comparison(
    input_dir: Path,
    output_dir: Path,
    dataset_id: str,
    dataset_label: str,
) -> None:
    benchmark_name = f"{dataset_id}_benchmark_vbm_msm.csv"
    benchmark = read_csv(input_dir, benchmark_name)

    required_cols = {"variable", "method", "lower", "upper"}
    missing = required_cols.difference(benchmark.columns)
    if missing:
        raise ValueError(f"{benchmark_name} is missing required columns: {sorted(missing)}")

    df = benchmark.copy()
    df["lower"] = pd.to_numeric(df["lower"], errors="coerce")
    df["upper"] = pd.to_numeric(df["upper"], errors="coerce")
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["lower", "upper"])

    method_map = {
        "VBM benchmark (bootstrap)": "VBM",
        "MSM benchmark (bootstrap)": "MSM",
    }
    df = df.loc[df["method"].astype(str).isin(method_map.keys())].copy()
    if df.empty:
        return

    df["method_label"] = df["method"].astype(str).map(method_map)
    df["label"] = df["variable"].astype(str).map(LABEL_MAP).fillna(df["variable"].astype(str))

    ordered_variables = [v for v in LABEL_MAP.keys() if v in set(df["variable"].astype(str))]
    ordered_variables += [
        v for v in df["variable"].astype(str).unique()
        if v not in ordered_variables
    ]

    x_positions = {v: i for i, v in enumerate(ordered_variables)}
    df["x_base"] = df["variable"].astype(str).map(x_positions)

    offsets = {"VBM": -0.13, "MSM": 0.13}

    fig, ax = plt.subplots(figsize=(16, 8.5))

    for method_label in ["VBM", "MSM"]:
        sub = df.loc[df["method_label"].eq(method_label)].copy()
        if sub.empty:
            continue

        sub = sub.sort_values("x_base")
        y_mid = (sub["lower"].to_numpy(dtype=float) + sub["upper"].to_numpy(dtype=float)) / 2.0
        yerr_low = y_mid - sub["lower"].to_numpy(dtype=float)
        yerr_high = sub["upper"].to_numpy(dtype=float) - y_mid
        x = sub["x_base"].to_numpy(dtype=float) + offsets[method_label]

        ax.errorbar(
            x,
            y_mid,
            yerr=np.vstack([yerr_low, yerr_high]),
            fmt="o",
            linewidth=2.2,
            capsize=4,
            markersize=6.5,
            label=method_label,
        )

    if "tau_hat" in df.columns:
        tau_hat = pd.to_numeric(df["tau_hat"], errors="coerce").dropna()
        if not tau_hat.empty:
            ax.axhline(float(tau_hat.iloc[0]), linestyle=":", linewidth=1.8, label="Observed-only ATT")

    ax.axhline(0, linestyle="--", linewidth=1.8, label="Null effect")

    labels = [LABEL_MAP.get(v, v) for v in ordered_variables]
    ax.set_xticks(np.arange(len(ordered_variables)))
    ax.set_xticklabels(labels, rotation=35, ha="right")
    ax.set_ylabel("ATT interval", fontsize=15)
    ax.set_xlabel("Benchmarked covariate", fontsize=15, labelpad=22)
    ax.set_title(
        f"Synthetic benchmark comparison: {dataset_label}",
        fontsize=20,
        pad=22,
    )
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.34),
        ncol=4,
        frameon=False,
        fontsize=12.5,
        handlelength=2.6,
        columnspacing=1.6,
    )
    ax.grid(axis="y", alpha=0.25)
    fig.subplots_adjust(bottom=0.42, top=0.90)

    save_figures(fig, output_dir, f"synthetic_benchmark_vbm_msm_paper_style_{safe_stem(dataset_id)}")



def make_detection_power_plot(input_dir: Path, output_dir: Path, summary: pd.DataFrame) -> None:
    """Compare VBM and MSM null-crossing power across synthetic scenarios.

    The y-axis is a margin relative to the null at the oracle omitted-variable
    benchmark. Positive values mean the sensitivity interval has reached or
    crossed zero; negative values mean the interval remains away from zero.
    """
    path = input_dir / "synthetic_model_detection_power.csv"
    if not path.exists():
        return

    df = pd.read_csv(path)
    required = {
        "dataset_id",
        "dataset_label",
        "model",
        "omitted_confounding_strength",
        "detection_margin",
        "null_reached_or_crossed",
    }
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"synthetic_model_detection_power.csv is missing required columns: {sorted(missing)}")

    df = df.copy()
    df["omitted_confounding_strength"] = pd.to_numeric(df["omitted_confounding_strength"], errors="coerce")
    df["detection_margin"] = pd.to_numeric(df["detection_margin"], errors="coerce")
    if "scenario_order" in df.columns:
        df["scenario_order"] = pd.to_numeric(df["scenario_order"], errors="coerce")
    else:
        order_map = {str(r["dataset_id"]): i for i, r in ordered_summary(summary).iterrows()}
        df["scenario_order"] = df["dataset_id"].astype(str).map(order_map)

    df = df.replace([np.inf, -np.inf], np.nan).dropna(
        subset=["omitted_confounding_strength", "detection_margin", "scenario_order"]
    )
    df = df.loc[df["model"].astype(str).isin(MODEL_ORDER)].copy()
    if df.empty:
        return

    df["null_reached_or_crossed"] = as_bool_series(df["null_reached_or_crossed"])
    df = df.sort_values(["scenario_order", "model"])

    fig, ax = plt.subplots(figsize=(15.5, 8.5))

    for model in MODEL_ORDER:
        sub = df.loc[df["model"].astype(str).eq(model)].sort_values("scenario_order")
        if sub.empty:
            continue
        marker = "o" if model == "VBM" else "s"
        ax.plot(
            sub["omitted_confounding_strength"],
            sub["detection_margin"],
            marker=marker,
            linewidth=2.6,
            markersize=7,
            label=f"{model} detection margin",
        )

        reached = sub.loc[sub["null_reached_or_crossed"]]
        if not reached.empty:
            ax.scatter(
                reached["omitted_confounding_strength"],
                reached["detection_margin"],
                s=115,
                facecolors="none",
                linewidths=2.0,
                label=f"{model} reaches null",
            )

    ax.axhline(0, linestyle="--", linewidth=1.8, label="Detection threshold")

    labels_df = df.drop_duplicates("dataset_id").sort_values("scenario_order")
    tick_x = labels_df["omitted_confounding_strength"].to_numpy(dtype=float)
    tick_labels = [str(v).split(":", 1)[0] for v in labels_df["dataset_label"]]
    ax.set_xticks(tick_x)
    ax.set_xticklabels(tick_labels)

    for _, row in labels_df.iterrows():
        ax.text(
            float(row["omitted_confounding_strength"]),
            -0.12,
            str(row["dataset_label"]),
            rotation=35,
            ha="right",
            va="top",
            fontsize=10,
            transform=ax.get_xaxis_transform(),
            clip_on=False,
        )

    ax.set_xlabel("Omitted confounding strength (designed index)", fontsize=15, labelpad=18)
    ax.set_ylabel("Null-crossing detection margin", fontsize=15)
    ax.set_title("VBM vs MSM detection power as omitted confounding increases", fontsize=20, pad=24)
    ax.text(
        0.01,
        0.97,
        "Positive margin: oracle omitted-confounder benchmark interval reaches/crosses zero",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=12,
        bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="0.8", alpha=0.95),
    )
    ax.grid(axis="y", alpha=0.25)

    # Place the legend outside the plotting area on the right.
    # This keeps it away from the title and avoids overlap with rotated x-axis labels.
    ax.legend(
        loc="center left",
        bbox_to_anchor=(1.02, 0.50),
        ncol=1,
        frameon=False,
        fontsize=12.5,
        handlelength=2.6,
        labelspacing=1.0,
        borderaxespad=0.0,
    )
    fig.subplots_adjust(bottom=0.34, top=0.88, right=0.78)

    save_figures(fig, output_dir, "synthetic_vbm_msm_detection_power_by_strength")


def make_all_figures(input_dir: Path, output_dir: Path) -> None:
    summary = ordered_summary(read_csv(input_dir, "synthetic_dataset_summary.csv"))

    for _, row in summary.iterrows():
        dataset_id = str(row["dataset_id"])
        dataset_label = str(row.get("dataset_label", dataset_id))
        make_synthetic_figure1(input_dir, output_dir, summary, dataset_id, dataset_label)
        make_synthetic_benchmark_comparison(input_dir, output_dir, dataset_id, dataset_label)

    make_detection_power_plot(input_dir, output_dir, summary)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("output/synthetic_extension/tables"),
        help="Directory containing synthetic CSV files.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/synthetic_extension/figures_enhanced"),
        help="Directory to save enhanced synthetic figures.",
    )
    args = parser.parse_args()

    clear_output_dir(args.output)
    make_all_figures(args.input, args.output)

    print(f"Done. Enhanced synthetic PNG figures saved to: {args.output.resolve()}")


if __name__ == "__main__":
    main()
