#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Clean enhanced plots for the omitted-confounder + observed-PS-misspecification extension.

This script is designed for the LONG-format summary table produced by the R extension,
where each row corresponds to one (dataset/scenario, ps_spec) combination and the ATT
estimate is stored in `observed_IPW_ATT`.

Expected input directory:
    output/omitted_ps_misspec_extension/tables/

Expected key files:
    omitted_ps_misspec_summary.csv
    omitted_ps_misspec_vbm_curves.csv      # optional, used for the S3 lower-bound curve

Output directory:
    output/omitted_ps_misspec_extension/figures_enhanced/
"""

from __future__ import annotations

import argparse
import re
import shutil
import sys
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Rectangle


# -----------------------------
# Basic configuration
# -----------------------------

FIG_DPI = 300
FIG_EXT = ".png"

COLORS = {
    "True ATT": "#1f77b4",
    "Oracle full PS": "#ff7f0e",
    "Correct observed PS": "#2ca02c",
    "Misspecified linear PS": "#d62728",
    "Null": "#555555",
    "Detect": "#2ca02c",
    "Miss": "#d62728",
    "Neutral": "#bdbdbd",
}

MARKERS = {
    "True ATT": "o",
    "Oracle full PS": "s",
    "Correct observed PS": "^",
    "Misspecified linear PS": "D",
}

SPEC_ORDER = ["Correct observed PS", "Misspecified linear PS"]


# -----------------------------
# Utilities
# -----------------------------

def normalize_name(x: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(x).lower())


def find_col(df: pd.DataFrame, candidates: Iterable[str], *, required: bool = True, purpose: str = "column") -> Optional[str]:
    """Find a column using exact or normalized matching."""
    cols = list(df.columns)
    norm_map = {normalize_name(c): c for c in cols}

    for cand in candidates:
        if cand in cols:
            return cand
        nc = normalize_name(cand)
        if nc in norm_map:
            return norm_map[nc]

    if required:
        available = "\n  ".join(cols)
        raise KeyError(f"Cannot find {purpose}. Tried {list(candidates)}. Available columns are:\n  {available}")
    return None


def first_existing_file(directory: Path, names: Iterable[str], contains: Optional[str] = None) -> Optional[Path]:
    for name in names:
        p = directory / name
        if p.exists():
            return p
    if contains:
        matches = sorted(directory.glob(f"*{contains}*.csv"))
        if matches:
            return matches[0]
    return None


def project_root_from_script() -> Path:
    # If placed in scripts/, project root is the parent of scripts/.
    here = Path(__file__).resolve()
    if here.parent.name.lower() == "scripts":
        return here.parent.parent
    # Otherwise use current working directory. This makes the script easy to run from project root.
    return Path.cwd().resolve()


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def savefig(out_path: Path) -> None:
    ensure_dir(out_path.parent)
    plt.tight_layout()
    plt.savefig(out_path, dpi=FIG_DPI, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_path}")


def parse_scenario_order(row: pd.Series) -> float:
    for c in ["scenario_order", "_scenario_order"]:
        if c in row.index and pd.notna(row[c]):
            try:
                return float(row[c])
            except Exception:
                pass

    text = " ".join(str(row.get(c, "")) for c in ["dataset_id", "dataset_label"] if c in row.index)
    m = re.search(r"S\s*(\d+)", text, flags=re.I)
    if m:
        return float(m.group(1))
    return 999.0


def short_scenario_label(dataset_id: object, dataset_label: object) -> str:
    text = str(dataset_label) if pd.notna(dataset_label) else str(dataset_id)
    did = str(dataset_id) if pd.notna(dataset_id) else ""
    joined = f"{did} {text}"

    m = re.search(r"S\s*(\d+)", joined, flags=re.I)
    s = f"S{m.group(1)}" if m else text[:10]

    low = joined.lower()
    if "mild" in low:
        return f"{s}\nmild"
    if "strong" in low:
        return f"{s}\nstrong"
    if "severe" in low:
        return f"{s}\nsevere"
    return s


def classify_ps_spec(row: pd.Series) -> str:
    text = " ".join(
        str(row.get(c, ""))
        for c in ["ps_spec", "ps_spec_label", "model", "model_label"]
    ).lower()

    if "miss" in text or "linear" in text:
        return "Misspecified linear PS"
    if "correct" in text:
        return "Correct observed PS"

    # Fallback: keep original label if available.
    if "ps_spec_label" in row.index and pd.notna(row["ps_spec_label"]):
        return str(row["ps_spec_label"])
    if "ps_spec" in row.index and pd.notna(row["ps_spec"]):
        return str(row["ps_spec"])
    return "Unknown PS"


def bool_from_value(x: object) -> bool:
    if pd.isna(x):
        return False
    if isinstance(x, (bool, np.bool_)):
        return bool(x)
    if isinstance(x, (int, float, np.integer, np.floating)):
        return float(x) != 0.0
    return str(x).strip().lower() in {"true", "t", "yes", "y", "1", "detect", "detected"}


def prepare_summary(summary: pd.DataFrame) -> pd.DataFrame:
    required = [
        "dataset_id", "dataset_label", "true_ATT", "observed_IPW_ATT",
        "oracle_full_IPW_ATT", "VBM_R2_star", "hidden_U_benchmark_R2",
        "U_benchmark_says_omitted_can_overturn",
    ]
    for col in required:
        find_col(summary, [col], required=True, purpose=col)

    df = summary.copy()
    df["_scenario_order"] = df.apply(parse_scenario_order, axis=1)
    df["_scenario_short"] = df.apply(
        lambda r: short_scenario_label(r.get("dataset_id", np.nan), r.get("dataset_label", np.nan)), axis=1
    )
    df["_spec_group"] = df.apply(classify_ps_spec, axis=1)

    # Keep only the two PS specifications that matter for this extension.
    keep = df["_spec_group"].isin(SPEC_ORDER)
    if keep.any():
        df = df.loc[keep].copy()

    # Numeric coercion for plotting.
    num_cols = [
        "true_ATT", "unweighted_ATT", "observed_IPW_ATT", "oracle_full_IPW_ATT",
        "VBM_R2_star", "hidden_U_benchmark_R2", "total_error_benchmark_R2",
        "U_interval_lower", "U_interval_upper", "total_interval_lower", "total_interval_upper",
    ]
    for c in num_cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    df["_detect"] = df["U_benchmark_says_omitted_can_overturn"].apply(bool_from_value)
    df["_detection_margin"] = df["hidden_U_benchmark_R2"] - df["VBM_R2_star"]
    df = df.sort_values(["_scenario_order", "_spec_group"]).reset_index(drop=True)
    return df


def scenario_table(summary: pd.DataFrame) -> pd.DataFrame:
    return (
        summary.sort_values("_scenario_order")
        .groupby(["dataset_id", "_scenario_short", "_scenario_order"], as_index=False)
        .agg(
            true_ATT=("true_ATT", "first"),
            oracle_full_IPW_ATT=("oracle_full_IPW_ATT", "first"),
        )
        .sort_values("_scenario_order")
    )


def wide_by_spec(summary: pd.DataFrame, value_col: str) -> pd.DataFrame:
    cols = ["dataset_id", "_scenario_short", "_scenario_order", "_spec_group", value_col]
    sub = summary[cols].dropna(subset=[value_col]).copy()
    wide = (
        sub.pivot_table(
            index=["dataset_id", "_scenario_short", "_scenario_order"],
            columns="_spec_group",
            values=value_col,
            aggfunc="first",
        )
        .reset_index()
        .sort_values("_scenario_order")
    )
    wide.columns.name = None
    return wide


# -----------------------------
# Plot 1: ATT diagnostic
# -----------------------------

def plot_att_diagnostic(summary: pd.DataFrame, out_dir: Path) -> None:
    base = scenario_table(summary)
    att_wide = wide_by_spec(summary, "observed_IPW_ATT")

    labels = list(base["_scenario_short"])
    x = np.arange(len(labels))

    fig, ax = plt.subplots(figsize=(9.2, 5.2))

    # Reference quantities.
    ax.plot(
        x, base["true_ATT"], marker=MARKERS["True ATT"], linewidth=2.2,
        label="True ATT", color=COLORS["True ATT"]
    )
    ax.plot(
        x, base["oracle_full_IPW_ATT"], marker=MARKERS["Oracle full PS"], linewidth=2.2,
        label="Oracle full PS", color=COLORS["Oracle full PS"]
    )

    # Observed-only analyses under two PS specifications.
    for spec in SPEC_ORDER:
        if spec in att_wide.columns:
            y = att_wide[spec].reindex(base.index).to_numpy()
            ax.plot(
                x, y, marker=MARKERS.get(spec, "o"), linewidth=2.2,
                label=spec, color=COLORS.get(spec, None)
            )

    ax.axhline(0, color=COLORS["Null"], linewidth=1.1, linestyle="--", alpha=0.55)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=11)
    ax.set_ylabel("ATT estimate")
    ax.set_title("ATT estimates: omitted U under correct vs misspecified observed PS")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False, ncol=2, loc="best")

    # Small direct annotation of the main contrast in the last scenario if available.
    if len(x) > 0:
        last_i = len(x) - 1
        if "Correct observed PS" in att_wide.columns and "Misspecified linear PS" in att_wide.columns:
            y1 = att_wide.iloc[last_i]["Correct observed PS"]
            y2 = att_wide.iloc[last_i]["Misspecified linear PS"]
            if np.isfinite(y1) and np.isfinite(y2):
                ax.annotate(
                    "PS misspecification\nchanges observed-only ATT",
                    xy=(x[last_i], y2), xytext=(x[last_i] - 0.55, max(y1, y2) + 0.08),
                    arrowprops=dict(arrowstyle="->", lw=1.0),
                    fontsize=9, ha="right"
                )

    savefig(out_dir / f"omitted_ps_misspec_att_diagnostic{FIG_EXT}")


# -----------------------------
# Plot 2: Detection margin
# -----------------------------

def plot_detection_margin(summary: pd.DataFrame, out_dir: Path) -> None:
    wide = wide_by_spec(summary, "_detection_margin")
    labels = list(wide["_scenario_short"])
    x = np.arange(len(labels))
    width = 0.34

    fig, ax = plt.subplots(figsize=(9.0, 5.0))

    offsets = {
        "Correct observed PS": -width / 2,
        "Misspecified linear PS": width / 2,
    }

    for spec in SPEC_ORDER:
        if spec not in wide.columns:
            continue
        y = pd.to_numeric(wide[spec], errors="coerce").to_numpy()
        bars = ax.bar(
            x + offsets[spec], y, width=width,
            label=spec, color=COLORS.get(spec, None), alpha=0.88
        )
        for bar, val in zip(bars, y):
            if not np.isfinite(val):
                continue
            decision = "DETECT" if val >= 0 else "MISS"
            va = "bottom" if val >= 0 else "top"
            y_text = val + (0.015 if val >= 0 else -0.015)
            ax.text(
                bar.get_x() + bar.get_width() / 2, y_text,
                decision, ha="center", va=va, fontsize=8.5, rotation=0
            )

    ax.axhline(0, color="black", linewidth=1.2)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=11)
    ax.set_ylabel("Detection margin: hidden-U R² − R²*")
    ax.set_title("VBM benchmark decision for the omitted confounder")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False, loc="best")

    ax.text(
        0.01, 0.02,
        "Above 0: hidden U is benchmarked as strong enough to overturn.\n"
        "Below 0: hidden U is benchmarked as insufficient.",
        transform=ax.transAxes, fontsize=9, va="bottom",
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="#cccccc", alpha=0.85),
    )

    out = out_dir / f"omitted_ps_misspec_detection_margin{FIG_EXT}"
    savefig(out)

    # Backward-compatible filename used by earlier reports/scripts.
    compat = out_dir / f"omitted_ps_misspec_r2_detection{FIG_EXT}"
    try:
        shutil.copyfile(out, compat)
        print(f"Saved compatibility copy: {compat}")
    except Exception as exc:
        print(f"Warning: could not create compatibility copy {compat}: {exc}")


# -----------------------------
# Plot 3: Minimal conclusion matrix
# -----------------------------

def plot_conclusion_matrix(summary: pd.DataFrame, out_dir: Path) -> None:
    scenario_info = (
        summary[["dataset_id", "_scenario_short", "_scenario_order"]]
        .drop_duplicates()
        .sort_values("_scenario_order")
    )
    scenarios = list(scenario_info["dataset_id"])
    scenario_labels = list(scenario_info["_scenario_short"])

    mat = np.full((len(SPEC_ORDER), len(scenarios)), np.nan)
    for i, spec in enumerate(SPEC_ORDER):
        for j, ds in enumerate(scenarios):
            row = summary[(summary["dataset_id"] == ds) & (summary["_spec_group"] == spec)]
            if not row.empty:
                mat[i, j] = 1.0 if bool(row.iloc[0]["_detect"]) else 0.0

    # 0 = miss, 1 = detect. Neutral handles any missing cells.
    cmap = ListedColormap(["#f2b8b5", "#b7e1b2"])

    fig, ax = plt.subplots(figsize=(8.8, 3.6))
    display = np.where(np.isnan(mat), 0, mat)
    ax.imshow(display, cmap=cmap, vmin=0, vmax=1, aspect="auto")

    ax.set_xticks(np.arange(len(scenario_labels)))
    ax.set_xticklabels(scenario_labels, fontsize=11)
    ax.set_yticks(np.arange(len(SPEC_ORDER)))
    ax.set_yticklabels(SPEC_ORDER, fontsize=11)
    ax.set_title("Can the omitted hidden U overturn the conclusion?")

    for i in range(len(SPEC_ORDER)):
        for j in range(len(scenarios)):
            if np.isnan(mat[i, j]):
                label = "NA"
            else:
                label = "DETECT" if mat[i, j] == 1 else "MISS"
            ax.text(j, i, label, ha="center", va="center", fontsize=12, weight="bold")

    # Highlight scenarios where the two PS specifications lead to different decisions.
    if len(SPEC_ORDER) == 2:
        for j in range(len(scenarios)):
            if np.isfinite(mat[0, j]) and np.isfinite(mat[1, j]) and mat[0, j] != mat[1, j]:
                ax.add_patch(Rectangle((j - 0.5, -0.5), 1, 2, fill=False, lw=2.2, edgecolor="black"))
                ax.text(j, 1.53, "opposite\nconclusion", ha="center", va="top", fontsize=8.5)

    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks(np.arange(-0.5, len(scenario_labels), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(SPEC_ORDER), 1), minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=2)
    ax.tick_params(which="minor", bottom=False, left=False)

    savefig(out_dir / f"omitted_ps_misspec_conclusion_matrix{FIG_EXT}")


# -----------------------------
# Plot 4: S3 lower-bound curve
# -----------------------------

def find_curve_columns(curves: pd.DataFrame) -> tuple[str, str]:
    r2_col = find_col(
        curves,
        ["R2", "R2_grid", "R_squared", "R2_value", "sensitivity_R2", "r2"],
        required=False,
        purpose="R² grid column in curve table",
    )
    if r2_col is None:
        # More permissive fallback.
        for c in curves.columns:
            nc = normalize_name(c)
            if nc in {"r2", "r2grid", "rsquared", "r2value"} or ("r2" in nc and "star" not in nc):
                r2_col = c
                break
    if r2_col is None:
        available = "\n  ".join(curves.columns)
        raise KeyError(f"Cannot find R² grid column in curve table. Available columns are:\n  {available}")

    lower_col = find_col(
        curves,
        ["lower_bound", "lower_CI", "ci_lower", "interval_lower", "bootstrap_lower", "U_lower", "U_interval_lower"],
        required=False,
        purpose="lower-bound column in curve table",
    )
    if lower_col is None:
        candidates = []
        for c in curves.columns:
            nc = normalize_name(c)
            if ("lower" in nc or "lwr" in nc) and "upper" not in nc:
                candidates.append(c)
        if candidates:
            # Prefer a generic/VBM lower bound over total-error lower if both exist.
            non_total = [c for c in candidates if "total" not in normalize_name(c)]
            lower_col = non_total[0] if non_total else candidates[0]
    if lower_col is None:
        available = "\n  ".join(curves.columns)
        raise KeyError(f"Cannot find lower-bound column in curve table. Available columns are:\n  {available}")

    return r2_col, lower_col


def prepare_curves(curves: pd.DataFrame) -> pd.DataFrame:
    df = curves.copy()
    if "dataset_id" not in df.columns and "dataset_label" in df.columns:
        df["dataset_id"] = df["dataset_label"]
    if "dataset_label" not in df.columns and "dataset_id" in df.columns:
        df["dataset_label"] = df["dataset_id"]

    df["_scenario_order"] = df.apply(parse_scenario_order, axis=1)
    df["_scenario_short"] = df.apply(
        lambda r: short_scenario_label(r.get("dataset_id", np.nan), r.get("dataset_label", np.nan)), axis=1
    )
    df["_spec_group"] = df.apply(classify_ps_spec, axis=1)
    return df


def select_severe_dataset(df: pd.DataFrame) -> object:
    text = df[[c for c in ["dataset_id", "dataset_label", "_scenario_short"] if c in df.columns]].astype(str).agg(" ".join, axis=1)
    severe = df[text.str.contains("severe|S3", case=False, regex=True, na=False)]
    if not severe.empty:
        return severe.sort_values("_scenario_order").iloc[-1]["dataset_id"]
    return df.sort_values("_scenario_order").iloc[-1]["dataset_id"]


def interp_y(x_grid: np.ndarray, y_grid: np.ndarray, x: float) -> float:
    mask = np.isfinite(x_grid) & np.isfinite(y_grid)
    xg = x_grid[mask]
    yg = y_grid[mask]
    if len(xg) == 0 or not np.isfinite(x):
        return np.nan
    order = np.argsort(xg)
    xg, yg = xg[order], yg[order]
    if x <= xg.min():
        return float(yg[0])
    if x >= xg.max():
        return float(yg[-1])
    return float(np.interp(x, xg, yg))


def plot_lower_bound_curve(curves: Optional[pd.DataFrame], summary: pd.DataFrame, out_dir: Path) -> None:
    if curves is None or curves.empty:
        print("Skipping lower-bound curve: no curve table was found.")
        return

    curves = prepare_curves(curves)
    r2_col, lower_col = find_curve_columns(curves)
    curves[r2_col] = pd.to_numeric(curves[r2_col], errors="coerce")
    curves[lower_col] = pd.to_numeric(curves[lower_col], errors="coerce")

    severe_id = select_severe_dataset(curves)
    severe = curves[curves["dataset_id"] == severe_id].copy()
    if severe.empty:
        print("Skipping lower-bound curve: cannot identify severe/S3 dataset.")
        return

    # Corresponding summary rows for benchmark markers.
    severe_summary = summary[summary["dataset_id"] == severe_id].copy()
    if severe_summary.empty:
        # Fallback by maximum scenario order.
        max_order = summary["_scenario_order"].max()
        severe_summary = summary[summary["_scenario_order"] == max_order].copy()

    fig, ax = plt.subplots(figsize=(9.2, 5.2))

    all_y = []
    for spec in SPEC_ORDER:
        sub = severe[severe["_spec_group"] == spec].sort_values(r2_col)
        if sub.empty:
            continue
        x = sub[r2_col].to_numpy(dtype=float)
        y = sub[lower_col].to_numpy(dtype=float)
        all_y.extend(y[np.isfinite(y)].tolist())

        ax.plot(
            x, y, linewidth=2.4, marker="o", markersize=4.2,
            label=f"{spec}: lower bound", color=COLORS.get(spec, None)
        )

        sm = severe_summary[severe_summary["_spec_group"] == spec]
        if not sm.empty:
            sm_row = sm.iloc[0]
            r2_star = float(sm_row["VBM_R2_star"]) if pd.notna(sm_row["VBM_R2_star"]) else np.nan
            u_r2 = float(sm_row["hidden_U_benchmark_R2"]) if pd.notna(sm_row["hidden_U_benchmark_R2"]) else np.nan

            if np.isfinite(r2_star):
                ax.axvline(
                    r2_star, color=COLORS.get(spec, None), linestyle="--", linewidth=1.5, alpha=0.75
                )
                ax.text(
                    r2_star, 0.04, "R²*", rotation=90, ha="right", va="bottom",
                    color=COLORS.get(spec, None), fontsize=8.5
                )

            if np.isfinite(u_r2):
                y_u = interp_y(x, y, u_r2)
                ax.scatter(
                    [u_r2], [y_u], s=85, color=COLORS.get(spec, None),
                    edgecolor="black", linewidth=0.7, zorder=5
                )
                decision = "detect" if y_u <= 0 else "miss"
                ax.annotate(
                    f"hidden U\n{decision}",
                    xy=(u_r2, y_u), xytext=(u_r2 + 0.02, y_u + 0.55),
                    arrowprops=dict(arrowstyle="->", lw=0.9, color=COLORS.get(spec, None)),
                    color=COLORS.get(spec, None), fontsize=8.5, ha="left"
                )

    ax.axhline(0, color="black", linestyle="--", linewidth=1.2, label="Null effect")
    ax.set_xlabel("R²")
    ax.set_ylabel("Lower bound of VBM interval")
    severe_label = severe["_scenario_short"].iloc[0] if "_scenario_short" in severe.columns else "S3"
    ax.set_title(f"Lower-bound curve in {severe_label}: false-negative mechanism")
    ax.grid(axis="both", alpha=0.22)
    ax.legend(frameon=False, loc="best", fontsize=9)

    # Clip extreme negative tails so the crossing region remains readable.
    all_y = np.asarray(all_y, dtype=float)
    all_y = all_y[np.isfinite(all_y)]
    if len(all_y) > 0:
        ymin = max(float(np.nanpercentile(all_y, 2)) - 0.3, -8.0)
        ymax = max(float(np.nanpercentile(all_y, 98)) + 0.4, 1.25)
        if ymin < ymax:
            ax.set_ylim(ymin, ymax)

    savefig(out_dir / f"omitted_ps_misspec_lower_bound_curve_severe{FIG_EXT}")


# -----------------------------
# Main
# -----------------------------

def load_inputs(tables_dir: Path) -> tuple[pd.DataFrame, Optional[pd.DataFrame]]:
    summary_path = first_existing_file(
        tables_dir,
        [
            "omitted_ps_misspec_summary.csv",
            "ps_misspec_omitted_summary.csv",
            "summary.csv",
        ],
        contains="summary",
    )
    if summary_path is None:
        raise FileNotFoundError(f"Cannot find summary CSV under {tables_dir}")

    curves_path = first_existing_file(
        tables_dir,
        [
            "omitted_ps_misspec_vbm_curves.csv",
            "ps_misspec_omitted_vbm_curves.csv",
            "vbm_curves.csv",
        ],
        contains="curves",
    )

    print(f"Reading summary table: {summary_path}")
    summary = pd.read_csv(summary_path)

    curves = None
    if curves_path is not None:
        print(f"Reading curve table:   {curves_path}")
        curves = pd.read_csv(curves_path)
    else:
        print("Curve table not found; the S3 lower-bound curve will be skipped.")

    return summary, curves


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate clean enhanced plots for omitted-U + PS-misspecification VBM extension.")
    parser.add_argument("--project-root", type=str, default=None, help="Project root. Default: parent of scripts/ or current directory.")
    parser.add_argument("--tables-dir", type=str, default=None, help="Input tables directory. Default: output/omitted_ps_misspec_extension/tables")
    parser.add_argument("--figures-dir", type=str, default=None, help="Output figures directory. Default: output/omitted_ps_misspec_extension/figures_enhanced")
    args = parser.parse_args()

    project_root = Path(args.project_root).resolve() if args.project_root else project_root_from_script()
    tables_dir = Path(args.tables_dir).resolve() if args.tables_dir else project_root / "output" / "omitted_ps_misspec_extension" / "tables"
    figures_dir = Path(args.figures_dir).resolve() if args.figures_dir else project_root / "output" / "omitted_ps_misspec_extension" / "figures_enhanced"

    print(f"Project root: {project_root}")
    print(f"Reading tables from: {tables_dir}")
    print(f"Writing figures to:  {figures_dir}\n")

    try:
        raw_summary, curves = load_inputs(tables_dir)
        summary = prepare_summary(raw_summary)
        ensure_dir(figures_dir)

        plot_att_diagnostic(summary, figures_dir)
        plot_detection_margin(summary, figures_dir)
        plot_conclusion_matrix(summary, figures_dir)
        plot_lower_bound_curve(curves, summary, figures_dir)

        print("\nDone. Enhanced plots generated successfully.")
    except Exception as exc:
        print("\nERROR while generating enhanced plots:", file=sys.stderr)
        print(str(exc), file=sys.stderr)
        print(
            "\nThis script expects a LONG-format summary table with columns such as:\n"
            "  dataset_id, dataset_label, ps_spec, ps_spec_label, true_ATT,\n"
            "  observed_IPW_ATT, oracle_full_IPW_ATT, VBM_R2_star,\n"
            "  hidden_U_benchmark_R2, U_benchmark_says_omitted_can_overturn\n",
            file=sys.stderr,
        )
        raise


if __name__ == "__main__":
    main()
