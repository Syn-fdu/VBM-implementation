#!/usr/bin/env python3
"""
Enhanced plots for the current good-overlap instability extension.

Run from the project root:

    python scripts/enhanced_good_overlap_plots.py

Inputs:
    output/good_overlap_instability_extension/tables/good_overlap_instability_metrics.csv

Outputs:
    output/good_overlap_instability_extension/enhanced_figures/
        enhanced_denominator_instability.png
        enhanced_ratio_and_flag_summary.png

The script is intentionally read-only with respect to the main workflow:
it only reads the extension CSV and writes enhanced figures to a separate folder.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def _find_project_root() -> Path:
    here = Path(__file__).resolve()
    return here.parents[1]


def _as_float_array(series: pd.Series) -> np.ndarray:
    return pd.to_numeric(series, errors="coerce").to_numpy(dtype=float)


def _safe_display_ratio(values: np.ndarray, cap: float = 1e6) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    return np.where(np.isfinite(values), np.minimum(values, cap), cap)


def _add_instability_band(ax, x, flag_rate, threshold=0.5):
    flagged = np.asarray(flag_rate, dtype=float) >= threshold
    if not flagged.any():
        return
    # Draw vertical boundary markers instead of colored shading, so the plot
    # remains clean in both grayscale and color.
    flagged_x = np.asarray(x)[flagged]
    ax.axvline(flagged_x.min(), linestyle=":", linewidth=1)
    ax.text(
        flagged_x.min(),
        ax.get_ylim()[1],
        "unstable region",
        rotation=90,
        va="top",
        ha="right",
        fontsize=9,
    )


def plot_denominator_instability(metrics: pd.DataFrame, out_path: Path) -> None:
    x = _as_float_array(metrics["observed_strength"])
    var_obs = _as_float_array(metrics["var_observed_control_median"])
    var_ideal = _as_float_array(metrics["var_ideal_control_median"])
    epsilon = _as_float_array(metrics["epsilon_added_to_denominator_median"])
    flag_rate = _as_float_array(metrics["flag_rate"])

    order = np.argsort(x)
    x, var_obs, var_ideal, epsilon, flag_rate = (
        x[order], var_obs[order], var_ideal[order], epsilon[order], flag_rate[order]
    )

    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    ax.plot(x, var_obs, marker="o", linewidth=2, label="Observed weight variance")
    ax.plot(x, var_ideal, marker="s", linewidth=2, label="Ideal weight variance")
    ax.plot(x, epsilon, linestyle="--", linewidth=2, label="Stabilizing epsilon")

    ax.set_yscale("log")
    ax.set_xscale("symlog", linthresh=0.001)
    ax.set_xlabel("Observed propensity strength")
    ax.set_ylabel("Control-weight variance (log scale)")
    ax.set_title("Enhanced plot 1: near-zero denominator under good overlap")
    ax.grid(True, which="both", linewidth=0.6, alpha=0.35)

    ax.annotate(
        "Observed weights are almost constant;\nVBM denominator is near zero.",
        xy=(x[0], max(var_obs[0], np.nanmin(var_obs[var_obs > 0]))),
        xytext=(0.006, np.nanmax(var_ideal) / 2),
        arrowprops={"arrowstyle": "->"},
        fontsize=9,
    )

    ax.annotate(
        "As observed assignment strengthens,\nvariance recovers and the diagnostic stabilizes.",
        xy=(0.05, var_obs[np.where(x == 0.05)[0][0]] if np.any(x == 0.05) else np.nanmedian(var_obs)),
        xytext=(0.07, np.nanmin(var_ideal) * 2.0),
        arrowprops={"arrowstyle": "->"},
        fontsize=9,
    )

    _add_instability_band(ax, x, flag_rate)
    ax.legend(frameon=False, loc="best")
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def plot_ratio_and_flag_summary(metrics: pd.DataFrame, out_path: Path) -> None:
    x = _as_float_array(metrics["observed_strength"])
    original_ratio = _safe_display_ratio(_as_float_array(metrics["original_variance_ratio_median"]))
    stabilized_ratio = _safe_display_ratio(_as_float_array(metrics["stabilized_variance_ratio_median"]))
    flag_rate = _as_float_array(metrics["flag_rate"])

    order = np.argsort(x)
    x, original_ratio, stabilized_ratio, flag_rate = (
        x[order], original_ratio[order], stabilized_ratio[order], flag_rate[order]
    )

    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    ax.plot(x, original_ratio, marker="o", linewidth=2, label="Original variance ratio")
    ax.plot(x, stabilized_ratio, marker="s", linestyle="--", linewidth=2, label="Stabilized variance ratio")
    ax.axhline(100, linestyle=":", linewidth=1.5, label="Warning threshold: ratio = 100")

    ax.set_yscale("log")
    ax.set_xscale("symlog", linthresh=0.001)
    ax.set_xlabel("Observed propensity strength")
    ax.set_ylabel("Variance ratio (log scale)")
    ax.set_title("Enhanced plot 2: instability flag aligns with inflated variance ratio")
    ax.grid(True, which="both", linewidth=0.6, alpha=0.35)

    ax2 = ax.twinx()
    ax2.plot(x, flag_rate, marker="^", linestyle="-.", linewidth=1.8, label="Flag rate")
    ax2.set_ylabel("Instability flag rate")
    ax2.set_ylim(-0.02, 1.05)

    ax.annotate(
        "Near perfect overlap:\nratio explodes because denominator is tiny.",
        xy=(x[0], original_ratio[0]),
        xytext=(0.003, original_ratio[0] / 20),
        arrowprops={"arrowstyle": "->"},
        fontsize=9,
    )

    ax.annotate(
        "Stabilization prevents numerical explosion\nbut still preserves the warning signal.",
        xy=(0.01, stabilized_ratio[np.where(x == 0.01)[0][0]] if np.any(x == 0.01) else np.nanmedian(stabilized_ratio)),
        xytext=(0.03, 300),
        arrowprops={"arrowstyle": "->"},
        fontsize=9,
    )

    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, frameon=False, loc="best")

    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def main() -> None:
    project_root = _find_project_root()
    input_path = project_root / "output" / "good_overlap_instability_extension" / "tables" / "good_overlap_instability_metrics.csv"
    output_dir = project_root / "output" / "good_overlap_instability_extension" / "enhanced_figures"
    output_dir.mkdir(parents=True, exist_ok=True)

    metrics = pd.read_csv(input_path)
    plot_denominator_instability(metrics, output_dir / "enhanced_denominator_instability.png")
    plot_ratio_and_flag_summary(metrics, output_dir / "enhanced_ratio_and_flag_summary.png")

    print("Enhanced plots written to:", output_dir)


if __name__ == "__main__":
    main()
