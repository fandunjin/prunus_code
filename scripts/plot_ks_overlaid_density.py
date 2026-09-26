#!/usr/bin/env python3
"""Plot overlaid, boundary-corrected Ks density curves for Figure 2.

The curves are probability densities, so comparisons with different numbers of
syntenic anchors remain comparable. The script intentionally refuses the old
one-pair-per-block result unless --allow-legacy-input is supplied.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator
import numpy as np


plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["font.size"] = 7
plt.rcParams["axes.linewidth"] = 0.8
plt.rcParams["axes.spines.right"] = False
plt.rcParams["axes.spines.top"] = False
plt.rcParams["legend.frameon"] = False
plt.rcParams["xtick.major.width"] = 0.8
plt.rcParams["ytick.major.width"] = 0.8


LABELS = {
    "pzi_hapA1_vs_pzi_hapB1": "P. zippeliana A1-B1",
    "pzi_hapA2_vs_pzi_hapB2": "P. zippeliana A2-B2",
    "psp_hapA1_vs_psp_hapB1": "P. spinulosa A1-B1",
    "psp_hapA2_vs_psp_hapB2": "P. spinulosa A2-B2",
    "ppa_hapA_vs_ppa_hapB": "P. padus A-B",
    "pzi_hapA1_vs_psp_hapA1": "P. zippeliana A1-P. spinulosa A1",
    "pzi_hapA2_vs_psp_hapA2": "P. zippeliana A2-P. spinulosa A2",
    "pzi_hapB1_vs_psp_hapB1": "P. zippeliana B1-P. spinulosa B1",
    "pzi_hapB2_vs_psp_hapB2": "P. zippeliana B2-P. spinulosa B2",
    "pzi_hapA1_vs_ppa_hapA": "P. zippeliana A1-P. padus A",
    "pzi_hapB1_vs_ppa_hapB": "P. zippeliana B1-P. padus B",
    "pzi_hapA1_vs_pcam": "P. zippeliana A1-P. campanulata",
    "pzi_hapB1_vs_pcam": "P. zippeliana B1-P. campanulata",
    "pzi_hapA1_vs_ppe": "P. zippeliana A1-P. persica",
    "pzi_hapB1_vs_ppe": "P. zippeliana B1-P. persica",
    "pzi_hapA1_vs_gtr": "P. zippeliana A1-G. trifoliata",
    "pzi_hapB1_vs_gtr": "P. zippeliana B1-G. trifoliata",
}


PRESETS = {
    "focal": [
        "pzi_hapA1_vs_pzi_hapB1",
        "psp_hapA1_vs_psp_hapB1",
        "ppa_hapA_vs_ppa_hapB",
    ],
    "sensitivity": [
        "pzi_hapA1_vs_pzi_hapB1",
        "pzi_hapA2_vs_pzi_hapB2",
        "psp_hapA1_vs_psp_hapB1",
        "psp_hapA2_vs_psp_hapB2",
        "ppa_hapA_vs_ppa_hapB",
    ],
    "cross_species": [
        "pzi_hapA1_vs_psp_hapA1",
        "pzi_hapA2_vs_psp_hapA2",
        "pzi_hapB1_vs_psp_hapB1",
        "pzi_hapB2_vs_psp_hapB2",
        "pzi_hapA1_vs_ppa_hapA",
        "pzi_hapB1_vs_ppa_hapB",
    ],
    "diploid_references": [
        "pzi_hapA1_vs_pcam",
        "pzi_hapB1_vs_pcam",
        "pzi_hapA1_vs_ppe",
        "pzi_hapB1_vs_ppe",
        "pzi_hapA1_vs_gtr",
        "pzi_hapB1_vs_gtr",
    ],
    "all": list(LABELS),
}


PRESET_DISPLAY_MAX = {
    "focal": 0.50,
    "sensitivity": 0.50,
    "cross_species": 0.35,
    "diploid_references": 1.00,
    "all": 3.00,
}


PRESET_ROLES = {
    "focal": "main_figure_2_ks_panel",
    "sensitivity": "supplementary_haplotype_robustness",
    "cross_species": "supplementary_corresponding_subgenomes",
    "diploid_references": "supplementary_diploid_outgroup_context",
    "all": "quality_control_only",
}


COLORS = {
    "pzi_hapA1_vs_pzi_hapB1": "#B64342",
    "pzi_hapA2_vs_pzi_hapB2": "#E39A94",
    "psp_hapA1_vs_psp_hapB1": "#0F4D92",
    "psp_hapA2_vs_psp_hapB2": "#7FA6D6",
    "ppa_hapA_vs_ppa_hapB": "#42949E",
    "pzi_hapA1_vs_psp_hapA1": "#6F4C9B",
    "pzi_hapA2_vs_psp_hapA2": "#B89BD9",
    "pzi_hapB1_vs_psp_hapB1": "#2E7D32",
    "pzi_hapB2_vs_psp_hapB2": "#8BCF8B",
    "pzi_hapA1_vs_ppa_hapA": "#C98223",
    "pzi_hapB1_vs_ppa_hapB": "#6B6B6B",
    "pzi_hapA1_vs_pcam": "#C98223",
    "pzi_hapB1_vs_pcam": "#C98223",
    "pzi_hapA1_vs_ppe": "#6F4C9B",
    "pzi_hapB1_vs_ppe": "#6F4C9B",
    "pzi_hapA1_vs_gtr": "#4D4D4D",
    "pzi_hapB1_vs_gtr": "#4D4D4D",
}


LINESTYLES = {
    "pzi_hapA2_vs_pzi_hapB2": "--",
    "psp_hapA2_vs_psp_hapB2": "--",
    "pzi_hapA2_vs_psp_hapA2": "--",
    "pzi_hapB2_vs_psp_hapB2": "--",
    "pzi_hapB1_vs_pcam": "--",
    "pzi_hapB1_vs_ppe": "--",
    "pzi_hapB1_vs_gtr": "--",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Merged 4DTV/KaKs TSV")
    parser.add_argument("--outdir", required=True)
    parser.add_argument(
        "--presets",
        nargs="+",
        choices=sorted(PRESETS),
        default=["focal", "sensitivity", "cross_species", "diploid_references"],
    )
    parser.add_argument("--min-ks", type=float, default=0.0)
    parser.add_argument("--max-ks", type=float, default=3.0)
    parser.add_argument(
        "--bandwidth",
        type=float,
        default=0.02,
        help="Common absolute Gaussian bandwidth on the Ks scale",
    )
    parser.add_argument(
        "--bandwidth-sensitivity",
        nargs="*",
        type=float,
        default=[0.02, 0.04, 0.06],
        help="Bandwidths shown in the focal three-panel sensitivity figure; use no values to disable",
    )
    parser.add_argument("--grid-points", type=int, default=601)
    parser.add_argument("--width-mm", type=float, default=183.0)
    parser.add_argument("--height-mm", type=float, default=105.0)
    parser.add_argument("--allow-legacy-input", action="store_true")
    return parser.parse_args()


def parse_float(value: str) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def read_ks(path: Path, min_ks: float, max_ks: float):
    grouped: dict[str, list[float]] = {}
    dropped: dict[str, int] = {}
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            comparison = row.get("comparison", "")
            if not comparison:
                continue
            value = parse_float(row.get("Ks", ""))
            if row.get("kaks_status", "ok") != "ok" or value is None or not min_ks < value < max_ks:
                dropped[comparison] = dropped.get(comparison, 0) + 1
                continue
            grouped.setdefault(comparison, []).append(value)
    return {key: np.asarray(values, dtype=float) for key, values in grouped.items()}, dropped


def reflected_gaussian_kde(values: np.ndarray, grid: np.ndarray, bandwidth: float) -> np.ndarray:
    """Gaussian KDE with reflection at the biological boundary Ks=0."""
    density = np.zeros_like(grid, dtype=float)
    denominator = values.size * bandwidth * math.sqrt(2.0 * math.pi)
    for start in range(0, values.size, 4096):
        chunk = values[start : start + 4096, None]
        z_direct = (grid[None, :] - chunk) / bandwidth
        z_reflect = (grid[None, :] + chunk) / bandwidth
        density += np.exp(-0.5 * z_direct * z_direct).sum(axis=0)
        density += np.exp(-0.5 * z_reflect * z_reflect).sum(axis=0)
    return density / denominator


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def plot_preset(
    preset: str,
    comparisons: list[str],
    values_by_comparison: dict[str, np.ndarray],
    display_min: float,
    display_max: float,
    grid_points: int,
    bandwidth: float,
    outdir: Path,
    width_mm: float,
    height_mm: float,
    filter_min: float,
    filter_max: float,
    output_stem: str | None = None,
) -> None:
    available = [comparison for comparison in comparisons if comparison in values_by_comparison]
    if not available:
        raise SystemExit(f"No valid Ks values were found for preset: {preset}")

    grid = np.linspace(display_min, display_max, grid_points)
    densities = {
        comparison: reflected_gaussian_kde(values_by_comparison[comparison], grid, bandwidth)
        for comparison in available
    }
    fig, ax = plt.subplots(figsize=(width_mm / 25.4, height_mm / 25.4))
    fill_alpha = 0.13 if len(available) <= 5 else 0.045
    line_width = 1.5 if len(available) <= 5 else 1.1

    summary_rows = []
    source_rows = []
    for index, x_value in enumerate(grid):
        row: dict[str, object] = {"Ks": f"{x_value:.6f}"}
        for comparison in available:
            row[comparison] = f"{densities[comparison][index]:.10f}"
        source_rows.append(row)

    for comparison in available:
        values = values_by_comparison[comparison]
        density = densities[comparison]
        color = COLORS.get(comparison, "#555555")
        linestyle = LINESTYLES.get(comparison, "-")
        peak_index = int(np.argmax(density))
        peak = float(grid[peak_index])
        median = float(np.median(values))
        label = f"{LABELS.get(comparison, comparison)}  n={values.size:,}"
        ax.fill_between(grid, density, color=color, alpha=fill_alpha, linewidth=0)
        ax.plot(grid, density, color=color, linestyle=linestyle, linewidth=line_width, label=label)
        median_height = float(np.interp(median, grid, density))
        ax.vlines(median, 0, median_height, color=color, linestyle=":", linewidth=0.75, alpha=0.8)
        summary_rows.append(
            {
                "preset": preset,
                "comparison": comparison,
                "display_label": LABELS.get(comparison, comparison),
                "n_valid_ks": values.size,
                "ks_median": f"{median:.8f}",
                "kde_peak": f"{peak:.8f}",
                "bandwidth": f"{bandwidth:.8f}",
                "ks_filter": f"{filter_min:g}<Ks<{filter_max:g}",
                "display_range": f"{display_min:g}-{display_max:g}",
            }
        )

    ax.set_xlim(float(grid[0]), float(grid[-1]))
    ax.set_ylim(bottom=0)
    if display_max <= 0.5:
        ax.xaxis.set_major_locator(MultipleLocator(0.1))
    elif display_max <= 1.2:
        ax.xaxis.set_major_locator(MultipleLocator(0.2))
    ax.set_xlabel(r"Synonymous substitutions per synonymous site ($K_{\mathrm{s}}$)")
    ax.set_ylabel("Probability density")
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.45, alpha=0.65)
    ax.set_axisbelow(True)
    legend_columns = 3 if len(available) > 12 else (2 if len(available) > 5 else 1)
    handles, labels = ax.get_legend_handles_labels()
    handles.append(Line2D([0], [0], color="#555555", linestyle=":", linewidth=0.9))
    labels.append("Median")
    ax.legend(
        handles,
        labels,
        loc="upper right",
        ncol=legend_columns,
        fontsize=6.2,
        handlelength=2.6,
        columnspacing=1.2,
        labelspacing=0.55,
    )
    fig.tight_layout(pad=0.8)

    base = outdir / (output_stem or f"ks_density_{preset}")
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(
        base.with_suffix(".tiff"),
        dpi=600,
        bbox_inches="tight",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(fig)

    write_tsv(
        base.with_suffix(".summary.tsv"),
        summary_rows,
        [
            "preset",
            "comparison",
            "display_label",
            "n_valid_ks",
            "ks_median",
            "kde_peak",
            "bandwidth",
            "ks_filter",
            "display_range",
        ],
    )
    write_tsv(
        base.with_suffix(".source_data.tsv"),
        source_rows,
        ["Ks"] + available,
    )


def plot_focal_with_references(
    values_by_comparison: dict[str, np.ndarray],
    bandwidth: float,
    outdir: Path,
    filter_min: float,
    filter_max: float,
    grid_points: int,
) -> None:
    """Place homeolog and diploid-reference distributions in separate panels."""
    panel_specs = [
        ("a", PRESETS["focal"], 0.50),
        ("b", PRESETS["diploid_references"], 1.00),
    ]
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(183.0 / 25.4, 78.0 / 25.4),
        sharey=False,
    )
    summary_rows = []
    source_rows = []

    for ax, (panel, comparisons, display_max) in zip(axes, panel_specs):
        available = [comparison for comparison in comparisons if comparison in values_by_comparison]
        if not available:
            raise SystemExit(f"No valid Ks values were found for combined panel {panel}")
        grid = np.linspace(0.0, display_max, grid_points)
        fill_alpha = 0.10 if panel == "a" else 0.045
        line_width = 1.35 if panel == "a" else 1.05

        for comparison in available:
            values = values_by_comparison[comparison]
            density = reflected_gaussian_kde(values, grid, bandwidth)
            color = COLORS.get(comparison, "#555555")
            linestyle = LINESTYLES.get(comparison, "-")
            median = float(np.median(values))
            peak = float(grid[int(np.argmax(density))])
            median_height = float(np.interp(median, grid, density))
            ax.fill_between(grid, density, color=color, alpha=fill_alpha, linewidth=0)
            ax.plot(
                grid,
                density,
                color=color,
                linestyle=linestyle,
                linewidth=line_width,
                label=f"{LABELS.get(comparison, comparison)}  n={values.size:,}",
            )
            ax.vlines(median, 0, median_height, color=color, linestyle=":", linewidth=0.70, alpha=0.8)
            summary_rows.append(
                {
                    "panel": panel,
                    "comparison_class": "within_species_homeolog" if panel == "a" else "interspecies_reference",
                    "comparison": comparison,
                    "display_label": LABELS.get(comparison, comparison),
                    "n_valid_ks": values.size,
                    "ks_median": f"{median:.8f}",
                    "kde_peak": f"{peak:.8f}",
                    "bandwidth": f"{bandwidth:.8f}",
                    "ks_filter": f"{filter_min:g}<Ks<{filter_max:g}",
                    "display_range": f"0-{display_max:g}",
                }
            )
            for index, x_value in enumerate(grid):
                source_rows.append(
                    {
                        "panel": panel,
                        "comparison": comparison,
                        "Ks": f"{x_value:.6f}",
                        "density": f"{density[index]:.10f}",
                    }
                )

        ax.set_xlim(0.0, display_max)
        ax.set_ylim(bottom=0)
        ax.xaxis.set_major_locator(MultipleLocator(0.1 if display_max <= 0.5 else 0.2))
        ax.set_xlabel(r"Synonymous substitutions per synonymous site ($K_{\mathrm{s}}$)")
        ax.grid(axis="y", color="#D9D9D9", linewidth=0.4, alpha=0.6)
        ax.set_axisbelow(True)
        ax.text(
            -0.12,
            1.02,
            panel,
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=8,
            fontweight="bold",
        )
        handles, labels = ax.get_legend_handles_labels()
        handles.append(Line2D([0], [0], color="#555555", linestyle=":", linewidth=0.9))
        labels.append("Median")
        ax.legend(
            handles,
            labels,
            loc="upper right",
            ncol=1,
            fontsize=5.7,
            handlelength=2.5,
            columnspacing=0.8,
            labelspacing=0.45,
        )
    axes[0].set_ylabel("Probability density")
    axes[1].set_ylabel("Probability density")
    fig.tight_layout(pad=0.7, w_pad=1.7)

    base = outdir / "ks_density_focal_with_references"
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(base.with_suffix(".tiff"), dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    write_tsv(
        base.with_suffix(".summary.tsv"),
        summary_rows,
        [
            "panel",
            "comparison_class",
            "comparison",
            "display_label",
            "n_valid_ks",
            "ks_median",
            "kde_peak",
            "bandwidth",
            "ks_filter",
            "display_range",
        ],
    )
    write_tsv(
        base.with_suffix(".source_data.tsv"),
        source_rows,
        ["panel", "comparison", "Ks", "density"],
    )


def plot_bandwidth_sensitivity(
    values_by_comparison: dict[str, np.ndarray],
    bandwidths: list[float],
    outdir: Path,
    filter_min: float,
    filter_max: float,
    grid_points: int,
) -> None:
    comparisons = [comparison for comparison in PRESETS["focal"] if comparison in values_by_comparison]
    if not comparisons or not bandwidths:
        return
    if any(value <= 0 for value in bandwidths):
        raise SystemExit("All bandwidth-sensitivity values must be positive")

    display_min = 0.0
    display_max = PRESET_DISPLAY_MAX["focal"]
    grid = np.linspace(display_min, display_max, grid_points)
    fig, axes = plt.subplots(
        1,
        len(bandwidths),
        figsize=(183.0 / 25.4, 62.0 / 25.4),
        sharex=True,
        sharey=False,
        squeeze=False,
    )
    axes = axes[0]
    source_rows = []
    summary_rows = []
    for panel_index, (ax, bandwidth) in enumerate(zip(axes, bandwidths)):
        for comparison in comparisons:
            values = values_by_comparison[comparison]
            density = reflected_gaussian_kde(values, grid, bandwidth)
            color = COLORS[comparison]
            linestyle = LINESTYLES.get(comparison, "-")
            label = f"{LABELS[comparison]}  n={values.size:,}"
            ax.fill_between(grid, density, color=color, alpha=0.11, linewidth=0)
            ax.plot(grid, density, color=color, linestyle=linestyle, linewidth=1.25, label=label)
            median = float(np.median(values))
            median_height = float(np.interp(median, grid, density))
            ax.vlines(median, 0, median_height, color=color, linestyle=":", linewidth=0.65, alpha=0.8)
            summary_rows.append(
                {
                    "panel": chr(ord("a") + panel_index),
                    "comparison": comparison,
                    "n_valid_ks": values.size,
                    "ks_median": f"{median:.8f}",
                    "kde_peak": f"{float(grid[int(np.argmax(density))]):.8f}",
                    "bandwidth": f"{bandwidth:.8f}",
                    "ks_filter": f"{filter_min:g}<Ks<{filter_max:g}",
                    "display_range": f"{display_min:g}-{display_max:g}",
                }
            )
            for index, x_value in enumerate(grid):
                source_rows.append(
                    {
                        "panel": chr(ord("a") + panel_index),
                        "bandwidth": f"{bandwidth:.8f}",
                        "comparison": comparison,
                        "Ks": f"{x_value:.6f}",
                        "density": f"{density[index]:.10f}",
                    }
                )
        ax.set_xlim(display_min, display_max)
        ax.set_ylim(bottom=0)
        ax.xaxis.set_major_locator(MultipleLocator(0.1))
        ax.grid(axis="y", color="#D9D9D9", linewidth=0.4, alpha=0.6)
        ax.set_axisbelow(True)
        ax.text(0.02, 0.98, chr(ord("a") + panel_index), transform=ax.transAxes, ha="left", va="top", fontsize=8, fontweight="bold")
        ax.text(0.98, 0.98, f"h = {bandwidth:g}", transform=ax.transAxes, ha="right", va="top", fontsize=6.5)
        ax.set_xlabel(r"$K_{\mathrm{s}}$")
        if panel_index == 0:
            ax.set_ylabel("Probability density")
    handles, labels = axes[0].get_legend_handles_labels()
    handles.append(Line2D([0], [0], color="#555555", linestyle=":", linewidth=0.9))
    labels.append("Median")
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 1.02), ncol=3, fontsize=6.2, handlelength=2.4)
    fig.tight_layout(rect=(0, 0, 1, 0.88), pad=0.7, w_pad=1.2)

    base = outdir / "ks_density_bandwidth_sensitivity"
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(base.with_suffix(".tiff"), dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    write_tsv(
        base.with_suffix(".summary.tsv"),
        summary_rows,
        ["panel", "comparison", "n_valid_ks", "ks_median", "kde_peak", "bandwidth", "ks_filter", "display_range"],
    )
    write_tsv(
        base.with_suffix(".source_data.tsv"),
        source_rows,
        ["panel", "bandwidth", "comparison", "Ks", "density"],
    )
def main() -> int:
    args = parse_args()
    input_path = Path(args.input).resolve()
    if not input_path.is_file() or input_path.stat().st_size == 0:
        raise SystemExit(f"Missing input: {input_path}")
    if "block_representative" in str(input_path) and not args.allow_legacy_input:
        raise SystemExit(
            "Refusing legacy one-pair-per-block input. Use the filtered all-anchor result, "
            "or add --allow-legacy-input only for a QC preview."
        )
    if args.min_ks < 0 or args.max_ks <= args.min_ks:
        raise SystemExit("Require 0 <= min-ks < max-ks")
    if args.bandwidth <= 0 or args.grid_points < 101:
        raise SystemExit("Bandwidth must be positive and grid-points must be at least 101")

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    values_by_comparison, dropped = read_ks(input_path, args.min_ks, args.max_ks)
    manifest_rows = []
    for preset in args.presets:
        display_max = PRESET_DISPLAY_MAX[preset]
        plot_preset(
            preset,
            PRESETS[preset],
            values_by_comparison,
            args.min_ks,
            display_max,
            args.grid_points,
            args.bandwidth,
            outdir,
            args.width_mm,
            args.height_mm,
            args.min_ks,
            args.max_ks,
        )
        manifest_rows.append(
            {
                "output_stem": f"ks_density_{preset}",
                "preset": preset,
                "manuscript_role": PRESET_ROLES[preset],
                "bandwidth": f"{args.bandwidth:.8f}",
                "ks_filter": f"{args.min_ks:g}<Ks<{args.max_ks:g}",
                "display_range": f"{args.min_ks:g}-{display_max:g}",
            }
        )

    plot_bandwidth_sensitivity(
        values_by_comparison,
        args.bandwidth_sensitivity,
        outdir,
        args.min_ks,
        args.max_ks,
        args.grid_points,
    )
    if args.bandwidth_sensitivity:
        manifest_rows.append(
            {
                "output_stem": "ks_density_bandwidth_sensitivity",
                "preset": "focal",
                "manuscript_role": "supplementary_kde_bandwidth_robustness",
                "bandwidth": ",".join(f"{value:g}" for value in args.bandwidth_sensitivity),
                "ks_filter": f"{args.min_ks:g}<Ks<{args.max_ks:g}",
                "display_range": f"{args.min_ks:g}-{PRESET_DISPLAY_MAX['focal']:g}",
            }
        )

    plot_focal_with_references(
        values_by_comparison,
        args.bandwidth,
        outdir,
        args.min_ks,
        args.max_ks,
        args.grid_points,
    )
    manifest_rows.append(
        {
            "output_stem": "ks_density_focal_with_references",
            "preset": "focal,diploid_references",
            "manuscript_role": "main_or_supplementary_contextual_comparison",
            "bandwidth": f"{args.bandwidth:.8f}",
            "ks_filter": f"{args.min_ks:g}<Ks<{args.max_ks:g}",
            "display_range": "panel_a=0-0.5;panel_b=0-1",
        }
    )

    audit_rows = []
    for comparison in sorted(set(values_by_comparison) | set(dropped)):
        audit_rows.append(
            {
                "comparison": comparison,
                "n_plotted": len(values_by_comparison.get(comparison, [])),
                "n_dropped_or_outside_filter": dropped.get(comparison, 0),
            }
        )
    write_tsv(outdir / "ks_density_input_audit.tsv", audit_rows, ["comparison", "n_plotted", "n_dropped_or_outside_filter"])
    write_tsv(
        outdir / "formal_figure_manifest.tsv",
        manifest_rows,
        ["output_stem", "preset", "manuscript_role", "bandwidth", "ks_filter", "display_range"],
    )
    (outdir / "plot_parameters.txt").write_text(
        "\n".join(
            [
                f"input={input_path.name}",
                f"presets={','.join(args.presets)}",
                f"ks_filter={args.min_ks}<Ks<{args.max_ks}",
                f"kde=Gaussian with reflection at Ks=0",
                f"bandwidth={args.bandwidth}",
                f"bandwidth_sensitivity={','.join(str(value) for value in args.bandwidth_sensitivity)}",
                f"grid_points={args.grid_points}",
                "curve_scaling=probability_density",
                "vertical_dotted_line=median",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    (outdir / "COMPLETE").write_text("Formal Ks density figures completed.\n", encoding="ascii")
    print(f"[OK] Overlaid Ks density figures written to: {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
