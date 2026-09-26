#!/usr/bin/env python3
"""Render three complementary Figure 2 Ks distributions on one shared axis."""

from __future__ import annotations

import argparse
import csv
import math
from collections import OrderedDict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 7,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 0.8,
    "legend.frameon": False,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
})


PANEL_A = OrderedDict((key, label) for key, label in [
    ("PZI", r"$P.$ zippeliana A1+B1"),
    ("PSP", r"$P.$ spinulosa A1+B1"),
    ("PPA", r"$P.$ padus A+B"),
])
PANEL_B = OrderedDict((key, label) for key, label in [
    ("PPE", r"$P.$ persica"),
    ("PCAM", r"$P.$ campanulata"),
    ("GTR", r"$G.$ trifoliata"),
])
PANEL_C = OrderedDict((key, label) for key, label in [
    ("PZI", r"$P.$ zippeliana (A+B)-$G.$ trifoliata"),
    ("PSP", r"$P.$ spinulosa (A+B)-$G.$ trifoliata"),
    ("PPA", r"$P.$ padus (A+B)-$G.$ trifoliata"),
    ("PCAM", r"$P.$ campanulata-$G.$ trifoliata"),
    ("PPE", r"$P.$ persica-$G.$ trifoliata"),
])
COLORS = {
    "PZI": "#B64342", "PSP": "#0F4D92", "PPA": "#3D8C77",
    "PPE": "#76539B", "PCAM": "#C47B24", "GTR": "#4A4A4A",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("No header in {}".format(path))
        return list(reader)


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def parse_ks(row: dict[str, str], path: Path) -> float:
    try:
        value = float(row["Ks"])
    except (KeyError, TypeError, ValueError):
        raise ValueError("Invalid Ks row in {}: {}".format(path, row))
    if not math.isfinite(value) or not 0.0 < value < 5.0:
        raise ValueError("Unfiltered Ks value in {}: {}".format(path, value))
    return value


def reflected_kde(values: list[float], grid: np.ndarray, bandwidth: float) -> np.ndarray:
    sample = np.asarray(values, dtype=float)
    density = np.zeros_like(grid, dtype=float)
    denominator = len(sample) * bandwidth * math.sqrt(2.0 * math.pi)
    for start in range(0, len(sample), 4096):
        block = sample[start:start + 4096, None]
        direct = (grid[None, :] - block) / bandwidth
        reflected = (grid[None, :] + block) / bandwidth
        density += np.exp(-0.5 * direct * direct).sum(axis=0)
        density += np.exp(-0.5 * reflected * reflected).sum(axis=0)
    return density / denominator


def grouped_self(path: Path, expected_panel: str, labels: OrderedDict[str, str]) -> dict[str, list[float]]:
    grouped = {key: [] for key in labels}
    for row in read_tsv(path):
        if row.get("panel") != expected_panel:
            continue
        comparison = row.get("comparison", "")
        if comparison not in grouped:
            raise ValueError("Unexpected comparison {} in {}".format(comparison, path))
        grouped[comparison].append(parse_ks(row, path))
    for comparison, values in grouped.items():
        if not values:
            raise ValueError("No {} Ks values for {}".format(expected_panel, comparison))
    return grouped


def grouped_orthologue(path: Path) -> dict[str, list[float]]:
    grouped = {key: [] for key in PANEL_C}
    for row in read_tsv(path):
        comparison = row.get("comparison", "")
        if comparison not in grouped:
            raise ValueError("Unexpected orthologue comparison {} in {}".format(comparison, path))
        grouped[comparison].append(parse_ks(row, path))
    for comparison, values in grouped.items():
        if not values:
            raise ValueError("No orthologue Ks values for {}".format(comparison))
    return grouped


def mode_label(mode: float) -> str:
    return r"$K_s\approx0$" if mode == 0.0 else r"$K_s$={:.3f}".format(mode)


def plot_panel(ax: plt.Axes, panel_id: str, labels: OrderedDict[str, str], values_by_name: dict[str, list[float]],
               ylabel: str, grid: np.ndarray, bandwidth: float, annotation_offsets: dict[str, tuple[int, int]],
               density_rows: list[dict[str, object]], peak_rows: list[dict[str, object]],
               summary_rows: list[dict[str, object]]) -> None:
    curves: list[tuple[str, np.ndarray, int]] = []
    ymax = 0.0
    for name, display in labels.items():
        values = values_by_name[name]
        density = reflected_kde(values, grid, bandwidth)
        mode_index = int(np.argmax(density))
        ymax = max(ymax, float(density[mode_index]))
        curves.append((name, density, mode_index))
        ax.fill_between(grid, density, color=COLORS[name], alpha=0.08, linewidth=0)
        ax.plot(grid, density, color=COLORS[name], lw=1.25, label=display)
        for x, y in zip(grid, density):
            density_rows.append({
                "panel": panel_id, "comparison": name,
                "Ks_grid": "{:.8f}".format(float(x)), "density": "{:.8f}".format(float(y)),
                "n": len(values), "bandwidth": "{:.4f}".format(bandwidth),
            })
        peak_rows.append({
            "panel": panel_id, "comparison": name, "n": len(values),
            "KDE_mode_Ks": "{:.8f}".format(float(grid[mode_index])),
            "mode_density": "{:.8f}".format(float(density[mode_index])),
        })
        summary_rows.append({
            "panel": panel_id, "comparison": name, "n": len(values),
            "KDE_mode_Ks": "{:.8f}".format(float(grid[mode_index])),
            "mode_density": "{:.8f}".format(float(density[mode_index])),
            "bandwidth": "{:.4f}".format(bandwidth),
        })

    ax.set_ylim(0.0, ymax * 1.23)
    ax.set_ylabel(ylabel)
    ax.legend(loc="upper right", fontsize=5.5, handlelength=2.0, borderpad=0.18, labelspacing=0.28)
    ax.tick_params(direction="out", length=3, pad=2)
    for name, density, mode_index in curves:
        dx, dy = annotation_offsets[name]
        mode = float(grid[mode_index])
        ax.annotate(
            mode_label(mode),
            xy=(mode, float(density[mode_index])), xytext=(dx, dy), textcoords="offset points",
            ha="center", va="bottom", fontsize=5.2, color=COLORS[name],
            arrowprops={"arrowstyle": "-", "lw": 0.42, "color": COLORS[name], "alpha": 0.85},
            annotation_clip=False,
        )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-ks", type=Path, required=True)
    parser.add_argument("--gtr-orthologue-ks", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--bandwidth", type=float, default=0.035)
    parser.add_argument("--orthologue-bandwidth", type=float, default=0.040,
                        help="KDE bandwidth for the GTR-anchored orthologue panel")
    parser.add_argument("--xmax", type=float, default=3.0)
    parser.add_argument("--grid-points", type=int, default=3001,
                        help="Shared grid points; 3001 gives 0.001 Ks spacing at xmax=3")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.bandwidth <= 0.0 or args.orthologue_bandwidth <= 0.0 or args.xmax <= 0.0:
        raise SystemExit("[ERROR] bandwidths and --xmax must be positive")
    if args.outdir.exists() and any(args.outdir.iterdir()) and not args.overwrite:
        raise SystemExit("[ERROR] Output directory is not empty: {}".format(args.outdir))
    args.outdir.mkdir(parents=True, exist_ok=True)

    panel_a = grouped_self(args.self_ks, "upper", PANEL_A)
    panel_b = grouped_self(args.self_ks, "lower", PANEL_B)
    panel_c = grouped_orthologue(args.gtr_orthologue_ks)
    if args.grid_points < 101:
        raise SystemExit("[ERROR] --grid-points must be at least 101")
    grid = np.linspace(0.0, args.xmax, args.grid_points)
    density_rows: list[dict[str, object]] = []
    peak_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []

    fig, axes = plt.subplots(3, 1, figsize=(3.55, 6.75), sharex=True)
    plot_panel(
        axes[0], "a", PANEL_A, panel_a, "Paralogue density", grid, args.bandwidth,
        {"PZI": (0, 8), "PSP": (20, 22), "PPA": (-24, 8)},
        density_rows, peak_rows, summary_rows,
    )
    plot_panel(
        axes[1], "b", PANEL_B, panel_b, "Paralogue density", grid, args.bandwidth,
        {"PPE": (26, 16), "PCAM": (54, 30), "GTR": (0, 10)},
        density_rows, peak_rows, summary_rows,
    )
    plot_panel(
        axes[2], "c", PANEL_C, panel_c, "Orthologue density", grid, args.orthologue_bandwidth,
        {"PZI": (24, 10), "PSP": (44, 24), "PPA": (-44, 10), "PCAM": (-26, 24), "PPE": (-6, 38)},
        density_rows, peak_rows, summary_rows,
    )
    for ax, panel_id in zip(axes, ("a", "b", "c")):
        ax.text(-0.13, 1.02, panel_id, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom")
    axes[0].tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    axes[1].tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    axes[-1].set_xlim(0.0, args.xmax)
    axes[-1].set_xlabel(r"Synonymous substitutions per synonymous site ($K_s$)")
    fig.subplots_adjust(left=0.17, right=0.98, top=0.985, bottom=0.08, hspace=0.38)

    prefix = args.outdir / "figure2_integrated_ks"
    fig.savefig(str(prefix) + ".svg", bbox_inches="tight")
    fig.savefig(str(prefix) + ".pdf", bbox_inches="tight")
    fig.savefig(str(prefix) + ".tiff", dpi=600, bbox_inches="tight")
    fig.savefig(str(prefix) + ".png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    write_tsv(args.outdir / "integrated_density_grid.tsv",
              ["panel", "comparison", "Ks_grid", "density", "n", "bandwidth"], density_rows)
    write_tsv(args.outdir / "integrated_peak_labels.tsv",
              ["panel", "comparison", "n", "KDE_mode_Ks", "mode_density"], peak_rows)
    write_tsv(args.outdir / "integrated_plot_summary.tsv",
              ["panel", "comparison", "n", "KDE_mode_Ks", "mode_density", "bandwidth"], summary_rows)
    (args.outdir / "plot_parameters.tsv").write_text(
        "backend\tPython matplotlib\n"
        "kde\tGaussian KDE with reflection at Ks=0\n"
        "self_paralogue_bandwidth\t{:.4f}\n"
        "orthologue_bandwidth\t{:.4f}\n"
        "filter\t0<Ks<5 before density estimation\n"
        "shared_xmax\t{:.3f}\n"
        "median_lines\tnot_drawn; density modes are annotated\n"
        "panel_a\trepresentative A+B self-syntenic paralogues\n"
        "panel_b\tgenome-native cross-chromosome self-syntenic paralogues\n"
        "panel_c\tGTR-anchored collinear orthologues\n"
        "grid_points\t{}\n".format(args.bandwidth, args.orthologue_bandwidth, args.grid_points),
        encoding="ascii",
    )
    (args.outdir / "COMPLETE").write_text("Integrated three-panel Ks figure plotted.\n", encoding="ascii")
    print("[OK] Integrated Ks figure written: {}".format(args.outdir))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
