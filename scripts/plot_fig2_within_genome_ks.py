#!/usr/bin/env python3
"""Render a compact two-row within-genome Ks diagnostic figure in Python."""

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


UPPER = OrderedDict((key, value) for key, value in [
    ("PZI", r"$P.$ zippeliana A1+B1"),
    ("PSP", r"$P.$ spinulosa A1+B1"),
    ("PPA", r"$P.$ padus A+B"),
])
LOWER = OrderedDict((key, value) for key, value in [
    ("PPE", r"$P.$ persica"),
    ("PCAM", r"$P.$ campanulata"),
    ("GTR", r"$G.$ trifoliata"),
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


def reflected_kde(values: list[float], grid: np.ndarray, bandwidth: float) -> np.ndarray:
    sample = np.asarray(values, dtype=float)
    density = np.zeros_like(grid, dtype=float)
    denom = len(sample) * bandwidth * math.sqrt(2.0 * math.pi)
    for start in range(0, len(sample), 4096):
        block = sample[start:start + 4096, None]
        direct = (grid[None, :] - block) / bandwidth
        reflected = (grid[None, :] + block) / bandwidth
        density += np.exp(-0.5 * direct * direct).sum(axis=0)
        density += np.exp(-0.5 * reflected * reflected).sum(axis=0)
    return density / denom


def grouped_values(source: Path) -> dict[str, dict[str, list[float]]]:
    grouped = {"upper": {key: [] for key in UPPER}, "lower": {key: [] for key in LOWER}}
    for row in read_tsv(source):
        panel = row.get("panel", "")
        comparison = row.get("comparison", "")
        if panel not in grouped or comparison not in grouped[panel]:
            raise SystemExit("[ERROR] Unexpected source row: {}".format(row))
        try:
            value = float(row["Ks"])
        except (KeyError, TypeError, ValueError):
            raise SystemExit("[ERROR] Invalid Ks source row: {}".format(row))
        if not math.isfinite(value) or not 0.0 < value < 5.0:
            raise SystemExit("[ERROR] Source contains unfiltered Ks: {}".format(value))
        grouped[panel][comparison].append(value)
    for panel, values_by_name in grouped.items():
        for name, values in values_by_name.items():
            if not values:
                raise SystemExit("[ERROR] No values for {} {}".format(panel, name))
    return grouped


def draw_panel(ax, values_by_name: dict[str, list[float]], labels: OrderedDict[str, str], panel: str,
               grid: np.ndarray, bandwidth: float, density_rows: list[dict[str, object]]) -> None:
    ymax = 0.0
    curves: list[tuple[str, np.ndarray, int, float]] = []
    for name, display in labels.items():
        values = values_by_name[name]
        density = reflected_kde(values, grid, bandwidth)
        ymax = max(ymax, float(density.max()))
        ax.fill_between(grid, density, color=COLORS[name], alpha=0.09, linewidth=0)
        ax.plot(grid, density, color=COLORS[name], lw=1.25, label=display)
        mode_index = int(np.argmax(density))
        curves.append((name, density, mode_index, float(grid[mode_index])))
        for x, y in zip(grid, density):
            density_rows.append({"panel": panel, "comparison": name, "Ks_grid": "{:.8f}".format(float(x)),
                                 "density": "{:.8f}".format(float(y)), "n": len(values),
                                 "bandwidth": "{:.4f}".format(bandwidth)})
    ax.set_ylabel("Paralogue density")
    ax.set_ylim(0.0, ymax * 1.10)
    # Label the density modes directly.  Offsets are deliberately staggered so
    # the closely spaced A-B peaks and the boundary modes remain readable.
    offsets = {
        "upper": {"PZI": (0, 10), "PSP": (18, 22), "PPA": (-22, 10)},
        "lower": {"PPE": (26, 16), "PCAM": (52, 30), "GTR": (0, 10)},
    }
    for name, density, mode_index, mode in curves:
        dx, dy = offsets.get(panel, {}).get(name, (0, 10))
        mode_label = r"$K_s\approx0$" if mode == 0.0 else r"$K_s$={:.3f}".format(mode)
        ax.annotate(
            mode_label,
            xy=(float(grid[mode_index]), float(density[mode_index])),
            xytext=(dx, dy),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=5.5,
            color=COLORS[name],
            arrowprops={"arrowstyle": "-", "lw": 0.45, "color": COLORS[name], "alpha": 0.85},
            annotation_clip=False,
        )
    ax.legend(loc="upper right", fontsize=5.8, handlelength=2.0, borderpad=0.15, labelspacing=0.34)
    ax.tick_params(direction="out", length=3, pad=2)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--bandwidth", type=float, default=0.035)
    parser.add_argument("--xmax", type=float, default=1.0)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.bandwidth <= 0 or args.xmax <= 0:
        raise SystemExit("[ERROR] --bandwidth and --xmax must be positive")
    if args.outdir.exists() and any(args.outdir.iterdir()) and not args.overwrite:
        raise SystemExit("[ERROR] Output directory is not empty: {}".format(args.outdir))
    args.outdir.mkdir(parents=True, exist_ok=True)
    grouped = grouped_values(args.source)
    grid = np.linspace(0.0, args.xmax, 1201)
    density_rows: list[dict[str, object]] = []

    fig, axes = plt.subplots(2, 1, figsize=(3.35, 4.45), sharex=True)
    draw_panel(axes[0], grouped["upper"], UPPER, "upper", grid, args.bandwidth, density_rows)
    draw_panel(axes[1], grouped["lower"], LOWER, "lower", grid, args.bandwidth, density_rows)
    axes[0].tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    axes[1].set_xlabel(r"$K_s$")
    axes[0].text(-0.13, 1.02, "a", transform=axes[0].transAxes, fontsize=9, fontweight="bold", va="bottom")
    axes[1].text(-0.13, 1.02, "b", transform=axes[1].transAxes, fontsize=9, fontweight="bold", va="bottom")
    axes[1].set_xlim(0.0, args.xmax)
    fig.subplots_adjust(left=0.16, right=0.98, top=0.98, bottom=0.11, hspace=0.34)
    prefix = args.outdir / "figure2_within_genome_ks"
    fig.savefig(str(prefix) + ".svg", bbox_inches="tight")
    fig.savefig(str(prefix) + ".pdf", bbox_inches="tight")
    fig.savefig(str(prefix) + ".tiff", dpi=600, bbox_inches="tight")
    fig.savefig(str(prefix) + ".png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    write_tsv(args.outdir / "within_genome_density_grid.tsv",
              ["panel", "comparison", "Ks_grid", "density", "n", "bandwidth"], density_rows)
    summary = []
    for panel, mapping in (("upper", UPPER), ("lower", LOWER)):
        for name in mapping:
            values = grouped[panel][name]
            density = reflected_kde(values, grid, args.bandwidth)
            summary.append({"panel": panel, "comparison": name, "n": len(values),
                            "Ks_median": "{:.8f}".format(float(np.median(values))),
                            "KDE_mode": "{:.8f}".format(float(grid[int(np.argmax(density))])),
                            "bandwidth": "{:.4f}".format(args.bandwidth), "xmax": "{:.3f}".format(args.xmax)})
    write_tsv(args.outdir / "within_genome_plot_summary.tsv",
              ["panel", "comparison", "n", "Ks_median", "KDE_mode", "bandwidth", "xmax"], summary)
    (args.outdir / "plot_parameters.tsv").write_text(
        "backend\tPython matplotlib\n"
        "kde\tGaussian KDE with reflection at Ks=0\n"
        "bandwidth\t{:.4f}\n"
        "filter\t0<Ks<5 before density estimation\n"
        "xmax\t{:.3f}\n"
        "median_lines\tnot_drawn; peak location is the relevant visual diagnostic\n".format(args.bandwidth, args.xmax),
        encoding="ascii")
    (args.outdir / "COMPLETE").write_text("Within-genome Ks density figure plotted.\n", encoding="ascii")
    print("[OK] Within-genome Ks plot written: {}".format(args.outdir))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
