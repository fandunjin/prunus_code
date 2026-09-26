#!/usr/bin/env python3
"""Plot A/B subgenome Ks distributions against PCAM, PPE and GTR."""

from __future__ import annotations

import argparse
import csv
import math
from collections import OrderedDict, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
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


COLORS = {"PZI": "#B64342", "PSP": "#0F4D92", "PPA": "#42949E"}
DISPLAY = {"PZI": "P. zippeliana", "PSP": "P. spinulosa", "PPA": "P. padus"}
QUERY = OrderedDict((
    ("PZI", ("pzi_hapA1", "pzi_hapB1")),
    ("PSP", ("psp_hapA1", "psp_hapB1")),
    ("PPA", ("ppa_hapA", "ppa_hapB")),
))
REFERENCES = OrderedDict((
    ("pcam", "P. campanulata"),
    ("ppe", "P. persica"),
    ("gtr", "G. trifoliata"),
))
STYLES = {"A": "-", "B": (0, (4.0, 2.0))}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def reflected_gaussian_kde(values: np.ndarray, grid: np.ndarray, bandwidth: float) -> np.ndarray:
    density = np.zeros_like(grid, dtype=float)
    denominator = len(values) * bandwidth * math.sqrt(2.0 * math.pi)
    for start in range(0, len(values), 4096):
        chunk = values[start:start + 4096, None]
        direct = (grid[None, :] - chunk) / bandwidth
        reflected = (grid[None, :] + chunk) / bandwidth
        density += np.exp(-0.5 * direct * direct).sum(axis=0)
        density += np.exp(-0.5 * reflected * reflected).sum(axis=0)
    return density / denominator


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--bandwidth", type=float, default=0.04)
    parser.add_argument("--xmax", type=float, default=1.0)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.outdir.exists() and any(args.outdir.iterdir()) and not args.overwrite:
        raise SystemExit(f"[ERROR] Output directory is not empty: {args.outdir}")
    args.outdir.mkdir(parents=True, exist_ok=True)
    obsolete_medians = args.outdir / "median_lines.tsv"
    if args.overwrite and obsolete_medians.is_file():
        obsolete_medians.unlink()

    grouped: dict[str, list[float]] = defaultdict(list)
    for row in read_tsv(args.input):
        comparison = row.get("comparison", "")
        try:
            value = float(row.get("Ks", ""))
        except (TypeError, ValueError):
            continue
        if math.isfinite(value) and 0.0 < value < 5.0:
            grouped[comparison].append(value)

    expected = {
        f"{query}_vs_{reference}"
        for reference in REFERENCES
        for query_pair in QUERY.values()
        for query in query_pair
    }
    missing = sorted(expected - set(grouped))
    empty = sorted(comparison for comparison in expected if not grouped.get(comparison))
    if missing or empty:
        raise SystemExit(f"[ERROR] Missing or empty comparisons: {sorted(set(missing + empty))}")

    grid = np.linspace(0.0, args.xmax, 1001)
    density_by_comparison: dict[str, np.ndarray] = {}
    for comparison in sorted(expected):
        density_by_comparison[comparison] = reflected_gaussian_kde(
            np.asarray(grouped[comparison], dtype=float), grid, args.bandwidth
        )
    global_peak = max(float(np.max(density)) for density in density_by_comparison.values())
    y_max = max(4.0, math.ceil(global_peak * 1.12 * 2.0) / 2.0)

    fig, axes = plt.subplots(3, 1, figsize=(3.6, 8.3), sharex=True, sharey=True)
    source_rows: list[dict[str, object]] = []
    peak_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []

    for panel_index, (reference, reference_label) in enumerate(REFERENCES.items()):
        ax = axes[panel_index]
        ax.set_xlim(0.0, args.xmax)
        ax.set_ylim(0.0, y_max)
        ax.grid(axis="y", color="#D9D9D9", lw=0.35, alpha=0.55)
        ax.set_axisbelow(True)
        ax.set_ylabel("Orthologue density")
        ax.text(-0.19, 1.015, chr(ord("a") + panel_index), transform=ax.transAxes,
                fontsize=10, fontweight="bold", va="bottom")
        ax.text(0.02, 0.96, reference_label, transform=ax.transAxes,
                fontsize=7, fontstyle="italic", va="top")
        curves = []

        for species, (query_a, query_b) in QUERY.items():
            for subgenome, query in (("A", query_a), ("B", query_b)):
                comparison = f"{query}_vs_{reference}"
                values = np.asarray(grouped[comparison], dtype=float)
                density = density_by_comparison[comparison]
                color = COLORS[species]
                label = f"{DISPLAY[species]} {subgenome}  n={len(values):,}"
                ax.fill_between(grid, density, color=color, alpha=0.035, linewidth=0)
                ax.plot(grid, density, color=color, lw=1.05, ls=STYLES[subgenome], label=label)
                mode_index = int(np.argmax(density))
                curves.append((species, subgenome, comparison, color, density, mode_index))
                peak_rows.append({
                    "panel": chr(ord("a") + panel_index),
                    "reference": reference.upper(),
                    "comparison": comparison,
                    "species": species,
                    "subgenome": subgenome,
                    "KDE_mode_Ks": f"{float(grid[mode_index]):.8f}",
                    "mode_density": f"{float(density[mode_index]):.8f}",
                    "n": len(values),
                })
                summary_rows.append({
                    "reference": reference.upper(),
                    "comparison": comparison,
                    "species": species,
                    "subgenome": subgenome,
                    "n": len(values),
                    "KDE_mode_Ks": f"{float(grid[mode_index]):.8f}",
                    "mode_density": f"{float(density[mode_index]):.8f}",
                    "bandwidth": f"{args.bandwidth:.4f}",
                })
                for x, y in zip(grid, density):
                    source_rows.append({
                        "panel": chr(ord("a") + panel_index),
                        "reference": reference.upper(),
                        "comparison": comparison,
                        "species": species,
                        "subgenome": subgenome,
                        "Ks_grid": f"{float(x):.8f}",
                        "density": f"{float(y):.8f}",
                        "n": len(values),
                        "bandwidth": f"{args.bandwidth:.4f}",
                    })

        # All six peak locations are shown with short, colour-matched leader
        # labels. A/B line style identifies the subgenome while the text gives
        # the KDE mode used as the plotted peak estimate.
        offsets = {
            "PZI_A": (-34, 16), "PZI_B": (-26, 32),
            "PSP_A": (0, 46), "PSP_B": (18, 60),
            "PPA_A": (34, 16), "PPA_B": (44, 32),
        }
        for species, subgenome, comparison, color, density, mode_index in curves:
            dx, dy = offsets[species + "_" + subgenome]
            ax.annotate(
                r"$K_s$={:.3f}".format(float(grid[mode_index])),
                xy=(float(grid[mode_index]), float(density[mode_index])),
                xytext=(dx, dy), textcoords="offset points",
                ha="center", va="bottom", fontsize=4.5, color=color,
                arrowprops={"arrowstyle": "-", "lw": 0.40, "color": color, "alpha": 0.82},
                annotation_clip=False,
            )

        ax.legend(loc="upper right", fontsize=4.9,
                  handlelength=2.0, borderpad=0.2, labelspacing=0.3)
        if panel_index < len(axes) - 1:
            ax.tick_params(axis="x", which="both", bottom=False, labelbottom=False)

    axes[-1].set_xlabel(r"Synonymous substitutions per synonymous site ($K_s$)")
    fig.subplots_adjust(left=0.18, right=0.98, top=0.985, bottom=0.075, hspace=0.34)
    prefix = args.outdir / "figure2a_subgenome_reference_ks"
    fig.savefig(str(prefix) + ".svg", bbox_inches="tight")
    fig.savefig(str(prefix) + ".pdf", bbox_inches="tight")
    fig.savefig(str(prefix) + ".tiff", dpi=600, bbox_inches="tight")
    fig.savefig(str(prefix) + ".png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    write_tsv(args.outdir / "density_grid.tsv",
              ["panel", "reference", "comparison", "species", "subgenome", "Ks_grid", "density", "n", "bandwidth"],
              source_rows)
    write_tsv(args.outdir / "peak_labels.tsv",
              ["panel", "reference", "comparison", "species", "subgenome", "KDE_mode_Ks", "mode_density", "n"],
              peak_rows)
    write_tsv(args.outdir / "plot_summary.tsv",
              ["reference", "comparison", "species", "subgenome", "n", "KDE_mode_Ks", "mode_density", "bandwidth"],
              summary_rows)
    (args.outdir / "plot_parameters.tsv").write_text(
        "backend\tPython matplotlib\n"
        "kde\tGaussian KDE with reflection at Ks=0\n"
        f"bandwidth\t{args.bandwidth:.4f}\n"
        "filter\t0<Ks<5 before density estimation\n"
        "median_lines\tnot_drawn; medians are descriptive statistics only\n"
        f"shared_xmax\t{args.xmax:.3f}\n"
        f"shared_ymax\t{y_max:.3f}\n"
        "line_encoding\tA solid; B dashed\n",
        encoding="ascii",
    )
    (args.outdir / "COMPLETE").write_text("Three-reference subgenome Ks figure completed.\n", encoding="ascii")
    print(f"[OK] Wrote three-reference A/B Ks figure: {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
