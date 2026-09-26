#!/usr/bin/env python3
"""Plot neutral internal synteny controls for the Figure 2 chromosome models.

The figure intentionally uses analysis labels only (copy1, copy2, hap01-hap04).
It visualizes (a) copy1-copy2 chromosome anchoring, (b) the independent evidence
used to reject a copy-specific fusion/fission call, and (c) all six anonymous
haplotype-pair comparisons.  It does not infer biological subgenome labels.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
import numpy as np


plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams.update({
    "font.size": 7,
    "axes.linewidth": 0.7,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
})


BLUE = "#0F4D92"
BLUE_LIGHT = "#8DB5D9"
GRAY = "#747474"
LIGHT = "#ECEFF2"
OUTLINE = "#303030"
HEATMAP = LinearSegmentedColormap.from_list("neutral_anchor", ["#FFFFFF", "#C8DDEC", BLUE])


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[1]
    default_summary = root / "downstream_analysis" / "figure2_polyploid" / "neutral_internal_synteny" / "summary"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", default=default_summary / "chromosome_matrix.tsv", type=Path)
    parser.add_argument("--resolution", default=default_summary / "copy1_copy2_gtr_resolution.tsv", type=Path)
    parser.add_argument("--hap-pairs", default=default_summary / "independent_haplotype_pairs.tsv", type=Path)
    parser.add_argument("--outdir", default=root / "downstream_analysis" / "figure2_polyploid" / "fig2_neutral_internal_synteny", type=Path)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def suffix(value: str) -> str:
    return value.split("_", 1)[1] if "_" in value else value


def natural_key(value: str) -> tuple:
    chunks: list[object] = []
    number = ""
    for char in value:
        if char.isdigit():
            number += char
        else:
            if number:
                chunks.append(int(number))
                number = ""
            chunks.append(char)
    if number:
        chunks.append(int(number))
    return tuple(chunks)


def label_panel(ax, label: str) -> None:
    ax.text(-0.16, 1.06, label, transform=ax.transAxes, fontsize=10,
            fontweight="bold", ha="left", va="bottom")


def main() -> int:
    args = parse_args()
    for path in (args.matrix, args.resolution, args.hap_pairs):
        if not path.is_file():
            raise SystemExit(f"Missing input: {path}")

    matrix_rows = [row for row in read_tsv(args.matrix) if row["comparison"] == "copy1_vs_copy2"]
    resolution = sorted(read_tsv(args.resolution), key=lambda row: natural_key(row["copy1_chr"]))
    hap_pairs = read_tsv(args.hap_pairs)
    if len(resolution) != 16:
        raise SystemExit(f"Expected 16 copy1 resolution rows, found {len(resolution)}")
    if len(hap_pairs) != 6:
        raise SystemExit(f"Expected six anonymous haplotype comparisons, found {len(hap_pairs)}")

    chromosomes = [row["copy1_chr"] for row in resolution]
    index = {chrom: i for i, chrom in enumerate(chromosomes)}
    heatmap = np.zeros((len(chromosomes), len(chromosomes)), dtype=float)
    for row in matrix_rows:
        query = suffix(row["query_chr"])
        target = suffix(row["reference_chr"])
        if query in index and target in index:
            heatmap[index[query], index[target]] = float(row["query_fraction"])

    primary_anchor = np.array([float(row["expected_target_anchor_fraction"]) for row in resolution])
    genome_alignment = np.array([float(row["same_chr_genome_alignment_fraction"]) for row in resolution])
    multiblock = np.array([len([block for block in row["gtr_supported_blocks"].split(",") if block and block != "NA"]) > 1
                           for row in resolution])
    breakpoint = np.array([int(row["copy_specific_breakpoint_candidate"]) for row in resolution])

    hap_pairs_data = []
    for row in hap_pairs:
        specificity = (float(row["query_top_specificity"]) + float(row["reference_top_specificity"])) / 2.0
        hap_pairs_data.append({"label": row["comparison"].replace("_vs_", "-"), "specificity": specificity})
    hap_pairs_data.sort(key=lambda row: row["specificity"])

    fig = plt.figure(figsize=(7.1, 3.9), constrained_layout=True)
    layout = fig.add_gridspec(1, 3, width_ratios=[1.22, 1.12, 0.86], wspace=0.42)
    ax_a = fig.add_subplot(layout[0, 0])
    ax_b = fig.add_subplot(layout[0, 1])
    ax_c = fig.add_subplot(layout[0, 2])

    # Panel a: full neutral chromosome anchor matrix.
    image = ax_a.imshow(heatmap, cmap=HEATMAP, vmin=0, vmax=max(0.70, float(heatmap.max())), aspect="equal", interpolation="nearest")
    ax_a.set_xticks(np.arange(16), [suffix(value) for value in chromosomes], rotation=90, fontsize=5.5)
    ax_a.set_yticks(np.arange(16), [suffix(value) for value in chromosomes], fontsize=5.5)
    ax_a.set_xlabel("copy2 chromosome")
    ax_a.set_ylabel("copy1 chromosome")
    ax_a.set_title("Neutral copy1-copy2 anchors", fontsize=8, pad=8)
    for i in range(16):
        ax_a.add_patch(plt.Rectangle((i - 0.5, i - 0.5), 1, 1, fill=False, ec=OUTLINE, lw=0.45))
    colorbar = fig.colorbar(image, ax=ax_a, fraction=0.043, pad=0.025)
    colorbar.set_label("row anchor fraction", fontsize=6)
    colorbar.ax.tick_params(labelsize=5.5, length=2)
    label_panel(ax_a, "a")

    # Panel b: both criteria used for a structural-breakpoint candidate.
    y = np.arange(16)
    ax_b.axvline(0.50, color="#B7B7B7", linestyle="--", lw=0.8, zorder=0)
    ax_b.scatter(genome_alignment, y, marker="s", s=21, color=BLUE_LIGHT, edgecolor="white", linewidth=0.35, zorder=2)
    ax_b.scatter(primary_anchor, y, marker="o", s=22, color=BLUE, edgecolor="white", linewidth=0.35, zorder=3)
    ax_b.scatter(primary_anchor[multiblock], y[multiblock], facecolors="none", edgecolors=OUTLINE,
                 marker="o", s=52, linewidth=0.75, zorder=4)
    if breakpoint.any():
        ax_b.scatter(primary_anchor[breakpoint], y[breakpoint], marker="x", s=32, color="#B64342", linewidth=1.1, zorder=5)
    ax_b.set_xlim(0.45, 1.02)
    ax_b.set_ylim(15.6, -0.6)
    ax_b.set_yticks(y, [suffix(value) for value in chromosomes], fontsize=5.5)
    ax_b.set_xlabel("fraction")
    ax_b.set_ylabel("copy1 chromosome")
    ax_b.set_title("Internal validation of GTR multi-block hits", fontsize=8, pad=8)
    ax_b.tick_params(axis="x", labelsize=6)
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=BLUE, markeredgecolor="white", markersize=5.5, label="primary anchors"),
        Line2D([0], [0], marker="s", color="none", markerfacecolor=BLUE_LIGHT, markeredgecolor="white", markersize=5.5, label="genome alignment"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="none", markeredgecolor=OUTLINE, markersize=7, label="GTR multi-block"),
    ]
    ax_b.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.0, 1.01), ncol=1,
                fontsize=5.6, handletextpad=0.35, borderaxespad=0)
    label_panel(ax_b, "b")

    # Panel c: data-driven pair strengths among independently assembled haplotypes.
    labels = [row["label"] for row in hap_pairs_data]
    values = np.array([row["specificity"] for row in hap_pairs_data])
    y_c = np.arange(len(labels))
    colors = [BLUE if value >= 0.90 else GRAY for value in values]
    ax_c.axvline(0.90, color="#B7B7B7", linestyle="--", lw=0.8, zorder=0)
    ax_c.hlines(y_c, 0.75, values, color="#D4D9DD", lw=1.1, zorder=1)
    ax_c.scatter(values, y_c, s=29, color=colors, edgecolor="white", linewidth=0.35, zorder=2)
    ax_c.set_xlim(0.75, 1.01)
    ax_c.set_ylim(len(labels) - 0.55, -0.55)
    ax_c.set_yticks(y_c, labels, fontsize=5.8)
    ax_c.set_xlabel("mean chromosome specificity")
    ax_c.set_title("Independent haplotype pairs", fontsize=8, pad=8)
    ax_c.tick_params(axis="x", labelsize=6)
    label_panel(ax_c, "c")

    for axis in (ax_a, ax_b, ax_c):
        axis.tick_params(length=2, width=0.6)
    fig.text(0.50, 0.002,
             "Neutral labels only. A breakpoint candidate requires both a disrupted primary internal target and weak same-chromosome alignment.",
             ha="center", va="bottom", fontsize=5.6, color="#444444")

    args.outdir.mkdir(parents=True, exist_ok=True)
    base = args.outdir / "figure2_neutral_internal_synteny_validation"
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    fig.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)

    with (args.outdir / "source_data_manifest.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["panel", "source_table", "quantity"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows([
            {"panel": "a", "source_table": str(args.matrix), "quantity": "row-normalized copy1-copy2 anchor fractions"},
            {"panel": "b", "source_table": str(args.resolution), "quantity": "primary-copy2 anchor fraction, genome alignment fraction and GTR multiblock flag"},
            {"panel": "c", "source_table": str(args.hap_pairs), "quantity": "mean of reciprocal top chromosome specificity"},
        ])
    print(f"[OK] Wrote neutral internal synteny validation figure: {base}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
