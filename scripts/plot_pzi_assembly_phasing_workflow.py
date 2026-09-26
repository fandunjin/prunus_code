#!/usr/bin/env python3
"""Draw the P. zippeliana chromosome assembly and haplotype phasing workflow."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle
from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures" / "pzi_assembly_phasing_workflow"

COLORS = {
    "ink": "#18212B",
    "muted": "#5F6B76",
    "line": "#AAB4BC",
    "pale": "#F3F6F7",
    "ont": "#177E89",
    "hic": "#D9822B",
    "assembly": "#637381",
    "phase1": "#2A6FBB",
    "phase2": "#D05A65",
    "a1": "#2A6FBB",
    "a2": "#43A6A5",
    "b1": "#E28A32",
    "b2": "#C94C78",
    "green": "#3B8E5A",
    "gold": "#C99A2E",
}


mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 9,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "axes.linewidth": 0.6,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    }
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def arrow(ax, x1, y1, x2, y2, color=None, lw=1.4, style="-|>", zorder=3):
    ax.add_patch(
        FancyArrowPatch(
            (x1, y1),
            (x2, y2),
            arrowstyle=style,
            mutation_scale=11,
            linewidth=lw,
            color=color or COLORS["ink"],
            shrinkA=0,
            shrinkB=0,
            zorder=zorder,
        )
    )


def rounded_box(ax, x, y, w, h, edge, face="white", lw=1.2, radius=0.7, zorder=1):
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle=f"round,pad=0.12,rounding_size={radius}",
        linewidth=lw,
        edgecolor=edge,
        facecolor=face,
        zorder=zorder,
    )
    ax.add_patch(patch)
    return patch


def label(ax, x, y, text, size=8, weight="normal", color=None, ha="center", va="center", **kwargs):
    return ax.text(
        x,
        y,
        text,
        fontsize=size,
        fontweight=weight,
        color=color or COLORS["ink"],
        ha=ha,
        va=va,
        linespacing=1.15,
        **kwargs,
    )


def draw_reads(ax, x, y, color, n=9, width=5.8, seed=7, lw=2.2):
    rng = np.random.default_rng(seed)
    for i in range(n):
        yy = y + i * 0.72
        start = x + rng.uniform(0, 1.8)
        length = rng.uniform(width * 0.48, width)
        ax.plot([start, start + length], [yy, yy], color=color, lw=lw, solid_capstyle="round")


def draw_contigs(ax, x, y, colors, n=10, width=6.4, seed=17):
    rng = np.random.default_rng(seed)
    for i in range(n):
        yy = y + i * 0.64
        start = x + rng.uniform(0, 2.2)
        length = rng.uniform(width * 0.45, width)
        ax.plot(
            [start, start + length],
            [yy, yy],
            color=colors[i % len(colors)],
            lw=3.0,
            solid_capstyle="butt",
        )


def draw_chromosome_stack(ax, x, y, colors, rows=8, copies=4, maxw=8.0, row_gap=0.9, copy_gap=0.2):
    rng = np.random.default_rng(11)
    widths = np.linspace(maxw, maxw * 0.58, rows)
    rng.shuffle(widths)
    for r in range(rows):
        for c in range(copies):
            yy = y + r * row_gap + c * copy_gap
            ax.plot(
                [x, x + widths[r]],
                [yy, yy],
                color=colors[c],
                lw=2.4,
                solid_capstyle="round",
            )


def draw_final_haplotypes(ax, x, y, w=10.0):
    sets = [
        ("A1", COLORS["a1"]),
        ("A2", COLORS["a2"]),
        ("B1", COLORS["b1"]),
        ("B2", COLORS["b2"]),
    ]
    for i, (name, color) in enumerate(sets):
        xx = x + i * (w / 4)
        label(ax, xx + 0.75, y + 6.0, name, 7.2, "bold", color=color)
        for r in range(8):
            length = 1.55 - 0.06 * ((r * 3) % 6)
            yy = y + 5.2 - r * 0.58
            ax.plot([xx, xx + length], [yy, yy], color=color, lw=2.0, solid_capstyle="round")


def add_inset(fig, parent_ax, rect_data):
    """Create an inset using coordinates from the 0-100 x 0-56 parent axis."""
    x, y, w, h = rect_data
    bbox = parent_ax.get_position()
    return fig.add_axes(
        [
            bbox.x0 + bbox.width * x / 100,
            bbox.y0 + bbox.height * y / 56,
            bbox.width * w / 100,
            bbox.height * h / 56,
        ]
    )


def reference_alignment_matrix():
    rows = read_tsv(ROOT / "haploid" / "evidence" / "pzi_vs_pcam_alignment_summary.tsv")
    queries = [f"chr{i:02d}" for i in range(1, 33)]
    refs = sorted({row["reference_chr"] for row in rows})
    qindex = {name: i for i, name in enumerate(queries)}
    rindex = {name: i for i, name in enumerate(refs)}
    matrix = np.zeros((len(refs), len(queries)), dtype=float)
    for row in rows:
        if row["query_chr"] in qindex:
            matrix[rindex[row["reference_chr"]], qindex[row["query_chr"]]] = float(row["aligned_block_bp"])
    return np.log10(matrix + 1)


def summary_metrics():
    summary = {row["metric"]: row["value"] for row in read_tsv(ROOT / "final" / "tables" / "pzi_2n32.final_genome_summary.tsv")}
    validation = {row["metric"]: row for row in read_tsv(ROOT / "haploid" / "final_haplotypes" / "FINAL_VALIDATION.tsv")}
    telomeres = read_tsv(ROOT / "final" / "tables" / "tidk_chr_end_summary.tsv")
    telomere_ends = sum(int(row["telomere_ends"]) for row in telomeres)
    hap_bp = sum(
        int(validation[f"hap{name}_fasta_bp"]["value"])
        for name in ("A1", "A2", "B1", "B2")
    )
    return summary, validation, telomere_ends, hap_bp


def draw():
    fig = plt.figure(figsize=(15.2, 8.4))
    ax = fig.add_axes([0.025, 0.035, 0.95, 0.93])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 56)
    ax.axis("off")

    # Title and phase bands.
    label(
        ax,
        2.0,
        54.5,
        "Chromosome-scale assembly and haplotype phasing of Prunus zippeliana",
        14.5,
        "bold",
        ha="left",
    )
    label(
        ax,
        98.0,
        54.5,
        "ONT + Hi-C  |  x = 8  |  2n = 4x = 32",
        8.4,
        color=COLORS["muted"],
        ha="right",
    )
    bands = [
        (1.0, 29.0, "I  Assembly", "#EAF3F3", COLORS["ont"]),
        (30.0, 31.0, "II  Chromosome scaffolding", "#F8F0E6", COLORS["hic"]),
        (62.0, 23.0, "III  Evidence-based phasing", "#EEF2F8", COLORS["phase1"]),
        (86.0, 13.0, "IV  Freeze and QC", "#EEF5EF", COLORS["green"]),
    ]
    for x, w, text, face, edge in bands:
        ax.add_patch(Rectangle((x, 50.7), w, 1.7, facecolor=face, edgecolor="none", zorder=0))
        label(ax, x + 0.8, 51.55, text, 7.2, "bold", color=edge, ha="left")

    base_y = 31.8
    box_h = 16.7
    stages = [
        (1.0, 12.8, COLORS["ont"]),
        (16.1, 12.8, COLORS["assembly"]),
        (31.3, 13.2, COLORS["hic"]),
        (47.0, 13.2, COLORS["hic"]),
        (62.8, 20.8, COLORS["phase1"]),
        (86.0, 13.0, COLORS["green"]),
    ]
    for x, w, edge in stages:
        rounded_box(ax, x, base_y, w, box_h, edge=edge, face="white", lw=1.15, radius=0.55)

    # Step 1: sequencing.
    label(ax, 7.4, 46.9, "1  Sequencing", 9.0, "bold")
    draw_reads(ax, 2.3, 40.0, COLORS["ont"], n=7, width=8.2, seed=3, lw=2.1)
    label(ax, 7.3, 38.7, "ONT long reads", 7.4, "bold", color=COLORS["ont"])
    for i in range(5):
        yy = 34.2 + i * 0.58
        ax.plot([2.8, 5.0], [yy, yy], color=COLORS["hic"], lw=2.0)
        ax.plot([9.4, 11.6], [yy, yy], color=COLORS["hic"], lw=2.0)
        ax.plot([5.0, 9.4], [yy, yy + (0.6 if i % 2 == 0 else -0.5)], color="#E7C59E", lw=0.8)
    label(ax, 7.3, 33.0, "Hi-C read pairs", 7.4, "bold", color=COLORS["hic"])

    # Step 2: hifiasm.
    label(ax, 22.5, 46.9, "2  Primary assembly", 9.0, "bold")
    label(ax, 22.5, 44.7, "hifiasm --ont + Hi-C", 7.5, "bold", color=COLORS["assembly"])
    draw_contigs(ax, 17.5, 37.1, [COLORS["phase1"], "#89A9CF"], n=8, width=9.0, seed=9)
    label(ax, 22.5, 35.9, "hap1 contigs", 7.2, color=COLORS["phase1"])
    draw_contigs(ax, 17.5, 32.7, [COLORS["phase2"], "#DB9CA3"], n=4, width=9.0, seed=13)
    label(ax, 22.5, 32.2, "hap2 contigs", 7.2, color=COLORS["phase2"])

    # Step 3: merge and HapHiC.
    label(ax, 37.9, 46.9, "3  Merge and scaffold", 9.0, "bold")
    label(ax, 37.9, 44.8, "merge hap1 + hap2", 7.5, "bold")
    arrow(ax, 34.0, 43.1, 41.6, 43.1, color=COLORS["line"], lw=1.0)
    ax.plot([34.2, 36.6], [42.2, 42.2], color=COLORS["phase1"], lw=3.2)
    ax.plot([39.0, 41.4], [42.2, 42.2], color=COLORS["phase2"], lw=3.2)
    label(ax, 37.9, 39.9, "HapHiC (n = 32)", 8.0, "bold", color=COLORS["hic"])
    for i in range(10):
        y = 37.9 - i * 0.48
        length = 8.3 - (i % 4) * 0.55
        ax.plot([33.6, 33.6 + length], [y, y], color="#7F8A93", lw=1.7, solid_capstyle="round")
    label(ax, 37.9, 32.6, "Hi-C clustering, ordering\nand orientation", 7.1, color=COLORS["muted"])

    # Step 4: manual curation and chromosome pairs.
    label(ax, 53.6, 46.9, "4  Contact-map curation", 9.0, "bold")
    label(ax, 53.6, 44.7, "Juicebox / HapHiC review", 7.4, "bold", color=COLORS["hic"])
    for p in range(8):
        yy = 42.4 - p * 1.08
        ax.plot([49.0, 52.7], [yy, yy], color=COLORS["phase1"], lw=2.3, solid_capstyle="round")
        ax.plot([54.4, 58.1], [yy, yy], color=COLORS["phase2"], lw=2.3, solid_capstyle="round")
        ax.plot([52.85, 54.25], [yy, yy], color=COLORS["gold"], lw=0.8, ls=(0, (2, 2)))
    label(ax, 53.6, 33.3, "32 chromosomes arranged as\n16 homolog-linked pairs", 7.4, "bold")

    # Step 5: multi-evidence phasing.
    label(ax, 73.2, 46.9, "5  Multi-evidence phasing", 9.0, "bold")
    bullets = [
        ("map chr32 back to hifiasm hap1 / hap2", COLORS["phase1"]),
        ("corrected homolog-pair support", COLORS["hic"]),
        ("self-alignment and segment profiles", COLORS["gold"]),
        ("multi-reference and outgroup synteny", COLORS["green"]),
    ]
    for i, (text, color) in enumerate(bullets):
        yy = 43.8 - i * 1.75
        ax.add_patch(Rectangle((64.4, yy - 0.25), 0.48, 0.48, facecolor=color, edgecolor="none"))
        label(ax, 65.4, yy, text, 7.1, color=COLORS["ink"], ha="left")
    label(ax, 73.2, 35.9, "Global exact-cover selection", 7.7, "bold", color=COLORS["phase1"])
    for i, color in enumerate([COLORS["a1"], COLORS["a2"], COLORS["b1"], COLORS["b2"]]):
        ax.add_patch(Rectangle((67.0 + i * 3.2, 33.1), 2.6, 1.25, facecolor=color, edgecolor="white", lw=0.4))
    label(ax, 73.2, 32.4, "8 homologous quartets", 7.2, "bold")

    # Step 6: final four haplotype sets.
    label(ax, 92.5, 46.9, "6  Final phased sets", 9.0, "bold")
    draw_final_haplotypes(ax, 87.5, 37.8, w=10.0)
    label(ax, 92.5, 36.3, "4 phased sets x 8 chr", 7.6, "bold")
    label(ax, 92.5, 34.1, "A1 / A2 / B1 / B2", 7.4, color=COLORS["muted"])
    rounded_box(ax, 88.0, 32.3, 9.0, 1.1, edge=COLORS["green"], face="#EEF5EF", lw=0.8, radius=0.3)
    label(ax, 92.5, 32.85, "FINAL FREEZE", 6.7, "bold", color=COLORS["green"])

    # Main arrows and branch to evidence.
    for x1, x2 in [(13.8, 16.1), (28.9, 31.3), (44.5, 47.0), (60.2, 62.8), (83.6, 86.0)]:
        arrow(ax, x1, 40.15, x2, 40.15, color=COLORS["ink"], lw=1.25)

    # Evidence strip background and headings.
    ax.add_patch(Rectangle((1.0, 2.0), 98.0, 26.7, facecolor="#FAFBFB", edgecolor="#D7DEE2", lw=0.8, zorder=0))
    label(ax, 2.0, 27.2, "Representative evidence used for curation, phasing and release", 8.8, "bold", ha="left")

    # Inset 1: real Hi-C matrix.
    hic_ax = add_inset(fig, ax, (3.0, 7.0, 19.2, 18.2))
    hic_img = Image.open(ROOT / "final" / "qc" / "hic" / "chr_pair_contact_heatmap.png").convert("RGB")
    hic_ax.imshow(hic_img)
    hic_ax.set_xlim(300, 2700)
    hic_ax.set_ylim(2250, 120)
    hic_ax.axis("off")
    label(ax, 12.6, 5.4, "Hi-C: strong chromosome-wise cis signal", 7.4, "bold")
    arrow(ax, 53.6, 31.7, 12.6, 26.1, color=COLORS["hic"], lw=0.85, style="-|>", zorder=2)

    # Inset 2: final quartet assignment table, generated directly from the scored output.
    self_ax = add_inset(fig, ax, (26.0, 7.0, 16.5, 18.2))
    quartet_rows = read_tsv(ROOT / "haploid" / "final_haplotypes" / "FINAL_HAPLOTYPE_OVERVIEW.tsv")
    quartet_cols = [("hapA1", "A1", COLORS["a1"]), ("hapA2", "A2", COLORS["a2"]), ("hapB1", "B1", COLORS["b1"]), ("hapB2", "B2", COLORS["b2"])]
    self_ax.set_xlim(0, 4)
    self_ax.set_ylim(8, 0)
    for col, (_, short, color) in enumerate(quartet_cols):
        self_ax.text(col + 0.5, -0.35, short, ha="center", va="center", fontsize=5.7, fontweight="bold", color=color, clip_on=False)
    for row_idx, row in enumerate(quartet_rows):
        self_ax.text(-0.15, row_idx + 0.5, row["x_chr"].replace("chr", "Q"), ha="right", va="center", fontsize=4.9, color=COLORS["muted"])
        for col, (key, _, color) in enumerate(quartet_cols):
            self_ax.add_patch(Rectangle((col + 0.08, row_idx + 0.10), 0.84, 0.8, facecolor=color, edgecolor="white", lw=0.6))
            self_ax.text(col + 0.5, row_idx + 0.5, row[key].replace("chr", ""), ha="center", va="center", fontsize=5.0, color="white", fontweight="bold")
    self_ax.set_xticks([])
    self_ax.set_yticks([])
    for spine in self_ax.spines.values():
        spine.set_color("#AAB4BC")
        spine.set_linewidth(0.6)
    self_ax.axis("off")
    label(ax, 34.25, 5.4, "Evidence-scored assignment: eight haplotype quartets", 7.4, "bold")
    arrow(ax, 68.7, 31.7, 34.2, 26.1, color=COLORS["gold"], lw=0.85, style="-|>", zorder=2)

    # Inset 3: cross-species alignment matrix from the actual PAF summary.
    ref_ax = add_inset(fig, ax, (46.2, 7.0, 16.5, 18.2))
    ref_mat = reference_alignment_matrix()
    ref_ax.imshow(ref_mat, cmap="YlGnBu", aspect="auto", interpolation="nearest")
    ref_ax.set_xticks([])
    ref_ax.set_yticks([])
    for spine in ref_ax.spines.values():
        spine.set_color("#AAB4BC")
        spine.set_linewidth(0.6)
    label(ax, 54.45, 5.4, "Reference synteny: chromosome correspondence", 7.4, "bold")
    arrow(ax, 77.4, 31.7, 54.4, 26.1, color=COLORS["green"], lw=0.85, style="-|>", zorder=2)

    # Inset 4: real SyRI/plotsr panel.
    syri_ax = add_inset(fig, ax, (67.5, 7.0, 14.0, 18.2))
    syri_img = Image.open(ROOT / "final" / "qc" / "syri_plotsr" / "ppa_vs_pzi_copy1_chr16_clean.png").convert("RGB")
    syri_ax.imshow(syri_img)
    syri_ax.set_xlim(80, syri_img.width - 80)
    syri_ax.set_ylim(syri_img.height - 70, 120)
    syri_ax.axis("off")
    label(ax, 74.5, 5.4, "SyRI / plotsr: macro-collinearity", 7.4, "bold")

    # QC summary uses current release tables.
    summary, validation, telomere_ends, hap_bp = summary_metrics()
    qc_ax = add_inset(fig, ax, (84.2, 7.0, 13.4, 18.2))
    qc_ax.axis("off")
    qc_ax.set_xlim(0, 1)
    qc_ax.set_ylim(0, 1)
    qc_ax.add_patch(FancyBboxPatch((0.02, 0.02), 0.96, 0.96, boxstyle="round,pad=0.02,rounding_size=0.035", facecolor="white", edgecolor="#B9C8BC", lw=0.9))
    qc_ax.text(0.08, 0.91, "Release checks", fontsize=8.0, fontweight="bold", color=COLORS["green"], va="top")
    metric_lines = [
        ("Chromosomes", validation["rename_map_chromosomes"]["value"] + " / 32"),
        ("Phased sets", "4 x 8 chr"),
        ("Quartets", validation["overview_quartets"]["value"] + " / 8"),
        ("Chr32 span", f"{hap_bp / 1e9:.3f} Gb"),
        ("Assembly N50", f"{int(summary['n50_bp']) / 1e6:.1f} Mb"),
        ("Telomere ends", f"{telomere_ends} / 64"),
    ]
    for i, (name, value) in enumerate(metric_lines):
        yy = 0.78 - i * 0.115
        qc_ax.text(0.08, yy, name, fontsize=6.8, color=COLORS["muted"], va="center")
        qc_ax.text(0.92, yy, value, fontsize=7.0, fontweight="bold", color=COLORS["ink"], ha="right", va="center")
        qc_ax.plot([0.08, 0.92], [yy - 0.052, yy - 0.052], color="#E4E9E5", lw=0.5)
    qc_ax.text(0.08, 0.075, "Hi-C | self alignment | synteny | SyRI | TIDK | QUAST", fontsize=5.5, color=COLORS["green"], va="bottom")
    label(ax, 90.9, 5.4, "Integrated assembly and phasing QC", 7.4, "bold")
    arrow(ax, 92.5, 31.7, 90.9, 26.1, color=COLORS["green"], lw=0.85, style="-|>", zorder=2)

    # Method note: preserve the biological interpretation of copy groups.
    label(
        ax,
        2.0,
        0.75,
        "Interpretation: hifiasm hap1/hap2 and the two 16-chromosome copy groups are phasing/provenance evidence; A1/A2/B1/B2 are four stable phased output sets, not parent-assigned labels.",
        6.6,
        color=COLORS["muted"],
        ha="left",
    )

    fig.savefig(OUT.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(OUT.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(OUT.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(OUT.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    draw()
