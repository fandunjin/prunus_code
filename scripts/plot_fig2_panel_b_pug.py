#!/usr/bin/env python3
"""Render an auditable Figure 2b PUG + ASTRAL panel in Python."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 7,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 0.75,
    "legend.frameon": False,
})


EVENTS = ["shared_candidate", "pzi_ab", "psp_ab", "ppa_ab"]
EVENT_LABELS = {
    "shared_candidate": "Pooled candidate",
    "pzi_ab": "PZI A-B",
    "psp_ab": "PSP A-B",
    "ppa_ab": "PPA A-B",
}
NODES = ["N1", "N2", "N3", "N4", "N5"]


def read_tsv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, fields, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def draw_tree(ax, node_counts):
    y = {"GTR": 0.0, "PPA": 1.0, "PSP": 2.0, "PZI": 3.0, "PCAM": 4.0, "PPE": 5.0}
    x = {
        "N5": 0.0,
        "N4": 1.0,
        "N2": 2.0,
        "N3": 2.0,
        "N1": 3.0,
        "GTR": 4.0,
        "PPA": 4.0,
        "PSP": 4.0,
        "PZI": 4.0,
        "PCAM": 4.0,
        "PPE": 4.0,
    }
    y["N1"] = (y["PSP"] + y["PZI"]) / 2.0
    y["N2"] = (y["PPA"] + y["N1"]) / 2.0
    y["N3"] = (y["PCAM"] + y["PPE"]) / 2.0
    y["N4"] = (y["N2"] + y["N3"]) / 2.0
    y["N5"] = (y["GTR"] + y["N4"]) / 2.0

    def connect(parent, children):
        ax.plot([x[parent], x[parent]], [min(y[c] for c in children), max(y[c] for c in children)], color="#3B3B3B", lw=1.1)
        for child in children:
            ax.plot([x[parent], x[child]], [y[child], y[child]], color="#3B3B3B", lw=1.1)

    connect("N5", ["GTR", "N4"])
    connect("N4", ["N2", "N3"])
    connect("N2", ["PPA", "N1"])
    connect("N1", ["PSP", "PZI"])
    connect("N3", ["PCAM", "PPE"])

    for label in ["GTR", "PPA", "PSP", "PZI", "PCAM", "PPE"]:
        ax.text(4.08, y[label], label, va="center", fontsize=7.2)
    for node in ["N1", "N2", "N3", "N4", "N5"]:
        count = node_counts.get(("shared_candidate", node), 0)
        color = "#B64342" if node == "N4" else "#6D6D6D"
        ax.scatter([x[node]], [y[node]], s=28 if node == "N4" else 18, color=color, zorder=3, edgecolor="white", linewidth=0.5)
        ax.text(x[node] + 0.09, y[node] + 0.15, f"{node}  {count:,}", color=color, fontsize=6.3, fontweight="bold" if node == "N4" else "normal")
    ax.text(0.0, 5.65, "PUG duplication placement", fontsize=8.3, fontweight="bold", ha="left")
    ax.text(0.0, 5.38, "Counts are distinct family trees with support >=95", fontsize=5.7, color="#555555")
    ax.text(0.0, -0.55, "N4 is the common stem of PPA, PSP, PZI, PCAM and PPE; it is not a PZI/PSP-specific node.", fontsize=5.6, color="#555555")
    ax.set_xlim(-0.2, 5.0)
    ax.set_ylim(-0.8, 5.9)
    ax.axis("off")


def draw_heatmap(ax, node_counts):
    matrix = np.asarray([[node_counts.get((event, node), 0) for node in NODES] for event in EVENTS], dtype=float)
    image = ax.imshow(matrix, cmap="Reds", vmin=0, vmax=max(1, matrix.max()), aspect="auto")
    ax.set_xticks(np.arange(len(NODES)), NODES)
    ax.set_yticks(np.arange(len(EVENTS)), [EVENT_LABELS[event] for event in EVENTS])
    ax.tick_params(axis="both", length=0, labelsize=6.1)
    for i in range(len(EVENTS)):
        for j in range(len(NODES)):
            value = int(matrix[i, j])
            ax.text(j, i, f"{value:,}" if value else "-", ha="center", va="center", fontsize=6.0, color="white" if value > matrix.max() * 0.45 else "#333333")
    ax.set_title("High-support duplicate mappings", loc="left", fontsize=8.3, fontweight="bold", pad=7)
    ax.set_xlabel("PUG species-tree node", fontsize=6.7)
    cbar = ax.figure.colorbar(image, ax=ax, fraction=0.045, pad=0.03)
    cbar.ax.tick_params(labelsize=5.6, length=2)
    cbar.set_label("Family trees", fontsize=6.1)
    return matrix


def get_topology_metrics(topology_rows, support_rows):
    topo = next(row for row in topology_rows if row["contrast"] == "pzi_psp_AB" and row["topology"] == "pzi_A+psp_A|pzi_B+psp_B")
    support = {}
    for row in support_rows:
        taxa = set(filter(None, row.get("displayed_descendant_taxa", "").split(",")))
        if taxa == {"psp_A", "pzi_A"}:
            support["A"] = row
        if taxa == {"psp_B", "pzi_B"}:
            support["B"] = row
    return topo, support


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pug", required=True, type=Path)
    parser.add_argument("--topology", required=True, type=Path)
    parser.add_argument("--support", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    pug_rows = read_tsv(args.pug)
    node_counts = {(row["event"], row["node"]): int(row["gene_tree_count"]) for row in pug_rows}
    topology_rows = read_tsv(args.topology)
    support_rows = read_tsv(args.support)
    topo, support = get_topology_metrics(topology_rows, support_rows)

    fig = plt.figure(figsize=(7.35, 5.1))
    grid = fig.add_gridspec(2, 2, width_ratios=[1.08, 1.0], height_ratios=[1.16, 0.84], wspace=0.44, hspace=0.68)
    ax_tree = fig.add_subplot(grid[0, 0])
    ax_heat = fig.add_subplot(grid[0, 1])
    ax_text = fig.add_subplot(grid[1, :])
    draw_tree(ax_tree, node_counts)
    matrix = draw_heatmap(ax_heat, node_counts)
    ax_text.axis("off")
    ax_text.text(0.0, 0.93, "Independent phylogenetic support for parallel A/B histories", fontsize=8.3, fontweight="bold", va="top")
    ax_text.text(0.0, 0.70, f"A-A | B-B topology: {100 * float(topo['frequency']):.2f}% ({topo['gene_count']}/{topo['decisive_genes']} decisive gene trees; 95% CI {100 * float(topo['gene_bootstrap_ci_low']):.2f}-{100 * float(topo['gene_bootstrap_ci_high']):.2f}%)", fontsize=6.8, color="#333333")
    if set(support) == {"A", "B"}:
        ax_text.text(0.0, 0.49, f"PZI-A + PSP-A: ASTRAL PP {float(support['A']['astral_pp']):.2f}, gCF {float(support['A']['gcf']):.1f}, sCF {float(support['A']['likelihood_scf']):.1f}", fontsize=6.8, color="#333333")
        ax_text.text(0.0, 0.29, f"PZI-B + PSP-B: ASTRAL PP {float(support['B']['astral_pp']):.2f}, gCF {float(support['B']['gcf']):.1f}, sCF {float(support['B']['likelihood_scf']):.1f}", fontsize=6.8, color="#333333")
    ax_text.text(0.0, 0.04, "Interpretation boundary: PUG places retained duplicate nodes on the tested species-tree backbone; it does not identify parents, prove maternal origin, or distinguish allopolyploidy from every alternative history.", fontsize=5.9, color="#555555")
    fig.text(0.015, 0.975, "b", fontsize=11, fontweight="bold", va="top")
    prefix = args.outdir / "figure2_panel_b_pug"
    fig.savefig(str(prefix) + ".svg", bbox_inches="tight")
    fig.savefig(str(prefix) + ".pdf", bbox_inches="tight")
    fig.savefig(str(prefix) + ".tiff", dpi=600, bbox_inches="tight")
    fig.savefig(str(prefix) + ".png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    source_rows = []
    for event in EVENTS:
        for node in NODES:
            source_rows.append({"event": event, "node": node, "gene_tree_count": node_counts.get((event, node), 0), "support_cutoff": 95})
    write_tsv(args.outdir / "panel_b_pug_source_data.tsv", ["event", "node", "gene_tree_count", "support_cutoff"], source_rows)
    write_tsv(args.outdir / "panel_b_topology_source_data.tsv", list(topo.keys()), [topo])
    if support:
        write_tsv(args.outdir / "panel_b_astral_support_source_data.tsv", list(next(iter(support.values())).keys()), list(support.values()))
    (args.outdir / "figure2_panel_b_manifest.tsv").write_text(
        "panel\tquestion\tboundary\n"
        "b\tWhere do high-support A/B duplicate nodes map on the Prunus backbone?\tN4 includes PCAM/PPE and is not a PZI/PSP-specific event; PUG does not establish parentage.\n",
        encoding="ascii",
    )
    print(f"[OK] Wrote Figure 2b PUG panel: {args.outdir}")


if __name__ == "__main__":
    main()
