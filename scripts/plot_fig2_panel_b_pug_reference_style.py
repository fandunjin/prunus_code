#!/usr/bin/env python3
"""Reference-style Figure 2b: species tree with PUG branch bars."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable


mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 7,
    "axes.linewidth": 0.75,
})

TAXA = ["PZI", "PSP", "PPA"]
COLORS = {"PZI": "#B64342", "PSP": "#2C64A7", "PPA": "#777777"}
NODES = ["N1", "N2", "N3", "N4", "N5"]


def read_tsv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, fields, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def tree_coordinates():
    y = {"GTR": 0.0, "PPA": 1.0, "PSP": 2.0, "PZI": 3.0, "PCAM": 4.0, "PPE": 5.0}
    x = {"N5": 0.0, "N4": 1.0, "N2": 2.0, "N3": 2.0, "N1": 3.0}
    x.update({taxon: 4.0 for taxon in ["GTR", "PPA", "PSP", "PZI", "PCAM", "PPE"]})
    y["N1"] = (y["PSP"] + y["PZI"]) / 2
    y["N2"] = (y["PPA"] + y["N1"]) / 2
    y["N3"] = (y["PCAM"] + y["PPE"]) / 2
    y["N4"] = (y["N2"] + y["N3"]) / 2
    y["N5"] = (y["GTR"] + y["N4"]) / 2
    return x, y


def connect(ax, x, y, parent, children):
    ax.plot([x[parent], x[parent]], [min(y[c] for c in children), max(y[c] for c in children)], color="black", lw=1.15)
    for child in children:
        ax.plot([x[parent], x[child]], [y[child], y[child]], color="black", lw=1.15)


def draw_branch_bar(ax, x, y, node, tax_counts, pooled_count, pooled_denominator, cmap, norm):
    if not any(tax_counts.get((taxon, node), 0) for taxon in TAXA):
        return
    bar_x = x[node] - 0.86
    bar_y = y[node] - 0.34
    cell_w = 0.24
    cell_h = 0.14
    for index, taxon in enumerate(TAXA):
        value = tax_counts.get((taxon, node), 0)
        face = cmap(norm(value)) if value else "white"
        ax.add_patch(Rectangle(
            (bar_x + index * cell_w, bar_y), cell_w - 0.015, cell_h,
            facecolor=face, edgecolor="black", linewidth=0.45,
            hatch="//" if not value else None,
        ))
    # The reference figure uses a node-specific mapped/supportable ratio,
    # rather than one global selected-family denominator.
    ax.text(
        bar_x + 1.5 * cell_w - 0.02,
        bar_y + cell_h + 0.07,
        f"{pooled_count:,} / {pooled_denominator:,}",
        ha="center", va="bottom", fontsize=6.1,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pug", required=True, type=Path)
    parser.add_argument("--selected-families", required=True, type=Path)
    parser.add_argument(
        "--mapped-rows", required=True, type=Path,
        help="Raw PUG Paralog_Pairs_Nodes_Bootstraps rows before the support cutoff; used for node-specific denominators.",
    )
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()

    pug_rows = read_tsv(args.pug)
    selected_rows = read_tsv(args.selected_families)
    mapped_rows = read_tsv(args.mapped_rows)
    pooled = {(row["event"], row["node"]): int(row["gene_tree_count"]) for row in pug_rows}
    pooled_denominators = {}
    for row in mapped_rows:
        node = row.get("node", row.get("Node", ""))
        if node:
            pooled_denominators[node] = pooled_denominators.get(node, 0) + 1
    tax_counts = {}
    tax_denominators = {taxon: {} for taxon in TAXA}
    for event, taxon in [("pzi_ab", "PZI"), ("psp_ab", "PSP"), ("ppa_ab", "PPA")]:
        for row in pug_rows:
            if row["event"] == event:
                tax_counts[(taxon, row["node"])] = int(row["gene_tree_count"])
        tax_rows_path = args.pug.parent / event / f"fig2b_{event}_Paralog_Pairs_Nodes_Bootstraps.txt"
        if tax_rows_path.exists():
            for row in read_tsv(tax_rows_path):
                tax_denominators[taxon][row["Node"]] = tax_denominators[taxon].get(row["Node"], 0) + 1

    max_tax_count = max(tax_counts.values(), default=1)
    cmap = mpl.colormaps.get_cmap("RdYlBu_r")
    norm = Normalize(vmin=0, vmax=max_tax_count)

    args.outdir.mkdir(parents=True, exist_ok=True)
    x, y = tree_coordinates()
    fig, ax = plt.subplots(figsize=(7.25, 4.95))
    connect(ax, x, y, "N5", ["GTR", "N4"])
    connect(ax, x, y, "N4", ["N2", "N3"])
    connect(ax, x, y, "N2", ["PPA", "N1"])
    connect(ax, x, y, "N1", ["PSP", "PZI"])
    connect(ax, x, y, "N3", ["PCAM", "PPE"])

    for taxon in ["GTR", "PPA", "PSP", "PZI", "PCAM", "PPE"]:
        ax.text(4.08, y[taxon], taxon, va="center", fontsize=8.0, fontstyle="italic")

    for node in NODES:
        draw_branch_bar(
            ax, x, y, node, tax_counts,
            pooled.get(("shared_candidate", node), 0),
            pooled_denominators.get(node, 0),
            cmap, norm,
        )

    ax.text(-0.13, 1.02, "b", transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")
    legend_x = -0.02
    legend_y = 5.86
    for index, taxon in enumerate(TAXA):
        ax.add_patch(Rectangle((legend_x + index * 0.62, legend_y), 0.18, 0.13, facecolor=COLORS[taxon], edgecolor="black", linewidth=0.4))
        ax.text(legend_x + index * 0.62 + 0.22, legend_y + 0.065, taxon, va="center", fontsize=6.4)

    ax.set_xlim(-1.05, 5.0)
    ax.set_ylim(-0.72, 6.12)
    ax.axis("off")
    cax = fig.add_axes([0.015, 0.30, 0.018, 0.34])
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    colorbar = fig.colorbar(sm, cax=cax)
    colorbar.set_ticks([0, max_tax_count])
    colorbar.ax.tick_params(labelsize=6, width=0.4, length=2)
    colorbar.outline.set_linewidth(0.4)
    fig.subplots_adjust(left=0.04, right=0.99, top=0.96, bottom=0.05)
    prefix = args.outdir / "figure2_panel_b_pug_reference_style"
    fig.savefig(str(prefix) + ".svg", bbox_inches="tight")
    fig.savefig(str(prefix) + ".pdf", bbox_inches="tight")
    fig.savefig(str(prefix) + ".tiff", dpi=600, bbox_inches="tight")
    fig.savefig(str(prefix) + ".png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    source_rows = []
    for node in NODES:
        source_rows.append({
            "node": node,
            "pooled_mapped": pooled.get(("shared_candidate", node), 0),
            "selected_families": len(selected_rows),
            "pooled_node_denominator": pooled_denominators.get(node, 0),
            "pooled_ratio": (
                f"{pooled.get(('shared_candidate', node), 0)}/{pooled_denominators.get(node, 0)}"
                if pooled_denominators.get(node, 0) else ""
            ),
            "PZI_mapped": tax_counts.get(("PZI", node), 0),
            "PSP_mapped": tax_counts.get(("PSP", node), 0),
            "PPA_mapped": tax_counts.get(("PPA", node), 0),
            "PZI_node_denominator": tax_denominators["PZI"].get(node, 0),
            "PSP_node_denominator": tax_denominators["PSP"].get(node, 0),
            "PPA_node_denominator": tax_denominators["PPA"].get(node, 0),
            "support_cutoff": 95,
        })
    write_tsv(args.outdir / "panel_b_reference_style_source_data.tsv", list(source_rows[0]), source_rows)
    print(f"[OK] Wrote reference-style Figure 2b: {args.outdir}")


if __name__ == "__main__":
    main()
