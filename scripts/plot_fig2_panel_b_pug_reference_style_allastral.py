#!/usr/bin/env python3
"""Plot the reference-style PUG panel using the complete ASTRAL taxon set."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle


mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 7,
    "axes.linewidth": 0.65,
})

TAXA = ["PZI", "PSP", "PPA"]
EVENTS = {"PZI": "pzi_ab", "PSP": "psp_ab", "PPA": "ppa_ab"}
LEAF_LABELS = {
    "ppa": "Prunus padus",
    "psp": "Prunus spinulosa",
    "pzi": "Prunus zippeliana",
    "pmal": "Prunus mahaleb",
    "pav": "Prunus avium",
    "pspe": "Prunus speciosa",
    "pcam": "Prunus campanulata",
    "psal": "Prunus salicina",
    "par": "Prunus armeniaca",
    "pmu": "Prunus mume",
    "pdu": "Prunus dulcis",
    "ppe": "Prunus persica",
    "gtr": "Gillenia trifoliata",
    "rch": "Rosa chinensis",
}


@dataclass
class TreeNode:
    name: str = ""
    children: list["TreeNode"] = field(default_factory=list)

    @property
    def is_leaf(self):
        return not self.children


class NewickParser:
    def __init__(self, text: str):
        self.text = text.strip()
        self.i = 0

    def parse(self):
        node = self._subtree()
        self._skip_space()
        if self.i < len(self.text) and self.text[self.i] == ";":
            self.i += 1
        return node

    def _skip_space(self):
        while self.i < len(self.text) and self.text[self.i].isspace():
            self.i += 1

    def _token(self, stops):
        start = self.i
        while self.i < len(self.text) and self.text[self.i] not in stops:
            self.i += 1
        return self.text[start:self.i].strip()

    def _branch_length(self):
        if self.i < len(self.text) and self.text[self.i] == ":":
            self.i += 1
            self._token(",();")

    def _subtree(self):
        self._skip_space()
        if self.text[self.i] == "(":
            self.i += 1
            children = [self._subtree()]
            self._skip_space()
            while self.i < len(self.text) and self.text[self.i] == ",":
                self.i += 1
                children.append(self._subtree())
                self._skip_space()
            if self.i >= len(self.text) or self.text[self.i] != ")":
                raise ValueError("Malformed Newick: missing closing parenthesis")
            self.i += 1
            self._skip_space()
            name = self._token(":,();")
            self._branch_length()
            return TreeNode(name=name, children=children)
        name = self._token(":,();")
        self._branch_length()
        if not name:
            raise ValueError("Malformed Newick: empty leaf")
        return TreeNode(name=name)


def read_tsv(path: Path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows):
    if not rows:
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def layout(root: TreeNode):
    x = {}
    y = {}
    leaves = []

    def visit(node, depth):
        x[id(node)] = depth
        if node.is_leaf:
            leaves.append(node)
            y[id(node)] = float(len(leaves) - 1)
            return y[id(node)]
        child_y = [visit(child, depth + 1) for child in node.children]
        y[id(node)] = sum(child_y) / len(child_y)
        return y[id(node)]

    visit(root, 0)
    return x, y, leaves


def all_internal(root):
    nodes = []

    def visit(node):
        if not node.is_leaf:
            nodes.append(node)
            for child in node.children:
                visit(child)

    visit(root)
    return nodes


def draw_tree(ax, root, x, y):
    def visit(node):
        for child in node.children:
            ax.plot([x[id(node)], x[id(child)]], [y[id(child)], y[id(child)]], color="black", lw=0.85)
        if not node.is_leaf:
            ax.plot(
                [x[id(node)], x[id(node)]],
                [min(y[id(child)] for child in node.children), max(y[id(child)] for child in node.children)],
                color="black", lw=0.85,
            )
            for child in node.children:
                visit(child)

    visit(root)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pug", required=True, type=Path)
    parser.add_argument("--mapped-rows", required=True, type=Path)
    parser.add_argument("--species-tree", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()

    pug_rows = read_tsv(args.pug)
    raw_rows = read_tsv(args.mapped_rows)
    root = NewickParser(args.species_tree.read_text(encoding="utf-8")).parse()
    x, y, leaves = layout(root)
    internals = all_internal(root)
    labeled_nodes = {node.name for node in internals if node.name}
    if not labeled_nodes:
        raise SystemExit(
            "Species tree has no internal PUG node labels; pass the PUG "
            "*_Labeled_Species_Tree.tre file, not the unlabelled input tree."
        )

    pooled = {(row["event"], row["node"]): int(row["gene_tree_count"]) for row in pug_rows}
    denominators = {}
    for row in raw_rows:
        node = row.get("Node", row.get("node", ""))
        if node:
            denominators[node] = denominators.get(node, 0) + 1
    tax_counts = {}
    tax_denominators = {taxon: {} for taxon in TAXA}
    for taxon, event in EVENTS.items():
        for row in pug_rows:
            if row["event"] == event:
                tax_counts[(taxon, row["node"])] = int(row["gene_tree_count"])
        raw_path = args.mapped_rows.parent.parent / event / f"fig2b_{event}_Paralog_Pairs_Nodes_Bootstraps.txt"
        if raw_path.exists():
            for row in read_tsv(raw_path):
                node = row.get("Node", row.get("node", ""))
                if node:
                    tax_denominators[taxon][node] = tax_denominators[taxon].get(node, 0) + 1

    max_count = max(tax_counts.values(), default=1)
    cmap = mpl.colormaps.get_cmap("RdYlBu_r")
    norm = Normalize(vmin=0, vmax=max_count)
    # Keep the standalone panel within a conventional 183-mm double-column width.
    fig, ax = plt.subplots(figsize=(7.2, 6.8))
    draw_tree(ax, root, x, y)
    max_x = max(x.values())
    for leaf in leaves:
        label = LEAF_LABELS.get(leaf.name, leaf.name)
        ax.text(max_x + 0.10, y[id(leaf)], label, va="center", fontsize=7.4,
                fontstyle="italic" if label.startswith(("Prunus", "Gillenia", "Rosa")) else "normal")

    cell_w, cell_h = 0.29, 0.22
    for node in internals:
        node_name = node.name
        denominator = denominators.get(node_name, 0)
        numerator = pooled.get(("shared_candidate", node_name), 0)
        if not denominator:
            continue
        bar_x = x[id(node)] - 0.52
        bar_y = y[id(node)] - 0.38
        for idx, taxon in enumerate(TAXA):
            count = tax_counts.get((taxon, node_name), 0)
            face = cmap(norm(count)) if count else "white"
            ax.add_patch(Rectangle(
                (bar_x + idx * cell_w, bar_y), cell_w - 0.018, cell_h,
                facecolor=face, edgecolor="black", linewidth=0.35,
                hatch="//" if not count else None,
            ))
        ax.text(
            bar_x + 1.5 * cell_w - 0.02,
            bar_y + cell_h + 0.06,
            f"{numerator:,} / {denominator:,}",
            ha="center", va="bottom", fontsize=6.2, zorder=5,
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.15},
        )

    # Species order under the cells follows the reference-figure convention.
    legend_y = -0.95
    for idx, taxon in enumerate(TAXA):
        ax.add_patch(Rectangle((idx * cell_w, legend_y), cell_w - 0.018, cell_h,
                               facecolor="white", edgecolor="black", linewidth=0.35))
        ax.text(idx * cell_w + (cell_w - 0.018) / 2, legend_y - 0.08, taxon,
                ha="center", va="top", fontsize=5.8)

    cax = fig.add_axes([0.025, 0.33, 0.016, 0.30])
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cb = fig.colorbar(sm, cax=cax)
    cb.set_ticks([0, max_count])
    cb.ax.tick_params(labelsize=6, width=0.4, length=2)
    cb.outline.set_linewidth(0.4)

    ax.text(-0.06, 1.01, "b", transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")
    ax.set_xlim(-0.65, max_x + 4.75)
    ax.set_ylim(legend_y - 0.55, len(leaves) - 0.15)
    ax.axis("off")
    fig.subplots_adjust(left=0.08, right=0.98, top=0.98, bottom=0.05)
    args.outdir.mkdir(parents=True, exist_ok=True)
    prefix = args.outdir / "figure2_panel_b_pug_reference_style_allastral"
    fig.savefig(str(prefix) + ".svg", bbox_inches="tight")
    fig.savefig(str(prefix) + ".pdf", bbox_inches="tight")
    fig.savefig(str(prefix) + ".tiff", dpi=600, bbox_inches="tight")
    fig.savefig(str(prefix) + ".png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    source = []
    for node in internals:
        name = node.name
        if name and denominators.get(name, 0):
            source.append({
                "node": name,
                "pooled_high_support": pooled.get(("shared_candidate", name), 0),
                "pooled_all_mapped": denominators[name],
                "pooled_ratio": f"{pooled.get(('shared_candidate', name), 0)}/{denominators[name]}",
                "PZI_high_support": tax_counts.get(("PZI", name), 0),
                "PSP_high_support": tax_counts.get(("PSP", name), 0),
                "PPA_high_support": tax_counts.get(("PPA", name), 0),
                "PZI_all_mapped": tax_denominators["PZI"].get(name, 0),
                "PSP_all_mapped": tax_denominators["PSP"].get(name, 0),
                "PPA_all_mapped": tax_denominators["PPA"].get(name, 0),
                "support_cutoff": 95,
            })
    write_tsv(args.outdir / "panel_b_allastral_source_data.tsv", source)
    print(f"[OK] Wrote all-ASTRAL-taxa reference-style panel b: {args.outdir}")


if __name__ == "__main__":
    main()
