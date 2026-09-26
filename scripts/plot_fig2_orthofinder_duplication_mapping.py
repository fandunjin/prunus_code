#!/usr/bin/env python3
"""Plot native OrthoFinder duplication mappings on its labelled species tree.

This plot is intentionally an audit of OrthoFinder reconciliation output.  It
does not classify duplications as WGD, identify parents, or estimate a
polyploidization time.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle


mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 7,
})


LEAF_LABELS = {
    "rch": "Rosa chinensis",
    "gtr": "Gillenia trifoliata",
    "pzi_A": "Prunus zippeliana A",
    "pzi_B": "Prunus zippeliana B",
    "psp_A": "Prunus spinulosa A",
    "psp_B": "Prunus spinulosa B",
    "ppa_A": "Prunus padus A",
    "ppa_B": "Prunus padus B",
    "pmal": "Prunus mahaleb",
    "pav": "Prunus avium",
    "pspe": "Prunus speciosa",
    "pcam": "Prunus campanulata",
    "psal": "Prunus salicina",
    "par": "Prunus armeniaca",
    "pmu": "Prunus mume",
    "pdu": "Prunus dulcis",
    "ppe": "Prunus persica",
}


@dataclass
class Node:
    name: str = ""
    children: list["Node"] = field(default_factory=list)

    @property
    def is_leaf(self) -> bool:
        return not self.children


class NewickParser:
    def __init__(self, text: str):
        self.text = text.strip()
        self.index = 0

    def parse(self) -> Node:
        node = self._subtree()
        self._skip_space()
        if self.index < len(self.text) and self.text[self.index] == ";":
            self.index += 1
        return node

    def _skip_space(self) -> None:
        while self.index < len(self.text) and self.text[self.index].isspace():
            self.index += 1

    def _token(self, stop: str) -> str:
        start = self.index
        while self.index < len(self.text) and self.text[self.index] not in stop:
            self.index += 1
        return self.text[start:self.index].strip()

    def _skip_branch_length(self) -> None:
        if self.index < len(self.text) and self.text[self.index] == ":":
            self.index += 1
            self._token(",();")

    def _subtree(self) -> Node:
        self._skip_space()
        if self.text[self.index] == "(":
            self.index += 1
            children = [self._subtree()]
            self._skip_space()
            while self.index < len(self.text) and self.text[self.index] == ",":
                self.index += 1
                children.append(self._subtree())
                self._skip_space()
            if self.index >= len(self.text) or self.text[self.index] != ")":
                raise ValueError("Malformed Newick: missing ')'")
            self.index += 1
            self._skip_space()
            node = Node(name=self._token(":,();"), children=children)
            self._skip_branch_length()
            return node
        node = Node(name=self._token(":,();"))
        if not node.name:
            raise ValueError("Malformed Newick: empty leaf")
        self._skip_branch_length()
        return node


def read_counts(path: Path) -> dict[str, dict[str, int]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    expected = {"Species Tree Node", "Duplications (all)", "Duplications (50% support)"}
    if not rows or not expected.issubset(rows[0]):
        raise ValueError(f"Unexpected OrthoFinder duplication table: {path}")
    return {
        row["Species Tree Node"]: {
            "all": int(row["Duplications (all)"]),
            "support50": int(row["Duplications (50% support)"]),
        }
        for row in rows
    }


def layout(root: Node):
    x, y, leaves, internal = {}, {}, [], []

    def visit(node: Node, depth: int) -> float:
        x[id(node)] = depth
        if node.is_leaf:
            leaves.append(node)
            y[id(node)] = float(len(leaves) - 1)
            return y[id(node)]
        internal.append(node)
        child_y = [visit(child, depth + 1) for child in node.children]
        y[id(node)] = sum(child_y) / len(child_y)
        return y[id(node)]

    visit(root, 0)
    return x, y, leaves, internal


def draw_tree(ax, node: Node, x, y) -> None:
    for child in node.children:
        ax.plot([x[id(node)], x[id(child)]], [y[id(child)], y[id(child)]], color="#2b2b2b", lw=0.8)
    if not node.is_leaf:
        ax.plot(
            [x[id(node)], x[id(node)]],
            [min(y[id(child)] for child in node.children), max(y[id(child)] for child in node.children)],
            color="#2b2b2b",
            lw=0.8,
        )
        for child in node.children:
            draw_tree(ax, child, x, y)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree", required=True, type=Path)
    parser.add_argument("--counts", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()

    root = NewickParser(args.tree.read_text(encoding="utf-8")).parse()
    counts = read_counts(args.counts)
    x, y, leaves, internal = layout(root)
    internal_names = {node.name for node in internal}
    missing = internal_names - set(counts)
    if missing:
        raise SystemExit(f"Missing internal-node counts: {', '.join(sorted(missing))}")

    internal_support = [counts[node.name]["support50"] for node in internal]
    positive = [value for value in internal_support if value > 0]
    norm = LogNorm(vmin=max(1, min(positive)), vmax=max(positive))
    cmap = mpl.colormaps.get_cmap("RdYlBu_r")

    fig, ax = plt.subplots(figsize=(7.3, 7.1))
    draw_tree(ax, root, x, y)
    right = max(x.values())
    for leaf in leaves:
        label = LEAF_LABELS.get(leaf.name, leaf.name)
        ax.text(right + 0.12, y[id(leaf)], label, va="center", fontsize=7.0, fontstyle="italic")

    # The coloured cell gives high-support duplications; the printed pair keeps
    # the complete native OrthoFinder count auditable without a second colour scale.
    cell_w, cell_h = 0.34, 0.20
    source_rows = []
    for node in internal:
        name = node.name
        record = counts[name]
        support = record["support50"]
        face = cmap(norm(max(1, support))) if support else "white"
        tile_x = x[id(node)] - 0.55
        tile_y = y[id(node)] - 0.10
        ax.add_patch(Rectangle((tile_x, tile_y), cell_w, cell_h, facecolor=face, edgecolor="#222222", linewidth=0.35))
        ax.text(tile_x + cell_w / 2, tile_y + cell_h + 0.055, f"{record['all']:,} / {support:,}", ha="center", va="bottom", fontsize=5.4)
        source_rows.append({
            "species_tree_node": name,
            "duplications_all": record["all"],
            "duplications_support_ge_50_percent": support,
            "displayed_label": f"{record['all']} / {support}",
        })

    ax.text(-0.06, 1.01, "b", transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")
    ax.text(0.0, 1.01, "OrthoFinder reconciliation: duplication mappings", transform=ax.transAxes, fontsize=7.6, va="bottom")
    ax.text(0.0, -0.08, "Labels: all duplications / duplications with support >=50%; colour encodes the latter.", transform=ax.transAxes, fontsize=5.8, color="#444444")
    ax.text(0.0, -0.115, "A/B are subgenome labels; this plot is a native duplication summary, not a WGD-event count.", transform=ax.transAxes, fontsize=5.8, color="#444444")
    ax.set_xlim(-0.72, right + 4.1)
    ax.set_ylim(-0.45, len(leaves) - 0.30)
    ax.axis("off")

    cax = fig.add_axes([0.028, 0.34, 0.015, 0.25])
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    colorbar = fig.colorbar(sm, cax=cax)
    colorbar.set_ticks([norm.vmin, norm.vmax])
    colorbar.ax.tick_params(labelsize=5.8, width=0.4, length=2)
    colorbar.outline.set_linewidth(0.4)
    colorbar.set_label("Duplications\n(support >=50%)", fontsize=5.8)

    args.outdir.mkdir(parents=True, exist_ok=True)
    prefix = args.outdir / "figure2b_orthofinder_native_duplication_mapping"
    fig.savefig(str(prefix) + ".svg", bbox_inches="tight")
    fig.savefig(str(prefix) + ".pdf", bbox_inches="tight")
    fig.savefig(str(prefix) + ".tiff", dpi=600, bbox_inches="tight")
    fig.savefig(str(prefix) + ".png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    with (args.outdir / "figure2b_orthofinder_native_duplication_source_data.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(source_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(source_rows)
    (args.outdir / "method_boundary.txt").write_text(
        "Input: native OrthoFinder Gene_Duplication_Events reconciliation summary.\n"
        "Tree: OrthoFinder rooted species tree with A/B subgenome leaves retained.\n"
        "Display: all duplications / duplications with OrthoFinder support >=50%.\n"
        "Boundary: counts include all gene-family duplications and are not WGD-event counts.\n",
        encoding="ascii",
    )
    print(f"[OK] Wrote native OrthoFinder duplication mapping: {args.outdir}")


if __name__ == "__main__":
    main()
