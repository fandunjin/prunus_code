#!/usr/bin/env python3
"""Render Figure 2b-style corrected-haplotype homeolog node mapping in Python."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.patches import Rectangle

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none", "pdf.fonttype": 42, "font.size": 7,
    "axes.linewidth": 0.7,
})

TREE = "(gtr,(rch,((ppa,(psp,pzi)),((pmal,(pav,(pspe,pcam))),((psal,(par,pmu)),(pdu,ppe))))));"
DISPLAY = {
    "gtr": "Gillenia trifoliata", "rch": "Rosa chinensis", "ppa": "Prunus padus",
    "psp": "Prunus spinulosa", "pzi": "Prunus zippeliana", "pmal": "Prunus mahaleb",
    "pav": "Prunus avium", "pspe": "Prunus speciosa", "pcam": "Prunus campanulata",
    "psal": "Prunus salicina", "par": "Prunus armeniaca", "pmu": "Prunus mume",
    "pdu": "Prunus dulcis", "ppe": "Prunus persica",
}
PAIR_ORDER = ["PZI", "PSP", "PPA"]
COLORS = {"PZI": "#427AA1", "PSP": "#5C946E", "PPA": "#C9913D"}


@dataclass
class Node:
    name: str = ""
    children: List["Node"] = field(default_factory=list)
    leaves: set = field(default_factory=set)

    @property
    def is_leaf(self): return not self.children


class Parser:
    def __init__(self, text): self.text, self.i = text.strip(), 0
    def token(self, stops):
        a = self.i
        while self.i < len(self.text) and self.text[self.i] not in stops: self.i += 1
        return self.text[a:self.i].strip()
    def branch(self):
        if self.i < len(self.text) and self.text[self.i] == ":": self.i += 1; self.token(",();")
    def subtree(self):
        if self.text[self.i] == "(":
            self.i += 1; children = [self.subtree()]
            while self.text[self.i] == ",": self.i += 1; children.append(self.subtree())
            self.i += 1; name = self.token(":,();"); self.branch()
            return Node(name, children)
        name = self.token(":,();"); self.branch(); return Node(name)
    def parse(self):
        root = self.subtree()
        def leaves(node):
            node.leaves = {node.name} if node.is_leaf else set().union(*(leaves(c) for c in node.children))
            return node.leaves
        leaves(root); return root


def read_tsv(path):
    with Path(path).open(encoding="utf-8", newline="") as f: return list(csv.DictReader(f, delimiter="\t"))


def labels(root):
    out = {}
    named = {
        frozenset({"ppa", "psp", "pzi"}): "Padus-Laurocerasus_crown",
        frozenset({"psp", "pzi"}): "PSP-PZI_crown",
        frozenset({"pmal", "pav", "pspe", "pcam"}): "Cerasus_crown",
        frozenset({"pav", "pspe", "pcam"}): "PAV-PSPE-PCAM_crown",
        frozenset({"pspe", "pcam"}): "PSPE-PCAM_crown",
        frozenset({"psal", "par", "pmu", "pdu", "ppe"}): "Amygdalus-Prunophora_crown",
        frozenset({"psal", "par", "pmu"}): "PSAL-PAR-PMU_crown",
        frozenset({"par", "pmu"}): "PAR-PMU_crown", frozenset({"pdu", "ppe"}): "PDU-PPE_crown",
    }
    def visit(node, parent=None):
        if node.is_leaf:
            out[id(node)] = f"{node.name.upper()}_tip"
        else:
            tips = frozenset(node.leaves)
            if parent is None: out[id(node)] = "Rosaceae_root"
            elif tips == frozenset(root.leaves - {"gtr"}): out[id(node)] = "Prunus-Rosa_crown"
            elif tips == frozenset(root.leaves - {"gtr", "rch"}): out[id(node)] = "Prunus_crown"
            else: out[id(node)] = named.get(tips, "_".join(sorted(tips)) + "_crown")
            for child in node.children: visit(child, node)
    visit(root); return out


def layout(root):
    x, y, leaves = {}, {}, []
    def visit(node, depth):
        x[id(node)] = depth
        if node.is_leaf:
            y[id(node)] = len(leaves); leaves.append(node); return y[id(node)]
        yy = [visit(child, depth + 1) for child in node.children]
        y[id(node)] = sum(yy) / len(yy); return y[id(node)]
    visit(root, 0); return x, y, leaves


def draw_tree(ax, node, x, y):
    for child in node.children:
        ax.plot([x[id(node)], x[id(child)]], [y[id(child)], y[id(child)]], color="#292929", lw=0.8)
        draw_tree(ax, child, x, y)
    if not node.is_leaf:
        yy = [y[id(child)] for child in node.children]
        ax.plot([x[id(node)], x[id(node)]], [min(yy), max(yy)], color="#292929", lw=0.8)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--summary", required=True, type=Path)
    ap.add_argument("--shared", required=True, type=Path)
    ap.add_argument("--outdir", required=True, type=Path)
    args = ap.parse_args()
    data = read_tsv(args.summary)
    shared = read_tsv(args.shared)
    root = Parser(TREE).parse(); names = labels(root); x, y, leaves = layout(root)
    internals = []
    def collect(node):
        if not node.is_leaf:
            internals.append(node)
            for child in node.children: collect(child)
    collect(root)
    lookup = {(r["taxon"], r["mapped_node"]): r for r in data}
    shared_counts = defaultdict(int)
    for r in shared: shared_counts[r["mapped_node"]] += 1
    def shade(taxon, fraction):
        """A taxon-specific white-to-colour strip, matching the three tiles."""
        base = to_rgb(COLORS[taxon])
        weight = min(1.0, max(0.0, fraction))
        return tuple(1.0 - weight * (1.0 - channel) for channel in base)
    fig, ax = plt.subplots(figsize=(7.20, 5.25))
    draw_tree(ax, root, x, y)
    max_x = max(x.values())
    for leaf in leaves:
        ax.text(max_x + .10, y[id(leaf)], DISPLAY[leaf.name], va="center", fontsize=6.5, fontstyle="italic")
    tile_w, tile_h = .27, .20
    source = []
    for node in internals:
        name = names[id(node)]
        values = [lookup.get((tax, name)) for tax in PAIR_ORDER]
        if not any(values): continue
        base_x, base_y = x[id(node)] - .45, y[id(node)] - .38
        for i, (taxon, rec) in enumerate(zip(PAIR_ORDER, values)):
            frac = float(rec["high_support_fraction"]) if rec else 0
            ax.add_patch(Rectangle((base_x + i * tile_w, base_y), tile_w - .012, tile_h,
                                   facecolor=shade(taxon, frac) if rec else "white", edgecolor="#222222", lw=.35,
                                   hatch="//" if not rec else None))
            if rec:
                source.append(dict(rec, shared_high_support_families=shared_counts[name]))
        total_high = sum(int(rec["ufboot_ge_cutoff_families"]) for rec in values if rec)
        total_all = sum(int(rec["all_testable_families"]) for rec in values if rec)
        ax.text(base_x + 1.5 * tile_w, base_y + tile_h + .045, f"{total_high}/{total_all}", ha="center", va="bottom", fontsize=5.7,
                bbox={"facecolor":"white", "edgecolor":"none", "pad":.1})
        if name in {"Padus-Laurocerasus_crown", "Prunus_crown", "PSP-PZI_crown"}:
            friendly = {"Padus-Laurocerasus_crown":"Padus-Laurocerasus", "Prunus_crown":"Prunus", "PSP-PZI_crown":"PSP-PZI"}[name]
            ax.text(base_x + .39, base_y - .11, friendly, ha="center", va="top", fontsize=5.4, color="#444444")
    # Directly label the three colour strips rather than add a legend block.
    key_y = -.80
    for i, taxon in enumerate(PAIR_ORDER):
        key_x = i * .76
        ax.add_patch(Rectangle((key_x, key_y), .13, .13, facecolor=COLORS[taxon], edgecolor="none"))
        ax.text(key_x + .17, key_y + .065, f"{taxon} A/B", va="center", fontsize=6.3)
    ax.text(0.0, key_y - .14, "Tile intensity: UFBoot >=95 family fraction", fontsize=5.7, va="top", color="#444444")
    ax.text(-.065, 1.01, "b", transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom")
    ax.set_xlim(-.65, max_x+3.65); ax.set_ylim(key_y-.38, len(leaves)-.2); ax.axis("off")
    fig.subplots_adjust(left=.07, right=.99, top=.985, bottom=.06)
    args.outdir.mkdir(parents=True, exist_ok=True)
    prefix = args.outdir / "figure2b_corrected_haplotype_homeolog_node_mapping"
    for ext, kwargs in (("svg", {}), ("pdf", {}), ("tiff", {"dpi":600}), ("png", {"dpi":300})):
        fig.savefig(str(prefix)+"."+ext, bbox_inches="tight", **kwargs)
    plt.close(fig)
    fields = list(source[0]) if source else ["taxon"]
    with (args.outdir / "figure2b_corrected_haplotype_source_data.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(source)
    print(f"[OK] Wrote corrected-haplotype Figure 2b mapping: {args.outdir}")

if __name__ == "__main__": main()
