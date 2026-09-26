#!/usr/bin/env python3
"""Render a compact copyAB OrthoFinder duplication-node mapping."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"], "font.size": 7, "svg.fonttype": "none", "pdf.fonttype": 42})


def parse_newick(s: str):
    i = 0
    def subtree():
        nonlocal i
        if s[i] == "(":
            i += 1; kids = [subtree()]
            while s[i] == ",": i += 1; kids.append(subtree())
            if s[i] != ")": raise ValueError("bad Newick")
            i += 1
            start = i
            while i < len(s) and s[i] not in ":,();": i += 1
            name = s[start:i].strip()
            if i < len(s) and s[i] == ":":
                i += 1
                while i < len(s) and s[i] not in ",();": i += 1
            return {"name": name, "kids": kids}
        start = i
        while i < len(s) and s[i] not in ",();:": i += 1
        name = s[start:i].strip()
        if i < len(s) and s[i] == ":":
            i += 1
            while i < len(s) and s[i] not in ",();": i += 1
        return {"name": name, "kids": []}
    return subtree()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", type=Path, required=True)
    ap.add_argument("--counts", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()
    root = parse_newick(args.tree.read_text(encoding="utf-8").strip().rstrip(";"))
    data = {}
    with args.counts.open(encoding="utf-8", newline="") as h:
        for r in csv.DictReader(h, delimiter="\t"):
            data[r["species_tree_node"]] = r

    x, y, leaves, internal = {}, {}, [], []
    def layout(n, depth=0):
        x[id(n)] = depth
        if not n["kids"]:
            leaves.append(n); y[id(n)] = len(leaves) - 1; return y[id(n)]
        internal.append(n); vals = [layout(k, depth + 1) for k in n["kids"]]; y[id(n)] = sum(vals) / len(vals); return y[id(n)]
    layout(root)
    fig, (ax, bar) = plt.subplots(1, 2, figsize=(8.2, max(4.8, len(leaves) * 0.24)), gridspec_kw={"width_ratios": [1.9, 1.0]})
    def draw(n):
        for k in n["kids"]:
            ax.plot([x[id(n)], x[id(k)]], [y[id(k)], y[id(k)]], color="#333", lw=.7); draw(k)
        if n["kids"]: ax.plot([x[id(n)], x[id(n)]], [min(y[id(k)] for k in n["kids"]), max(y[id(k)] for k in n["kids"])], color="#333", lw=.7)
    draw(root)
    xmax = max(x.values())
    for n in leaves:
        label = n["name"].replace("_copyAB", " copyAB")
        ax.text(xmax + .12, y[id(n)], label, va="center", fontsize=6.5)
    vals = [int(data[n["name"]]["duplication_events_support_ge_cutoff"]) for n in internal if n["name"] in data]
    vmax = max(vals) if vals else 1
    norm = Normalize(0, vmax)
    cmap = mpl.colormaps["viridis"]
    for n in internal:
        r = data.get(n["name"], {"duplication_events_all_rows": "0", "duplication_events_support_ge_cutoff": "0"})
        all_n = int(r["duplication_events_all_rows"]); sup = int(r["duplication_events_support_ge_cutoff"])
        ax.scatter(x[id(n)] - .12, y[id(n)], s=18 + 35 * norm(sup), c=[sup], cmap=cmap, norm=norm, edgecolor="white", linewidth=.3, zorder=3)
        ax.text(x[id(n)] - .12, y[id(n)] + .18, f"{all_n}/{sup}", ha="center", va="bottom", fontsize=5)
    ax.axis("off"); ax.set_xlim(-.35, xmax + 2.6); ax.set_ylim(-.7, len(leaves) - .3)
    bar.axis("off")
    sm = ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
    cb = fig.colorbar(sm, ax=bar, fraction=.24, pad=.05); cb.ax.tick_params(labelsize=6); cb.set_label("Supported duplication families", fontsize=6)
    fig.suptitle("CopyAB OrthoFinder duplication mapping", x=.04, ha="left", fontsize=8.5, fontweight="bold")
    args.outdir.mkdir(parents=True, exist_ok=True)
    prefix = args.outdir / "figure_copyAB_duplication_node_mapping"
    fig.savefig(str(prefix) + ".svg", bbox_inches="tight"); fig.savefig(str(prefix) + ".pdf", bbox_inches="tight"); fig.savefig(str(prefix) + ".png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] Wrote copyAB duplication mapping: {prefix}.pdf")


if __name__ == "__main__":
    main()
