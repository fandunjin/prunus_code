#!/usr/bin/env python3
"""Quick diagnostic: draw the five-species MADS tree with tips coloured by species/track."""
from __future__ import annotations

import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
TREE = ROOT / "analysis" / "inflorescence_family_pipeline_crossspecies5_20260912" / "trees" / "MADS.treefile"
OUT = ROOT / "analysis" / "inflorescence_family_pipeline_crossspecies5_20260912" / "trees" / "MADS_five_species_preview.png"

COLORS = {
    "pzi": "#0F4D92", "psp": "#B64342", "ppa": "#42949E",
    "pcam": "#4D4D4D", "ppe": "#8C8C8C",
}
PHASE_RE = re.compile(r"^(pzi|psp|ppa|pcam|ppe)_hap(A1|A2|B1|B2|A|B|S)_")


class Node:
    __slots__ = ("children", "length", "label", "x", "y")

    def __init__(self):
        self.children = []
        self.length = 0.0
        self.label = ""
        self.x = 0.0
        self.y = 0.0


re_label = re.compile(r"[^:,()]+")
re_len = re.compile(r":\s*([0-9.eE+-]+)")


def parse(text: str) -> Node:
    s = text.strip().rstrip(";")
    pos = 0

    def sub() -> Node:
        nonlocal pos
        n = Node()
        if pos < len(s) and s[pos] == "(":
            pos += 1
            while True:
                n.children.append(sub())
                if pos < len(s) and s[pos] == ",":
                    pos += 1
                    continue
                if pos < len(s) and s[pos] == ")":
                    pos += 1
                    break
            m = re_label.match(s, pos)
            if m:
                n.label = m.group(0)
                pos = m.end()
        else:
            m = re_label.match(s, pos)
            if m:
                n.label = m.group(0)
                pos = m.end()
        m = re_len.match(s, pos)
        if m:
            n.length = float(m.group(1))
            pos = m.end()
        return n

    return sub()


def layout(n: Node, x=0.0, y0=None):
    if y0 is None:
        y0 = [0.0]
    n.x = x + n.length
    if n.children:
        for c in n.children:
            layout(c, n.x, y0)
        n.y = sum(c.y for c in n.children) / len(n.children)
    else:
        n.y = y0[0]
        y0[0] += 1.0


root = parse(TREE.read_text(encoding="utf-8", errors="ignore"))
layout(root)
leaves = []


def collect(n):
    if not n.children:
        leaves.append(n)
        return
    for c in n.children:
        collect(c)


collect(root)
max_x = max(n.x for n in leaves) or 1.0

fig, ax = plt.subplots(figsize=(6.0, 11.0), dpi=160)


def draw(n):
    for c in n.children:
        ax.plot([n.x / max_x, n.x / max_x], [n.y, c.y], color="#B0B0B0", lw=0.25)
        ax.plot([n.x / max_x, c.x / max_x], [c.y, c.y], color="#B0B0B0", lw=0.25)
        draw(c)


draw(root)
for n in leaves:
    m = PHASE_RE.match(n.label)
    color = COLORS.get(m.group(1), "#CCCCCC") if m else "#CCCCCC"
    marker = "o"
    size = 3.0
    if m and m.group(1) in {"pcam", "ppe"}:
        marker = "D"
        size = 4.0
    ax.scatter([n.x / max_x], [n.y], s=size, marker=marker, color=color, zorder=3)
ax.set_xlim(-0.01, 1.01)
ax.set_ylim(len(leaves) + 1, -1)
ax.axis("off")
ax.text(0, -0.002, "MADS tree, five species: blue=PZI red=PSP teal=PPA; diamonds=PCAM/PPE (single track)",
        transform=ax.transAxes, fontsize=7)
fig.savefig(OUT, bbox_inches="tight")
print("saved", OUT)
