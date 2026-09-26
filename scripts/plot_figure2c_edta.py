#!/usr/bin/env python3
"""Figure 2c redrawn from the 0903 per-haplotype EDTA annotations.

Left: stacked TE length (Mb). Right: the same classes as genome fraction.
Display groups PZI-A/PZI-B and PSP-A/PSP-B are the mean of the two A or B
haplotypes; PPA-A/PPA-B, PPE and PCAM are single assemblies.
"""

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams.update({
    "font.size": 7,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
})

CLASSES = ["Copia", "Gypsy", "CACTA", "hAT", "Mutator", "PIF", "Tc1", "Helitron", "Other TEs"]
TE_COLORS = {
    "Copia": "#334E7D",
    "Gypsy": "#2E8A78",
    "CACTA": "#C94C43",
    "hAT": "#D28B2D",
    "Mutator": "#8065A1",
    "PIF": "#B65E8A",
    "Tc1": "#B8A14A",
    "Helitron": "#5F6268",
    "Other TEs": "#A3A3A3",
}
DISPLAY_ORDER = ["PZI-A", "PZI-B", "PSP-A", "PSP-B", "PPA-A", "PPA-B", "PPE", "PCAM"]
GROUP_MEMBERS = {
    "PZI-A": ["pziA1", "pziA2"],
    "PZI-B": ["pziB1", "pziB2"],
    "PSP-A": ["pspA1", "pspA2"],
    "PSP-B": ["pspB1", "pspB2"],
    "PPA-A": ["ppaA"],
    "PPA-B": ["ppaB"],
    "PPE": ["ppe"],
    "PCAM": ["pcam"],
}

BASE = Path("${PROJECT_ROOT}")
IN = BASE / "0903/figures/figure2_abc/figure2c_te_summary_edta.tsv"


def main():
    rows = list(csv.DictReader(open(IN, encoding="utf-8"), delimiter="\t"))
    bp = defaultdict(lambda: defaultdict(float))
    genome = defaultdict(float)
    for r in rows:
        bp[r["target"]][r["class"]] = float(r["bp"])
        genome[r["target"]] = float(r["genome_bp"])

    disp_bp = {}
    disp_genome = {}
    for d in DISPLAY_ORDER:
        members = GROUP_MEMBERS[d]
        disp_bp[d] = {c: float(np.mean([bp[m].get(c, 0.0) for m in members])) for c in CLASSES}
        disp_genome[d] = float(np.mean([genome[m] for m in members]))

    fig = plt.figure(figsize=(4.6, 3.5), facecolor="white")
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1], wspace=0.28, left=0.20, right=0.985, top=0.90, bottom=0.30)
    ax_len = fig.add_subplot(gs[0, 0])
    ax_pct = fig.add_subplot(gs[0, 1], sharey=ax_len)

    y = np.arange(len(DISPLAY_ORDER))
    left_len = np.zeros(len(DISPLAY_ORDER))
    left_pct = np.zeros(len(DISPLAY_ORDER))
    for c in CLASSES:
        vals_mb = np.array([disp_bp[d][c] / 1e6 for d in DISPLAY_ORDER])
        vals_frac = np.array([disp_bp[d][c] / disp_genome[d] for d in DISPLAY_ORDER])
        ax_len.barh(y, vals_mb, left=left_len, height=0.66, color=TE_COLORS[c], edgecolor="white", linewidth=0.25, label=c)
        ax_pct.barh(y, vals_frac, left=left_pct, height=0.66, color=TE_COLORS[c], edgecolor="white", linewidth=0.25)
        left_len += vals_mb
        left_pct += vals_frac

    ax_len.set_yticks(y, DISPLAY_ORDER)
    ax_len.invert_yaxis()
    ax_len.set_xlabel("TE length (Mb)", labelpad=4)
    ax_len.set_xlim(0, np.ceil(left_len.max() / 50.0) * 50)
    ax_pct.set_xlim(0, 0.8)
    ax_pct.set_xticks([0, 0.2, 0.4, 0.6, 0.8], ["0%", "20%", "40%", "60%", "80%"])
    ax_pct.set_xlabel("TE / genome", labelpad=4)
    ax_len.tick_params(labelsize=6.5)
    ax_pct.tick_params(labelsize=6.5)
    ax_len.grid(axis="x", color="#E7E7E7", lw=0.5)
    ax_pct.grid(axis="x", color="#E7E7E7", lw=0.5)
    ax_len.set_axisbelow(True)
    ax_pct.set_axisbelow(True)
    ax_len.set_title("TE composition", loc="left", pad=6, fontsize=8, fontweight="bold")

    handles = [plt.Rectangle((0, 0), 1, 1, color=TE_COLORS[c]) for c in CLASSES]
    fig.legend(handles, CLASSES, loc="lower center", ncol=5, fontsize=6, frameon=False,
               bbox_to_anchor=(0.5, 0.0), handlelength=1.0, handleheight=0.9, columnspacing=1.2)
    out = BASE / "0903/figures/figure2_abc/figure2c_te_composition_edta"
    for ext in ("png", "pdf", "svg", "tiff"):
        fig.savefig(f"{out}.{ext}", dpi=600, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    main()
