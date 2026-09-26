#!/usr/bin/env python3
"""Review figures for Fig. 5 (tasks 1+2):

  Fig5_ancestry_review  - junction ancestry classes (WGDI shared-breakpoint
                          logic) with gene-environment and TE summaries.
  Fig5_knox_review      - KNOX/BELL loci: recurrence across targets and local
                          WGDI collinearity around the four junctions.

Outputs SVG/PDF/PNG/TIFF into figure5/output/ plus a PNG preview in the temp
directory for inspection.
"""
import csv
import io
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

ROOT = Path("${PROJECT_ROOT}")
V4 = ROOT / "analysis" / "breakpoint_environment_20260912"
ANC = V4 / "ancestry"
KNOX = ROOT / "analysis" / "knox_breakpoint_review_20260912"
FIGDIR = ROOT / "0903" / "manuscript_pipeline" / "stage2_write" / "phase4_draft" / "figure5" / "output"
PREVIEW = Path("${HOME}/AppData/Local/Temp/opencode/figs")

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.7,
    "legend.frameon": False,
})

CLASS_ORDER = ["ancestral_diploid_shared", "ancestral_trio_only", "lineage_specific"]
CLASS_LABEL = {
    "ancestral_diploid_shared": "ancestral (shared, diploid-supported)",
    "ancestral_trio_only": "ancestral (racemose-shared)",
    "lineage_specific": "lineage-specific",
}
CLASS_COLOR = {
    "ancestral_diploid_shared": "#0F4D92",
    "ancestral_trio_only": "#7BA7D7",
    "lineage_specific": "#E08214",
}


def read_tsv(path):
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def save(fig, name):
    FIGDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGDIR / f"{name}.svg", bbox_inches="tight")
    fig.savefig(FIGDIR / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(FIGDIR / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGDIR / f"{name}.tiff", dpi=600, bbox_inches="tight")
    fig.savefig(PREVIEW / f"{name}_preview.png", dpi=180, bbox_inches="tight")
    print("saved", name)


def fig_ancestry():
    counts_rows = read_tsv(ANC / "candidate_ancestry_classes.tsv")
    enrich = read_tsv(ANC / "breakpoint_enrichment_by_ancestry.tsv")
    te = read_tsv(ANC / "breakpoint_te_by_ancestry.tsv")

    short = {"ancestral_diploid_shared": "diploid-shared",
             "ancestral_trio_only": "trio-only",
             "lineage_specific": "lineage-specific"}

    fig = plt.figure(figsize=(7.4, 3.3), dpi=300)
    gs = fig.add_gridspec(1, 3, width_ratios=[0.95, 1.25, 1.25], wspace=0.85,
                          left=0.06, right=0.985, top=0.85, bottom=0.17)

    # ---- a: counts by class and track
    ax = fig.add_subplot(gs[0, 0])
    counts = defaultdict(lambda: defaultdict(int))
    for r in counts_rows:
        counts[r["track"]][r["ancestry_class"]] += 1
    tracks = ["A", "B"]
    y = np.arange(len(tracks))
    left = np.zeros(len(tracks))
    for cls in CLASS_ORDER:
        vals = np.array([counts[t][cls] for t in tracks], dtype=float)
        ax.barh(y, vals, left=left, color=CLASS_COLOR[cls], height=0.55,
                edgecolor="white", linewidth=0.4)
        for yi, (v, l) in enumerate(zip(vals, left)):
            if v:
                ax.text(l + v / 2, yi, f"{int(v)}", ha="center", va="center",
                        fontsize=5.4, color="white" if cls == "ancestral_diploid_shared" else "#222222")
        left += vals
    ax.set_yticks(y)
    ax.set_yticklabels([f"{t}-associated" for t in tracks], fontsize=6.5)
    ax.set_xlabel("junction candidates", fontsize=6.5)
    ax.set_title("a  ancestry classes", loc="left", fontsize=7, fontweight="bold")
    ax.tick_params(labelsize=6)
    ax.legend(handles=[Line2D([], [], marker="s", ls="", color=CLASS_COLOR[c],
                              label=f"{short[c]}  ({counts['A'][c]+counts['B'][c]})", markersize=5)
                       for c in CLASS_ORDER],
              fontsize=5.0, loc="center", bbox_to_anchor=(1.0, 0.5), handletextpad=0.4)

    # ---- b: gene / candidate ratio dots
    ax = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[0, 2], sharey=ax)
    ypos = {}
    yi = 0
    for t in ["A", "B"]:
        for cls in CLASS_ORDER:
            ypos[(t, cls)] = yi
            yi += 1
    for r in enrich:
        t, cls = r["track"], r["ancestry"]
        yv = ypos[(t, cls)]
        for metric, off, color, marker in (("gene", -0.14, "#333333", "o"),
                                           ("candidate", 0.14, "#C03A2B", "s")):
            ratio = float(r[f"{metric}_ratio"] or "nan")
            if np.isfinite(ratio):
                ax.scatter([ratio], [yv + off], s=14, marker=marker, color=color, zorder=3)
    ax.axvline(1, color="#888888", lw=0.7, ls="--", zorder=1)
    ax.set_yticks(list(ypos.values()))
    ax.set_yticklabels([f"{t} · {short[c]}" for (t, c) in ypos.keys()], fontsize=5.8)
    ax.set_xlabel("observed / chromosome-matched null", fontsize=6.5)
    ax.set_xlim(0, 1.45)
    ax.set_ylim(-0.6, len(ypos) - 0.4)
    ax.set_title("b  window content", loc="left", fontsize=7, fontweight="bold")
    ax.tick_params(labelsize=6)
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color="#333333", label="all genes", markersize=4),
                       Line2D([], [], marker="s", ls="", color="#C03A2B", label="floral candidates", markersize=4)],
              fontsize=5.2, loc="upper left", handletextpad=0.3)

    # ---- c: TE ratios
    for r in te:
        if r["track_name"] not in ("copia_coverage", "gypsy_coverage"):
            continue
        ratio = float(r["mean_ratio"] or "nan")
        if not np.isfinite(ratio):
            continue
        yv = ypos[(r["track"], r["ancestry"])]
        off = -0.14 if r["track_name"] == "copia_coverage" else 0.14
        color = "#3F6FB5" if r["track_name"] == "copia_coverage" else "#5B9255"
        ax_c.scatter([ratio], [yv + off], s=14, color=color, zorder=3)
        if r["sign_p"] and float(r["sign_p"]) < 0.1:
            ax_c.text(ratio, yv + off + 0.18, f"p={r['sign_p']}", fontsize=4.4, ha="center", color=color)
    ax_c.axvline(1, color="#888888", lw=0.7, ls="--", zorder=1)
    ax_c.set_xlabel("TE coverage / matched null", fontsize=6.5)
    ax_c.set_xlim(0.5, 1.8)
    ax_c.set_title("c  repeat context", loc="left", fontsize=7, fontweight="bold")
    ax_c.tick_params(labelsize=6)
    ax_c.legend(handles=[Line2D([], [], marker="o", ls="", color="#3F6FB5", label="Copia", markersize=4),
                         Line2D([], [], marker="o", ls="", color="#5B9255", label="Gypsy", markersize=4)],
                fontsize=5.2, loc="lower right", handletextpad=0.3)

    save(fig, "Fig5_ancestry_review")
    plt.close(fig)


def fig_knox():
    hood = read_tsv(KNOX / "knox_locus_neighbourhood.tsv")
    clusters = read_tsv(ROOT / "0903" / "karyotype" / "fig4_breakpoint_candidates" / "junction_candidate_clusters.tsv")
    coll = read_tsv(KNOX / "knox_collinearity_local.tsv")

    cl_by_id = {c["cluster_id"]: c for c in clusters}
    targets = ["pzi_hapA1", "pzi_hapA2", "pzi_hapB1", "pzi_hapB2",
               "psp_hapA1", "psp_hapA2", "psp_hapB1", "psp_hapB2",
               "ppa_hapA", "ppa_hapB", "pcam", "ppe", "gtr", "rch"]
    tshort = {"pzi_hapA1": "PZI-A1", "pzi_hapA2": "PZI-A2", "pzi_hapB1": "PZI-B1", "pzi_hapB2": "PZI-B2",
              "psp_hapA1": "PSP-A1", "psp_hapA2": "PSP-A2", "psp_hapB1": "PSP-B1", "psp_hapB2": "PSP-B2",
              "ppa_hapA": "PPA-A", "ppa_hapB": "PPA-B", "pcam": "PCAM", "ppe": "PPE", "gtr": "GTR", "rch": "RCH"}
    loci = [
        ("C2A\n(AK5, A side)", {"J0009"}, {"C2A"}),
        ("C6B\n(AK5, B side)", {"J0009"}, {"C6B"}),
        ("C1A\n(AK9, A side)", {"J0017", "J0018", "J0045", "J0046"}, {"C1A"}),
        ("C7B\n(AK9, B side)", {"J0017", "J0018"}, {"C7B"}),
    ]
    RACEMOSE_PREFIX = ("pzi", "psp", "ppa")

    dist = {}
    for r in hood:
        if r["radius"] != "200kb":
            continue
        d = float(r["nearest_knox_distance_bp"]) if r["nearest_knox_distance_bp"] else None
        if d is not None:
            key = (r["target_label"], r["chromosome"], r["cluster_id"])
            dist[key] = min(dist.get(key, 1e12), d)

    membership = {c["cluster_id"]: set(c.get("targets", "").split(",")) for c in clusters}

    fig = plt.figure(figsize=(7.4, 7.8), dpi=300)
    gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 0.24, 2.05], hspace=0.30,
                          left=0.16, right=0.985, top=0.95, bottom=0.055)

    # ---- a: recurrence matrix
    ax = fig.add_subplot(gs[0, 0])
    for yi, (lname, clset, chroms) in enumerate(loci):
        for xi, t in enumerate(targets):
            active = {cid for cid in clset if t in membership.get(cid, set())}
            if not active:
                continue
            hit_d = None
            for cid in active:
                for chrom in chroms:
                    k = (t, chrom, cid)
                    if k in dist:
                        hit_d = dist[k] if hit_d is None else min(hit_d, dist[k])
            lineage_only = active <= {"J0045", "J0046"}
            if lineage_only and hit_d is not None:
                ax.scatter([xi], [yi], s=34, marker="D", color="#E08214", zorder=4,
                           edgecolor="white", linewidth=0.4)
                ax.text(xi, yi - 0.42, f"{hit_d/1e3:.0f}", ha="center", va="top", fontsize=4.2, color="#B36A00")
            elif hit_d is not None:
                ax.scatter([xi], [yi], s=42, color="#C03A2B", zorder=4,
                           edgecolor="white", linewidth=0.4)
                ax.text(xi, yi - 0.42, f"{hit_d/1e3:.0f}", ha="center", va="top", fontsize=4.2, color="#444444")
            elif t.startswith(RACEMOSE_PREFIX):
                ax.scatter([xi], [yi], s=12, color="#9E9E9E", zorder=3)
            else:
                ax.scatter([xi], [yi], s=26, facecolor="white", edgecolor="#777777", linewidth=0.7, zorder=3)
    ax.set_xlim(-0.6, len(targets) - 0.4)
    ax.set_ylim(len(loci) - 0.4, -0.85)
    ax.set_xticks(range(len(targets)))
    ax.set_xticklabels([tshort[t] for t in targets], rotation=45, ha="right", fontsize=5.6)
    ax.set_yticks(range(len(loci)))
    ax.set_yticklabels([l[0] for l in loci], fontsize=6)
    for x in (3.5, 7.5, 9.5):
        ax.axvline(x, color="#DDDDDD", lw=0.5, zorder=0)
    ax.set_title("a  KNOX/BELL junction loci across phase sets and outgroups", loc="left",
                 fontsize=6.8, fontweight="bold")
    ax.tick_params(length=2)

    ax_leg = fig.add_subplot(gs[1, 0])
    ax_leg.axis("off")
    ax_leg.legend(handles=[
        Line2D([], [], marker="o", ls="", color="#C03A2B", markersize=4.5,
               label="KNOX/BELL hit within ±200 kb (number = kb to junction)"),
        Line2D([], [], marker="D", ls="", color="#E08214", markersize=4,
               label="PSP-specific transition with KNOX/BELL hit (J0045/J0046)"),
        Line2D([], [], marker="o", ls="", color="#9E9E9E", markersize=3,
               label="signature present, no KNOX/BELL hit in window"),
        Line2D([], [], marker="o", ls="", markerfacecolor="white", markeredgecolor="#777777",
               color="none", label="outgroup signature (no annotation available)"),
    ], fontsize=4.8, loc="center", ncol=2, handletextpad=0.3, columnspacing=1.2, frameon=False)

    # ---- b: collinearity panels
    sel = ["pzi_hapA1_C1A_657_663", "pzi_hapA1_C2A_374_391",
           "pzi_hapB1_C7B_899_914", "pzi_hapB1_C6B_537_560"]
    inner = gs[2, 0].subgridspec(2, 2, hspace=0.42, wspace=0.32)
    by_cand = defaultdict(list)
    for r in coll:
        by_cand[r["candidate_id"]].append(r)
    unit_color = {"GTR03": "#B7B7B7", "GTR04": "#F2B134", "GTR05": "#5B9255",
                  "GTR07": "#C9595F", "GTR08": "#DF8442", "GTR09": "#3A5C93", "C5B": "#8E7CC3"}
    for si, cid in enumerate(sel):
        ax = fig.add_subplot(inner[si // 2, si % 2])
        rs = by_cand.get(cid, [])
        for unit in sorted({r["ref_unit"] for r in rs}):
            pts = [(int(r["target_gene_bp"]), int(r["ref_order"])) for r in rs if r["ref_unit"] == unit]
            ax.scatter([p[0] for p in pts], [p[1] for p in pts], s=3.2,
                       color=unit_color.get(unit, "#999999"), alpha=0.75, label=unit)
        knox = [(int(r["target_gene_bp"]), int(r["ref_order"])) for r in rs if r["knox"] == "yes"]
        if knox:
            ax.scatter([p[0] for p in knox], [p[1] for p in knox], s=64, marker="*",
                       facecolor="#C03A2B", edgecolor="black", linewidth=0.4, zorder=5, label="KNOX/BELL")
        center = int(rs[0]["junction_center_bp"]) if rs else 0
        ax.axvline(center, color="black", ls="--", lw=0.7)
        ax.set_title(cid.replace("pzi_", "").replace("_", "  "), fontsize=6.2)
        ax.set_xlabel("target bp (Mb)", fontsize=6)
        ax.set_ylabel("reference gene order", fontsize=6)
        ax.tick_params(labelsize=5.4)
        ax.xaxis.set_major_formatter(lambda v, _p: f"{v/1e6:.2f}")
        if si == 0:
            ax.legend(fontsize=4.8, loc="lower right", ncol=2, handletextpad=0.3, columnspacing=0.6)
    save(fig, "Fig5_knox_review")
    plt.close(fig)


def main():
    fig_ancestry()
    fig_knox()
    return 0


if __name__ == "__main__":
    sys.exit(main())
