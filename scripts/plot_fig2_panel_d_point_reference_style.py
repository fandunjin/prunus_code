#!/usr/bin/env python3
"""Plot formal Fig. 2d from a normalized, converged POInT likelihood export.

Candidate anchor states, JCVI anchors, and Dollo-parsimony summaries are
deliberately rejected. They are not substitutes for POInT.
"""

from __future__ import print_function

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["font.size"] = 7
plt.rcParams["axes.linewidth"] = 0.75
plt.rcParams["axes.spines.right"] = False
plt.rcParams["axes.spines.top"] = False
plt.rcParams["legend.frameon"] = False


BRANCHES = ["post_duplication", "prunus_stem", "psp_pzi_stem", "pzi_tip", "psp_tip", "ppa_tip"]
BRANCH_LABELS = {
    "post_duplication": "post-duplication",
    "prunus_stem": "Prunus stem",
    "psp_pzi_stem": "PSP-PZI stem",
    "pzi_tip": "PZI",
    "psp_tip": "PSP",
    "ppa_tip": "PPA",
}
Y_POS = {
    "post_duplication": 0.84,
    "prunus_stem": 0.67,
    "psp_pzi_stem": 0.48,
    "pzi_tip": 0.28,
    "psp_tip": 0.16,
    "ppa_tip": 0.04,
}
COLORS = {"A": "#5B7FCA", "B": "#E28E2C"}


def read_tsv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, fields, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--point-output", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()
    rows = read_tsv(args.point_output)
    required = {"branch_id", "subgenome", "surviving_genes", "losses", "retention_ratio", "posterior", "analysis_status"}
    if not rows or not required.issubset(rows[0]):
        raise SystemExit("ERROR: normalized output lacks required columns: {}".format(", ".join(sorted(required))))
    lowered_path = str(args.point_output).lower()
    forbidden = ("candidate", "alltaxa_fractionation_states", "jcvi", "dollo")
    if any(token in lowered_path for token in forbidden):
        raise SystemExit("ERROR: candidate/JCVI/Dollo data cannot be plotted as formal POInT output")
    statuses = {row["analysis_status"].lower() for row in rows}
    if statuses - {"complete", "converged"}:
        raise SystemExit("ERROR: POInT output is not complete/converged: {}".format(sorted(statuses)))
    for row in rows:
        if row["subgenome"] not in ("A", "B") or row["branch_id"] not in BRANCHES:
            raise SystemExit("ERROR: unexpected POInT branch or subgenome: {}".format(row))
        if float(row["retention_ratio"]) < 0 or float(row["retention_ratio"]) > 1:
            raise SystemExit("ERROR: invalid retention ratio: {}".format(row))
        if float(row["posterior"]) < 0 or float(row["posterior"]) > 1:
            raise SystemExit("ERROR: invalid posterior: {}".format(row))
    by_key = {(row["branch_id"], row["subgenome"]): row for row in rows}
    for branch in BRANCHES:
        for subgenome in ("A", "B"):
            if (branch, subgenome) not in by_key:
                raise SystemExit("ERROR: missing POInT row for {} {}".format(branch, subgenome))

    outdir = Path(args.outdir)
    if outdir.exists() and any(outdir.iterdir()):
        raise SystemExit("ERROR: output directory is not empty: {}".format(outdir))
    outdir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(4.3, 4.2))
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.04, 1.0)
    ax.axis("off")
    # Upper and lower colored fields mirror the two post-duplication copies in Fig. 3d.
    ax.add_patch(plt.Rectangle((0.40, 0.52), 0.58, 0.44, color="#E9ECF7", ec="none", zorder=0))
    ax.add_patch(plt.Rectangle((0.40, 0.02), 0.58, 0.44, color="#FBE9D8", ec="none", zorder=0))
    ax.text(0.96, 0.92, "A", ha="right", va="center", fontsize=8, color=COLORS["A"], fontweight="bold")
    ax.text(0.96, 0.08, "B", ha="right", va="center", fontsize=8, color=COLORS["B"], fontweight="bold")

    # A compact rooted guide tree: GTR is the low-copy outgroup and is not given A/B loss counts.
    x_root, x_stem, x_tip = 0.08, 0.23, 0.39
    y_root, y_gtr, y_prunus, y_ppa, y_psp, y_pzi = 0.50, 0.88, 0.35, 0.56, 0.28, 0.14
    ax.plot([x_root, x_root], [y_gtr, y_prunus], color="#272727", lw=1.0)
    ax.plot([x_root, x_stem], [y_gtr, y_gtr], color="#272727", lw=1.0)
    ax.plot([x_root, x_stem], [y_prunus, y_prunus], color="#272727", lw=1.0)
    ax.plot([x_stem, x_stem], [y_ppa, y_pzi], color="#272727", lw=1.0)
    ax.plot([x_stem, x_tip], [y_ppa, y_ppa], color="#272727", lw=1.0)
    ax.plot([x_stem, x_tip], [y_psp, y_psp], color="#272727", lw=1.0)
    ax.plot([x_stem, x_tip], [y_pzi, y_pzi], color="#272727", lw=1.0)
    for label, y in [("GTR", y_gtr), ("PPA", y_ppa), ("PSP", y_psp), ("PZI", y_pzi)]:
        ax.text(x_tip + 0.015, y, label, va="center", fontsize=6.8)
    wgd_x = (x_root + x_stem) / 2.0
    ax.plot([wgd_x, wgd_x], [y_prunus - 0.035, y_prunus + 0.035], color="#B64342", lw=2.0)
    ax.text(wgd_x, y_prunus - 0.055, "WGD", ha="center", va="top", fontsize=5.5, color="#B64342")

    for branch in BRANCHES:
        y = Y_POS[branch]
        for subgenome, offset in (("A", 0.028), ("B", -0.028)):
            row = by_key[(branch, subgenome)]
            text = "{} ({})".format(row["surviving_genes"], row["losses"])
            ax.text(0.44, y + offset, text, ha="left", va="center", fontsize=6.0, color=COLORS[subgenome])
        ax.text(0.40, y + 0.062, BRANCH_LABELS[branch], ha="left", va="bottom", fontsize=5.4, color="#4D4D4D")

    # Model summary is deliberately separate from the branch counts.
    first = rows[0]
    summary = []
    for key in ("log_likelihood", "copy2_bias", "switch_rate"):
        if key in first and first[key]:
            summary.append("{}={}".format(key, first[key]))
    if summary:
        ax.text(0.40, 0.485, "  ".join(summary), ha="left", va="center", fontsize=5.3, color="#4D4D4D")
    ax.text(0.40, 0.005, "surviving genes (losses)", ha="left", va="bottom", fontsize=5.3, color="#4D4D4D")
    ax.text(0.02, 0.98, "d", fontsize=10, fontweight="bold", ha="left", va="top")
    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    prefix = outdir / "figure2_panel_d_point_reference_style"
    fig.savefig(str(prefix) + ".svg", bbox_inches="tight")
    fig.savefig(str(prefix) + ".pdf", bbox_inches="tight")
    fig.savefig(str(prefix) + ".tiff", dpi=600, bbox_inches="tight")
    fig.savefig(str(prefix) + ".png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    write_tsv(outdir / "panel_d_source_data.tsv", list(rows[0].keys()), rows)
    (outdir / "COMPLETE").write_text("Formal POInT likelihood panel plotted.\n", encoding="ascii")
    print("[OK] Wrote panel d: {}".format(outdir))


if __name__ == "__main__":
    main()
