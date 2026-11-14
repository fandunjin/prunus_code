"""Render corrected Fig. 2B through native JCVI Karyotype with physical lengths.
Prepared anchors and positions retain the original seven-track workflow.
"""

from pathlib import Path
import csv, re, json, os
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from jcvi.graphics.karyotype import Karyotype

import argparse

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--project-root", type=Path, required=True)
parser.add_argument("--out", type=Path, required=True)
args = parser.parse_args()
O = args.out
D = O / "jcvi"
D.mkdir(exist_ok=True)
old = args.project_root / "0903/native_multispecies_A1B1_pcam_pairwise"

rows = list(csv.DictReader((O / "riparian_layout.tsv").open(), delimiter="\t"))
blocks = [
    r
    for r in csv.DictReader((O / "all_display_blocks.tsv").open(), delimiter="\t")
    if r["pair"] != "pzi-ppa"
]
ss = ["pzi", "psp", "ppa", "pcam", "ppe", "gtr", "rch"]
names = [
    "P. zippeliana",
    "P. spinulosa",
    "P. padus",
    "P. campanulata",
    "P. persica",
    "G. trifoliata",
    "R. chinensis",
]
colors = [
    "#0072B2",
    "#E69F00",
    "#009E73",
    "#D55E00",
    "#CC79A7",
    "#56B4E9",
    "#BEA800",
    "#505050",
]
sizes = {r["species"] + "_" + r["chromosome"]: int(r["length_bp"]) for r in rows}
allbeds = {}
nl = []
sl = []
for sp, name in zip(ss, names):
    src = O / f"{sp}.bed" if sp in ["pzi", "psp", "ppa"] else old / f"bed/{sp}.bed"
    data = {v[3]: v for v in [s.split() for s in src.read_text().splitlines()]}
    if sp == "ppa":
        for line in (old / "bed/ppa.bed").read_text().splitlines():
            v = line.split()
            data.setdefault(v[3], v)
    allbeds[sp] = data
    with (D / f"{sp}.bed").open("w") as f:
        for v in sorted(data.values(), key=lambda v: (v[0].split("__")[-1], int(v[1]))):
            c = v[0].split("__")[-1]
            key = sp + "_" + c
            if key not in sizes:
                continue
            assert int(v[2]) <= sizes[key], (key, v[2], sizes[key])
            f.write("\t".join([key] + v[1:]) + "\n")
    rs = [r for r in rows if r["species"] == sp]
    span = (
        0.74
        * sum(int(r["length_bp"]) for r in rs)
        / sum(int(r["length_bp"]) for r in rows if r["species"] == "pzi")
    )
    x0 = 0.43 - span / 2
    x1 = 0.43 + span / 2
    y = float(rs[0]["y"])
    nl.append(f"{y}, {x0}, {x1}, 0, #4C5660, {name}, center, {sp}.bed, center")
    sl.append(
        ",".join(
            sp + "_" + r["chromosome"] + ("-" if r["display_flipped"] == "True" else "")
            for r in rs
        )
    )
for i, (a, b) in enumerate(zip(ss[:-1], ss[1:])):
    bb = [r for r in blocks if r["pair"] == a + "-" + b]
    with (D / f"{a}.{b}.simple").open("w") as f:
        for r in bb:
            gs = [
                r[k]
                for k in [
                    "left_start_gene",
                    "left_end_gene",
                    "right_start_gene",
                    "right_end_gene",
                ]
            ]
            for gene, sp, c in [
                (gs[0], a, r["left_chromosome"]),
                (gs[1], a, r["left_chromosome"]),
                (gs[2], b, r["right_chromosome"]),
                (gs[3], b, r["right_chromosome"]),
            ]:
                assert allbeds[sp][gene][0].split("__")[-1] == c
            col = colors[(int(re.search(r"\d+", r["left_chromosome"]).group()) - 1) % 8]
            f.write(col + "*" + "\t".join(gs + [r["anchors"], r["orientation"]]) + "\n")
    nl.append(f"e, {i}, {i+1}, {a}.{b}.simple")
(D / "layout").write_text("\n".join(nl) + "\n")
(D / "seqids").write_text("\n".join(sl) + "\n")
(D / "sizes.json").write_text(json.dumps(sizes, indent=2))
os.chdir(D)
plt.rcParams.update(
    {
        "text.usetex": False,
        "font.family": "DejaVu Sans",
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
        "font.size": 18,
    }
)
fig = plt.figure(figsize=(18, 24))
ax = fig.add_axes([0, 0, 1, 1])
k = Karyotype(
    ax,
    "seqids",
    "layout",
    generank=False,
    sizes=sizes,
    gap=0.003,
    height=0.008,
    lw=0.6,
    plot_label=False,
    plot_circles=False,
    keep_chrlabels=True,
    chrstyle="roundrect",
    shadestyle="curve",
)
for tr, sp, name in zip(k.tracks, ss, names):
    ax.text(
        max(0.83, tr.xend + 0.012),
        tr.y,
        name,
        fontsize=17,
        fontstyle="italic",
        ha="left",
        va="center",
        clip_on=False,
    )
    if sp == "pzi":
        for c in tr.seqids:
            ax.text(
                tr.offsets[c] + tr.sizes[c] * tr.ratio / 2,
                tr.y + 0.012,
                c.removeprefix("pzi_"),
                fontsize=16,
                rotation=60,
                ha="left",
                va="bottom",
            )
for r in rows:
    tr = k.tracks[ss.index(r["species"])]
    c = r["species"] + "_" + r["chromosome"]
    r["xstart"] = tr.offsets[c]
    r["xend"] = tr.offsets[c] + tr.sizes[c] * tr.ratio
with (O / "riparian_layout.tsv").open("w") as f:
    wr = csv.DictWriter(f, fieldnames=rows[0].keys(), delimiter="\t")
    wr.writeheader()
    wr.writerows(rows)
ax.plot([0.27, 0.81], [0.60, 0.60], ls=(0, (6, 6)), lw=1.8, color="#555555")
ax.text(0.84, 0.60, "★", fontsize=30, color="#c33f45", ha="center", va="center")
# Leave room at right for species labels while retaining the 18x24 source canvas.
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis("off")
for fmt in ["pdf", "svg", "png"]:
    fig.savefig(O / f"Fig2b_correct_A1B1.{fmt}", dpi=150)
(O / "jcvi_render_validation.json").write_text(
    json.dumps(
        {
            "renderer": "jcvi.graphics.karyotype.Karyotype",
            "generank": False,
            "explicit_sizes": True,
            "blocks": len(blocks),
            "tracks": len(k.tracks),
            "all_endpoints_valid": True,
        },
        indent=2,
    )
)
print("Native JCVI rendered:", len(blocks), "ribbons;", len(sizes), "chromosomes")
