"""Decompose genome-size differences and plot the satellite contribution (SF17)."""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--composition", type=Path, required=True)
parser.add_argument("--summary", type=Path, required=True)
parser.add_argument("--out", type=Path, required=True)
args = parser.parse_args()
args.out.mkdir(parents=True, exist_ok=True)
data = pd.read_csv(args.composition, sep="\t")
summary = pd.read_csv(args.summary, sep="\t").set_index("target")
data["component"] = np.where(
    data.subfamily == "Prunus_170 homologues (163-170 bp)",
    "PziCent167 family",
    np.where(
        data.repeat_class == "Satellite", "Other tandem monomers", data.repeat_class
    ),
)
wide = data.groupby(["component", "target"]).bp.sum().unstack(fill_value=0)
wide.loc["Remaining sequence"] = summary.genome_bp - wide.sum(axis=0)
categories = [
    "PziCent167 family",
    "Other tandem monomers",
    "LTR/Gypsy",
    "LTR/Copia",
    "Other LTR",
    "TIR/DNA",
    "Helitron",
    "LINE",
    "SINE",
    "Other TE",
    "Remaining sequence",
]
categories = [x for x in categories if x in wide.index] + [
    x for x in wide.index if x not in categories
]
rows = []
for hap in ["A1", "B1"]:
    for component in categories:
        pzi = int(wide.loc[component, "PZI-" + hap])
        psp = int(wide.loc[component, "PSP-" + hap])
        rows.append(
            dict(
                haplotype=hap,
                component=component,
                PZI_bp=pzi,
                PSP_bp=psp,
                difference_bp=pzi - psp,
            )
        )
result = pd.DataFrame(rows)
for hap in ["A1", "B1"]:
    assert (
        result[result.haplotype.eq(hap)].difference_bp.sum()
        == summary.loc["PZI-" + hap, "genome_bp"]
        - summary.loc["PSP-" + hap, "genome_bp"]
    )
result.to_csv(
    args.out / "genome_size_increment_decomposition.tsv", sep="\t", index=False
)
groups = {
    "PziCent167 family": "PziCent167",
    "Other tandem monomers": "Other tandem",
    "LTR/Gypsy": "LTR",
    "LTR/Copia": "LTR",
    "Other LTR": "LTR",
    "Remaining sequence": "Remaining sequence",
}
result["plot_component"] = result.component.map(groups).fillna("Other TE")
values = result.groupby(["haplotype", "plot_component"]).difference_bp.sum() / 1e6
labels = ["PziCent167", "Other tandem", "LTR", "Other TE", "Remaining sequence"]
colors = ["#377eb8", "#9baec2", "#c9a77a", "#baa1c8", "#d6d6d6"]
plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 7,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
    }
)
fig, ax = plt.subplots(figsize=(6.2, 1.72))
fig.subplots_adjust(left=0.105, right=0.96, bottom=0.42, top=0.83)
for y, hap in enumerate(["A1", "B1"]):
    total = values.loc[hap].sum()
    left = 0
    for label, color in zip(labels, colors):
        amount = values.get((hap, label), 0)
        ax.barh(
            y,
            amount,
            left=left,
            height=0.5,
            color=color,
            edgecolor="white",
            linewidth=0.3,
            label=label if y == 0 else None,
        )
        if label == "PziCent167":
            ax.text(
                left + amount / 2,
                y,
                f"{amount:.1f} Mb ({100*amount/total:.1f}%)",
                ha="center",
                va="center",
                color="white",
            )
        left += amount
    ax.text(left + 2, y, f"{left:.1f} Mb", va="center")
ax.set_yticks([0, 1], ["A1", "B1"])
ax.invert_yaxis()
ax.set_xlim(0, 210)
ax.set_xticks([0, 50, 100, 150, 200])
ax.set_xlabel("Genome-size difference: P. zippeliana − P. spinulosa (Mb)")
ax.spines[["top", "right", "left"]].set_visible(False)
ax.tick_params(axis="y", length=0)
ax.legend(
    loc="upper center",
    bbox_to_anchor=(0.5, -0.65),
    ncol=3,
    frameon=False,
    fontsize=6.5,
    handlelength=1.2,
    columnspacing=1,
)
for extension in ["pdf", "svg"]:
    fig.savefig(
        args.out / f"Supplemental_Figure17.{extension}",
        bbox_inches="tight",
        pad_inches=0.06,
    )
plt.close(fig)
