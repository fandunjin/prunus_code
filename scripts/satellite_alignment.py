"""Five observed species representatives, never synthetic cross-species consensus."""

from pathlib import Path
from collections import Counter
import csv, json, subprocess, hashlib
import numpy as np, pandas as pd
from Bio import SeqIO, Align
from Bio.Seq import Seq
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch

import argparse

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--census-dir", type=Path, required=True)
parser.add_argument("--source", type=Path, required=True)
parser.add_argument("--out", type=Path, required=True)
args = parser.parse_args()
R = args.census_dir
S = args.source
O = args.out
O.mkdir(parents=True, exist_ok=True)

groups = {
    "pzi": ["pziA1", "pziA2", "pziB1", "pziB2"],
    "psp": ["pspA1", "pspA2", "pspB1", "pspB2"],
    "ppa": ["ppaA", "ppaB"],
    "pcam": ["pcam"],
    "ppe": ["ppe"],
}
labels = {
    "pzi": "P. zippeliana",
    "psp": "P. spinulosa",
    "ppa": "P. padus",
    "pcam": "P. campanulata",
    "ppe": "P. persica",
    "PlantSat_X82221": "PlantSat / X82221",
    "X82222": "X82222",
}
reps = {}
metadata = []
for sp, targets in groups.items():
    counts = Counter()
    total = 0
    per_target = {}
    for t in targets:
        c = Counter()
        for row in csv.DictReader(
            (R / t / "sequence_classification.tsv").open(), delimiter="\t"
        ):
            if row["published_homology"] == "True":
                c[row["sequence"]] += int(row["occurrences"])
        counts.update(c)
        total += sum(c.values())
        per_target[t] = c
    sequence, n = sorted(counts.items(), key=lambda x: (-x[1], x[0]))[0]
    reps[sp] = sequence
    t = max(targets, key=lambda x: per_target[x][sequence])
    example = None
    for row in csv.DictReader(
        (S / "trash" / t / f"{t}.fa_repeats_with_seq.csv").open()
    ):
        if row["sequence"].upper() == sequence:
            example = row
            break
    assert example
    metadata.append(
        dict(
            id=sp,
            species=labels[sp],
            length=len(sequence),
            exact_sequence_occurrences=n,
            all_supported_family_occurrences=total,
            exact_sequence_fraction=n / total,
            haplotypes_pooled=";".join(targets),
            example_target=t,
            example_chromosome=example["seqID"],
            example_start1=example["start"],
            example_end1=example["end"],
            original_sequence=sequence,
            sha256=hashlib.sha256(sequence.encode()).hexdigest(),
        )
    )
    print(sp, n, total, len(sequence), flush=True)
refs = {
    r.id: str(r.seq) for r in SeqIO.parse(S / "analysis/published_monomers.fa", "fasta")
}
reference = refs["ref_Prunus170_consensus|ref|ref|w166"]
reps["PlantSat_X82221"] = reference
reps["X82222"] = refs["ref_Prunus170_X82222|ref|ref|w166"]
a = Align.PairwiseAligner()
a.mode = "global"
a.match_score = 2
a.mismatch_score = -1
a.open_gap_score = -3
a.extend_gap_score = -0.5
oriented = {}
for name, s in reps.items():
    candidates = []
    for strand, seq in [("+", s), ("-", str(Seq(s).reverse_complement()))]:
        for shift in range(len(seq)):
            rotated = seq[shift:] + seq[:shift]
            candidates.append(
                (a.score(reference, rotated), strand == "+", -shift, rotated, strand)
            )
    score, _, negshift, rotated, strand = max(candidates)
    oriented[name] = rotated
    for row in metadata:
        if row["id"] == name:
            row.update(
                strand_for_display=strand,
                circular_left_shift=-negshift,
                phase_alignment_score=score,
                display_sequence=rotated,
            )
with (O / "representatives_original.fa").open("w") as f:
    for k, s in reps.items():
        f.write(f">{k}\n{s}\n")
with (O / "representatives_phase_normalized.fa").open("w") as f:
    for k, s in oriented.items():
        f.write(f">{k}\n{s}\n")
cmd = [
    "mafft",
    "--localpair",
    "--maxiterate",
    "1000",
    "--thread",
    "1",
    str(O / "representatives_phase_normalized.fa"),
]
with (O / "representatives_aligned.fa").open("w") as out, (O / "mafft.log").open(
    "w"
) as err:
    subprocess.run(cmd, check=True, stdout=out, stderr=err)
aln = {
    r.id: str(r.seq).upper()
    for r in SeqIO.parse(O / "representatives_aligned.fa", "fasta")
}
length = len(next(iter(aln.values())))
assert all(len(s) == length for s in aln.values())
rows = []
species = list(groups)
for i in range(length):
    letters = [aln[k][i] for k in species]
    c = Counter(x for x in letters if x in "ACGT")
    base, n = c.most_common(1)[0] if c else ("-", 0)
    rows.append(
        dict(
            alignment_column=i + 1,
            consensus_base=base,
            maximum_agreement=n,
            all_five_identical=n == 5,
            four_of_five_identical=n == 4,
            number_gaps=letters.count("-"),
        )
    )
cons = pd.DataFrame(rows)
pd.DataFrame(metadata).to_csv(O / "representative_selection.tsv", sep="\t", index=False)
cons.to_csv(O / "alignment_conservation.tsv", sep="\t", index=False)
# Contiguous identical blocks are descriptive, with no post-hoc minimum run filter.
runs = []
start = None
for i, passed in enumerate(list(cons.all_five_identical) + [False]):
    if passed and start is None:
        start = i
    if not passed and start is not None:
        runs.append(
            dict(
                start_column=start + 1,
                end_column=i,
                length=i - start,
                sequence="".join(cons.consensus_base.iloc[start:i]),
            )
        )
        start = None
pd.DataFrame(runs).to_csv(O / "fully_conserved_runs.tsv", sep="\t", index=False)
plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 7,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
    }
)
order = list(reps)
chunks = [(i, min(i + 60, length)) for i in range(0, length, 60)]
fig, axes = plt.subplots(len(chunks), 1, figsize=(10.5, 1.9 * len(chunks)))
fig.subplots_adjust(left=0.18, right=0.99, top=0.9, bottom=0.09, hspace=0.55)
for ax, (lo, hi) in zip(np.atleast_1d(axes), chunks):
    for j in range(lo, hi):
        n = int(cons.iloc[j].maximum_agreement)
        if n >= 4:
            ax.add_patch(
                Rectangle(
                    (j - lo - 0.47, -0.48),
                    0.94,
                    5,
                    facecolor="#b8dccc" if n == 5 else "#f7e5b3",
                    edgecolor="none",
                    zorder=0,
                )
            )
        for y, k in enumerate(order):
            ax.text(
                j - lo,
                y,
                aln[k][j],
                fontfamily="DejaVu Sans Mono",
                fontsize=7,
                ha="center",
                va="center",
                color="#1d1d1d" if y < 5 else "#666666",
            )
    ax.set_xlim(-0.7, 59.7)
    ax.set_ylim(len(order) - 0.4, -0.7)
    ax.set_yticks(range(len(order)), [labels[x] for x in order], fontsize=8)
    ax.set_xticks(
        [j - lo for j in range(lo, hi) if (j + 1) % 10 == 0 or j == lo],
        [str(j + 1) for j in range(lo, hi) if (j + 1) % 10 == 0 or j == lo],
        fontsize=7,
    )
    ax.xaxis.tick_top()
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
fig.legend(
    handles=[
        Patch(
            facecolor="#b8dccc", label="Identical in all five species representatives"
        ),
        Patch(facecolor="#f7e5b3", label="Identical in four of five"),
    ],
    ncol=2,
    frameon=False,
    loc="lower center",
    bbox_to_anchor=(0.55, 0.02),
    fontsize=8,
)
for ext in ["pdf", "svg", "png"]:
    fig.savefig(O / f"representatives.{ext}", dpi=180, bbox_inches="tight")
plt.close(fig)
summary = dict(
    selection="Most frequent exact observed homology-supported monomer per species pooled across available haplotypes; exact raw phases are not pooled across rotations",
    alignment_command=cmd,
    phase_normalization="Both strands and all circular shifts scored globally against PlantSat/X82221; score 2/-1, gap open -3, extend -0.5; deterministic ties prefer plus strand then smaller shift",
    n_species=5,
    alignment_columns=length,
    all_five_identical_columns=int(cons.all_five_identical.sum()),
    four_of_five_identical_columns=int(cons.four_of_five_identical.sum()),
    caveat="Agreement among five chosen representatives is not site conservation across all within-species copies, and not functional validation",
    reference_identity="PlantSat consensus and X82221 are identical in the archived reference FASTA",
)
(O / "alignment_summary.json").write_text(json.dumps(summary, indent=2))
print(summary, flush=True)
