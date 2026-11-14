from pathlib import Path
import csv, json, collections, re
import numpy as np
from Bio import SeqIO

import argparse

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--project-root", type=Path, required=True)
parser.add_argument("--out", type=Path, required=True)
args = parser.parse_args()
R = args.project_root / "0903"
O = args.out
C = O / "circos"
D = C / "data"
D.mkdir(parents=True, exist_ok=True)
(C / "figures").mkdir(exist_ok=True)

old = R / "haploid_repartition/circos/pzi_psp_ppa"
names = [
    "gene_density",
    "te_coverage",
    "gypsy_coverage",
    "copia_coverage",
    "tir_coverage",
    "helitron_coverage",
    "gc_fraction",
    "ncrna_density",
]
B = 1000000
allrows = []
chroms = []
checks = []
summ = []
pzi_ncrna = collections.defaultdict(lambda: [0, 0])
# Map ncRNAs from original chromosome coordinates to the corrected PZI phase assignment.
pzi_map = {
    v["original_chr32"]: v
    for v in csv.DictReader(
        (R / "haploid_repartition/maps/pzi_repartition_map.tsv").open(), delimiter="\t"
    )
}
nc = collections.defaultdict(list)
for line in (R / "annotation/ncrna_clean/pzi/final/pzi.ncRNA.final.gff3").open():
    v = line.rstrip().split("\t")
    if len(v) != 9 or v[0] not in pzi_map:
        continue
    m = pzi_map[v[0]]
    a, b = int(v[3]) - 1, int(v[4])
    attrs = dict(x.split("=", 1) for x in v[8].split(";") if "=" in x)
    cl = attrs.get("ncRNA_class", v[2])
    h = m["haplotype"]
    nc[h, m["new_seqid"]].append((a, b))
    pzi_ncrna[h, cl][0] += 1
    pzi_ncrna[h, cl][1] += b - a
for sp in ["pzi", "psp", "ppa"]:
    out = {n: [] for n in names}
    for h in ["hapA", "hapB"] if sp == "ppa" else ["hapA1", "hapB1"]:
        base = R / f"haploid_repartition/inputs/{sp}"
        fa = base / f"genome/{sp}_{h}.fa"
        seqs = {x.id: str(x.seq).upper() for x in SeqIO.parse(fa, "fasta")}
        lengths = {c: len(s) for c, s in seqs.items()}
        edges = {c: np.r_[np.arange(0, L, B), L] for c, L in lengths.items()}
        values = {n: {c: np.zeros(len(e) - 1) for c, e in edges.items()} for n in names}
        outside = []
        for c, s in seqs.items():
            for i, (a, b) in enumerate(zip(edges[c][:-1], edges[c][1:])):
                sub = s[a:b]
                valid = sum(sub.count(q) for q in "ACGT")
                values["gc_fraction"][c][i] = (
                    (sub.count("G") + sub.count("C")) / valid if valid else 0
                )
            chroms.append(
                dict(
                    species=sp,
                    haplotype=h,
                    chromosome=c,
                    length_bp=len(s),
                    source_fasta=str(fa),
                )
            )
        del seqs
        for line in (base / f"gff/{sp}_{h}.gff3").open():
            v = line.rstrip().split("\t")
            if len(v) != 9 or v[2] != "gene" or v[0] not in lengths:
                continue
            c = v[0]
            a, b = int(v[3]) - 1, int(v[4])
            assert 0 <= a < b <= lengths[c], (sp, h, v[:5])
            values["gene_density"][c][((a + b) // 2) // B] += 1
        for c, L in lengths.items():
            if sp == "pzi":
                intervals = nc[h, c]
            else:
                intervals = []
            if sp == "pzi":
                for a, b in intervals:
                    assert 0 <= a < b <= L, (sp, h, c, a, b, L)
                    values["ncrna_density"][c][((a + b) // 2) // B] += 1
        if sp != "pzi":
            nf = (
                R
                / f"annotation/ncrna_haploid_clean/{sp}/{h}/final/{sp}_{h}.ncRNA.final.gff3"
            )
            for line in nf.open():
                v = line.rstrip().split("\t")
                if len(v) != 9 or v[0] not in lengths:
                    continue
                c = v[0]
                a, b = int(v[3]) - 1, int(v[4])
                assert 0 <= a < b <= lengths[c], (sp, h, c, a, b, lengths[c])
                values["ncrna_density"][c][((a + b) // 2) // B] += 1
        intervals = collections.defaultdict(list)
        for line in (base / f"edta/{h}/{sp}_{h}.fa.mod.EDTA.TEanno.gff3").open():
            v = line.rstrip().split("\t")
            if len(v) != 9 or v[0] not in lengths:
                continue
            c = v[0]
            a, b = int(v[3]) - 1, int(v[4])
            L = lengths[c]
            if a < 0 or b > L:
                outside.append((c, a, b, L))
            a = max(0, a)
            b = min(L, b)
            if a >= b:
                continue
            ats = dict(t.split("=", 1) for t in v[8].split(";") if "=" in t)
            cls = ats.get("classification", "")
            cats = ["te_coverage"]
            if "Gypsy" in cls:
                cats.append("gypsy_coverage")
            elif "Copia" in cls:
                cats.append("copia_coverage")
            elif "Helitron" in cls:
                cats.append("helitron_coverage")
            elif cls.startswith("DNA/") or "TIR" in cls:
                cats.append("tir_coverage")
            for n in cats:
                intervals[n, c].append((a, b))
        for (n, c), ints in intervals.items():
            merged = []
            for a, b in sorted(ints):
                if merged and a <= merged[-1][1]:
                    merged[-1] = (merged[-1][0], max(b, merged[-1][1]))
                else:
                    merged.append((a, b))
            for a, b in merged:
                for k in range(a // B, (b - 1) // B + 1):
                    values[n][c][k] += min(b, edges[c][k + 1]) - max(a, edges[c][k])
            values[n][c] /= np.diff(edges[c])
        for n in names:
            for c, e in edges.items():
                vs = values[n][c]
                if n in ["gene_density", "ncrna_density"]:
                    vs = vs / (np.diff(e) / B)
                for a, b, v in zip(e[:-1], e[1:], vs):
                    out[n].append(f"{sp}_{c} {int(a)} {int(b)} {v:.8f}")
                    allrows.append(
                        dict(
                            species=sp,
                            haplotype=h,
                            chromosome=c,
                            track=n,
                            start0=int(a),
                            end=int(b),
                            value=float(v),
                            unit=(
                                "features_per_Mb"
                                if n in ["gene_density", "ncrna_density"]
                                else "fraction"
                            ),
                        )
                    )
                summ.append(
                    dict(
                        species=sp,
                        haplotype=h,
                        chromosome=c,
                        track=n,
                        windows=len(vs),
                        sequence_length_bp=lengths[c],
                        min=float(vs.min()),
                        max=float(vs.max()),
                        weighted_mean=float(np.average(vs, weights=np.diff(e))),
                    )
                )
        checks.append(
            dict(
                species=sp,
                haplotype=h,
                EDTA_boundary_clipped_records=len(outside),
                gene_records=int(sum(values["gene_density"][c].sum() for c in lengths)),
                ncrna_records=int(
                    sum(values["ncrna_density"][c].sum() for c in lengths)
                ),
            )
        )
    for n, lines in out.items():
        (D / f"{sp}_{n}.txt").write_text("\n".join(lines) + "\n")
    print("Finished tracks", sp, flush=True)
for name, rows in [
    ("track_values.tsv", allrows),
    ("track_summary.tsv", summ),
    ("chromosomes.tsv", chroms),
]:
    with (O / name).open("w") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        wr.writeheader()
        wr.writerows(rows)
with (O / "pzi_ncrna_corrected_counts.tsv").open("w") as f:
    wr = csv.writer(f, delimiter="\t")
    wr.writerow(["haplotype", "ncRNA_class", "count", "length_bp"])
    wr.writerows([h, cl, *v] for (h, cl), v in sorted(pzi_ncrna.items()))
(O / "track_validation.json").write_text(json.dumps(checks, indent=2))
lengths = {f"{x['species']}_{x['chromosome']}": x["length_bp"] for x in chroms}
kt = []
for line in (old / "data/karyotype.txt").read_text().splitlines():
    v = line.split()
    v[5] = str(lengths[v[2]])
    kt.append(" ".join(v))
(D / "karyotype.txt").write_text("\n".join(kt) + "\n")
config = (old / "conf/circos.conf").read_text()
m = re.search(r"^chromosomes_scale = (.*)$", config, re.M)
oldscales = {a: float(b[:-1]) for a, b in [x.split("=") for x in m.group(1).split(",")]}
total = {
    sp: sum(v for c, v in oldscales.items() if c.startswith(sp + "_"))
    for sp in ["pzi", "psp", "ppa"]
}
scale = {
    c: total[c.split("_")[0]]
    * L
    / sum(v for cc, v in lengths.items() if cc.startswith(c.split("_")[0] + "_"))
    for c, L in lengths.items()
}
config = (
    config[: m.start(1)]
    + ",".join(f"{c}={scale[c]:.10f}r" for c in oldscales)
    + config[m.end(1) :]
)
maxgene = max(x["value"] for x in allrows if x["track"] == "gene_density")
config = config.replace("max = 215.000000", f"max = {maxgene:.6f}")
(O / "max_gene_density.txt").write_text(str(maxgene))
(C / "circos.conf").write_text(config)
print("Max gene density", maxgene, checks)
