from pathlib import Path
import gzip, re
import pandas as pd

import argparse

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--canonical", type=Path, required=True)
parser.add_argument("--pairs", type=Path, required=True)
parser.add_argument("--out", type=Path, required=True)
args = parser.parse_args()
R = args.out
R.mkdir(parents=True, exist_ok=True)

maps = {}
for sp in ["pcam", "ppe"]:
    p = args.canonical / f"{sp}/primary.gff3"
    tx = {}
    prots = []
    with p.open("rb") as f:
        gz = f.read(2) == b"\x1f\x8b"
    with gzip.open(p, "rt") if gz else p.open() as f:
        for line in f:
            if line.startswith("#"):
                continue
            v = line.rstrip().split("\t")
            if len(v) != 9:
                continue
            a = dict(x.split("=", 1) for x in v[8].split(";") if "=" in x)
            if v[2].lower() in ["mrna", "transcript"]:
                tx[a["ID"]] = a.get("Parent", a["ID"])
            if v[2] == "CDS" and a.get("protein_id"):
                prots.append((a["protein_id"], a["Parent"]))
    maps[sp] = {pr: tx[t] for pr, t in prots}


def loc(sp, g):
    if sp in maps:
        return maps[sp][g]
    return re.sub(r"-mRNA-\d+$|\.t\d+$", "", g)


p = args.pairs
d = pd.read_csv(p, sep="\t")
d["query_locus"] = [loc(s, g) for s, g in zip(d.query_label, d.query_gene)]
d["ref_locus"] = [loc(s, g) for s, g in zip(d.ref_label, d.ref_gene)]
rows = []
blocks = []
for comp, z in d.groupby("comparison", sort=False):
    group = z.groupby(["query_locus", "ref_locus"], sort=False)
    collapsed = group.Ks.median().reset_index()
    for k in [
        "source_class",
        "panel",
        "curve",
        "comparison",
        "query_label",
        "ref_label",
    ]:
        collapsed[k] = z.iloc[0][k]
    collapsed["isoform_pair_records"] = group.size().values
    blocks.append(collapsed)
    rows.append(
        dict(
            comparison=comp,
            records=len(z),
            distinct_locus_pairs=len(collapsed),
            extra_records=len(z) - len(collapsed),
            fraction_extra=1 - len(collapsed) / len(z),
            old_median_Ks=z.Ks.median(),
            locus_pair_median_Ks=collapsed.Ks.median(),
        )
    )
pd.DataFrame(rows).to_csv(R / "archived_Ks_locus_pair_audit.tsv", sep="\t", index=False)
pd.concat(blocks).to_csv(R / "Ks_locus_pair_sensitivity.tsv", sep="\t", index=False)
print(pd.DataFrame(rows).to_string(index=False))
