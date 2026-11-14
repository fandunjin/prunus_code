"""Rebuild eight-track synteny inputs with explicit GFF parent-locus identity."""

from pathlib import Path
import csv, re, json, gzip
from collections import defaultdict
from primary_locus_proteomes import fasta, write, sha

import argparse

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--catalog", type=Path, required=True)
parser.add_argument("--canonical", type=Path, required=True)
parser.add_argument("--primary-inputs", type=Path, required=True)
parser.add_argument("--out", type=Path, required=True)
args = parser.parse_args()
src = args.catalog
out = args.out
if (out / "inventory.tsv").exists():
    raise SystemExit("Completed output exists")
out.mkdir(exist_ok=True)
labels = {
    "PZI_A": "pzi_hapA1",
    "PZI_B": "pzi_hapB1",
    "PSP_A": "psp_hapA1",
    "PSP_B": "psp_hapB1",
    "PPA_A": "ppa_hapA",
    "PPA_B": "ppa_hapB",
    "PCAM": "pcam",
    "PPE": "ppe",
    "GTR": "gtr",
}


def attrs(s):
    return dict(x.split("=", 1) for x in s.strip().split(";") if "=" in x)


def protein_parents(tr):
    p = args.canonical / f"{tr.lower()}/primary.gff3"
    tx = {}
    prots = []
    with p.open("rb") as probe:
        compressed = probe.read(2) == b"\x1f\x8b"
    with gzip.open(p, "rt") if compressed else p.open() as f:
        for l in f:
            if l.startswith("#"):
                continue
            v = l.rstrip().split("\t")
            if len(v) != 9:
                continue
            a = attrs(v[8])
            if v[2].lower() in ["mrna", "transcript"]:
                tx[a["ID"]] = a.get("Parent", a["ID"])
            elif v[2] == "CDS" and a.get("protein_id"):
                prots.append((a["protein_id"], a["Parent"]))
    result = {}
    for pr, t in prots:
        g = tx[t]
        assert pr not in result or result[pr] == g, (tr, pr)
        result[pr] = g
    return result


report = []
for tr, label in labels.items():
    pep = dict((h.split()[0], s) for h, s in fasta(src / f"{tr}.pep"))
    cds = dict((h.split()[0], s) for h, s in fasta(src / f"{tr}.cds"))
    bed = [l.split("\t") for l in (src / f"{tr}.bed").read_text().splitlines() if l]
    gmap = protein_parents(tr) if tr in ["PCAM", "PPE"] else {}
    grouped = defaultdict(list)
    invalid = 0
    for v in bed:
        key = v[3]
        raw = key.split("__", 1)[1]
        g = gmap[raw] if gmap else re.sub(r"(?:-mRNA-\d+|\.t\d+)$", "", raw)
        assert key in pep and key in cds, (tr, key)
        if len(cds[key]) < 300 or "*" in pep[key]:
            invalid += 1
            continue
        grouped[g].append((key, v))
    selected = []
    audit = []
    for g, items in grouped.items():
        key, v = min(items, key=lambda x: (-len(cds[x[0]]), -len(pep[x[0]]), x[0]))
        selected.append((key, v))
        audit.append(
            dict(
                track=tr,
                gene_locus=g,
                representative=key,
                cds_length=len(cds[key]),
                alternative_records=len(items) - 1,
            )
        )
    selected.sort(key=lambda x: (x[1][0], int(x[1][1]), int(x[1][2]), x[0]))
    dest = out / label
    dest.mkdir(exist_ok=True)
    with (dest / f"{label}.bed").open("w") as f:
        for key, v in selected:
            vv = list(v)
            vv[0] = vv[0].split("__", 1)[1]
            vv[3] = key.split("__", 1)[1]
            f.write("\t".join(vv) + "\n")
    for ext, records in [("pep", pep), ("cds", cds)]:
        with (dest / f"{label}.{ext}").open("w") as f:
            for key, v in selected:
                f.write(">" + key.split("__", 1)[1] + "\n" + records[key] + "\n")
    seqids = list(dict.fromkeys(v[0].split("__", 1)[1] for key, v in selected))
    (dest / f"{label}.seqids").write_text(",".join(seqids) + "\n")
    write(out / f"{tr}.locus_map.tsv", audit, list(audit[0]))
    if tr.startswith(("PZI", "PSP", "PPA")):
        tr2 = tr.lower().replace("_a", "_A").replace("_b", "_B")
        with (args.primary_inputs / f"id_maps/{tr2}.id_map.tsv").open() as f:
            fam = {r["original_id"] for r in csv.DictReader(f, delimiter="\t")}
        assert fam == {k.split("__", 1)[1] for k, v in selected}, (
            tr,
            "family/synteny representatives differ",
        )
    rr = dict(
        track=tr,
        input_records=len(bed),
        primary_loci=len(selected),
        short_or_invalid_records=invalid,
        pep_sha256=sha(dest / f"{label}.pep"),
    )
    report.append(rr)
    print(json.dumps(rr), flush=True)
write(out / "inventory.tsv", report, list(report[0]))
