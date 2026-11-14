from pathlib import Path
import csv, json
from collections import defaultdict

import argparse

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--source", type=Path, required=True)
parser.add_argument("--census-dir", type=Path, required=True)
parser.add_argument("--genomes", type=Path, required=True)
parser.add_argument("--centromeres", type=Path, required=True)
parser.add_argument("--out", type=Path, required=True)
args = parser.parse_args()
R = args.census_dir
T = args.source / "trash"
O = args.out
O.mkdir(parents=True, exist_ok=True)
gs = {
    r["target"]: sum(
        int(x.split("\t")[1]) for x in Path(r["fai"]).read_text().splitlines()
    )
    for r in csv.DictReader(args.genomes.open(), delimiter="\t")
}


def merge(xs):
    out = []
    for s, e in sorted(xs):
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(e, out[-1][1]))
        else:
            out.append((s, e))
    return out


def size(d):
    return sum(e - s for xs in d.values() for s, e in xs)


def intersect(a, b):
    out = {}
    for c, xs in a.items():
        ys = b.get(c, [])
        j = 0
        z = []
        for s, e in xs:
            while j < len(ys) and ys[j][1] <= s:
                j += 1
            k = j
            while k < len(ys) and ys[k][0] < e:
                lo = max(s, ys[k][0])
                hi = min(e, ys[k][1])
                if hi > lo:
                    z.append((lo, hi))
                k += 1
        out[c] = merge(z)
    return out


rs = []
for raw in csv.DictReader((R / "census.tsv").open(), delimiter="\t"):
    t = raw["target"]
    arr = defaultdict(list)
    mono = defaultdict(list)
    fam = defaultdict(list)
    boundary = 0
    for r in csv.DictReader((T / t / f"{t}.fa_arrays.csv").open()):
        lo, hi = sorted([int(float(r["start"])), int(float(r["end"]))])
        boundary += lo < 1
        arr[r["seqID"]].append((max(0, lo - 1), hi))
    for r in csv.DictReader((R / t / "monomer_intervals.tsv").open(), delimiter="\t"):
        pair = (int(r["start0"]), int(r["end"]))
        mono[r["chromosome"]].append(pair)
        if r["published_homology"] == "True":
            fam[r["chromosome"]].append(pair)
    arr = {k: merge(v) for k, v in arr.items()}
    mono = {k: merge(v) for k, v in mono.items()}
    fam = {k: merge(v) for k, v in fam.items()}
    cen = defaultdict(list)
    for r in csv.DictReader(
        (args.centromeres / t / "pericentromere_regions.tsv").open(), delimiter="\t"
    ):
        cen[r["chrom"]].append((int(r["start"]), int(r["end"])))
    cen = {k: merge(v) for k, v in cen.items()}
    amb = size(intersect(arr, mono))
    afb = size(intersect(arr, fam))
    cb = size(cen)
    fb = size(fam)
    gb = gs[t]
    assert 0 <= afb <= amb <= size(arr) <= gb
    r = {
        **raw,
        "genome_bp": gb,
        "array_union_bp": size(arr),
        "monomer_ge50_union_bp": size(mono),
        "family_monomer_union_bp": fb,
        "family_bp_within_arrays": afb,
        "other_ge50_monomer_bp_within_arrays": amb - afb,
        "remaining_array_bp": size(arr) - amb,
        "family_monomer_bp_outside_arrays": fb - afb,
        "pericentromere_proxy_bp": cb,
        "family_bp_in_pericentromere_proxy": size(intersect(fam, cen)),
        "array_bp_in_pericentromere_proxy": size(intersect(arr, cen)),
        "array_records_clipped_at_chr_start": boundary,
    }
    r["family_pericentromere_fraction"] = (
        r["family_bp_in_pericentromere_proxy"] / fb if fb else 0
    )
    r["family_proxy_enrichment"] = (
        r["family_pericentromere_fraction"] / (cb / gb) if cb else 0
    )
    rs.append(r)
    with (O / "census_validated.tsv").open("w") as f:
        wr = csv.DictWriter(f, fieldnames=list(r), delimiter="\t")
        wr.writeheader()
        wr.writerows(rs)
    print(
        t,
        "family union",
        fb,
        "outside arrays",
        fb - afb,
        "proxy enrichment",
        r["family_proxy_enrichment"],
        flush=True,
    )
(O / "CENSUS_VALIDATED.json").write_text(
    json.dumps(
        {
            "targets": len(rs),
            "unit": "coordinate-union bp; summed occurrence bp retained separately",
            "array_decomposition": "homology-supported family within arrays + other >=50 bp monomers within arrays + remaining array sequence; the last includes short monomers and gaps",
            "family_definition": "163-170 bp monomers with >=70% identity and >=100 aligned query bases to published references, E<=1e-5; conservative operational subset",
            "boundary_rule": "clip start to chromosome start; sequence length validates effective monomer span; raw files unchanged",
        },
        indent=2,
    )
)
