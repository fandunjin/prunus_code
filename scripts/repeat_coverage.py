"""Recompute exclusive repeat bp using validated monomer membership and coordinates."""

from pathlib import Path
import csv, json, heapq
from collections import defaultdict, Counter
import repeat_classes as prs

import argparse

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--manifest", type=Path, required=True)
parser.add_argument("--census-dir", type=Path, required=True)
parser.add_argument("--out", type=Path, required=True)
args = parser.parse_args()
R = args.census_dir
O = args.out
O.mkdir(parents=True, exist_ok=True)
mp = {
    "PZI-A1": "pziA1",
    "PZI-B1": "pziB1",
    "PSP-A1": "pspA1",
    "PSP-B1": "pspB1",
    "PPA-A": "ppaA",
    "PPA-B": "ppaB",
    "PCAM": "pcam",
    "PPE": "ppe",
}
records = []
summary = []
spatial = []
for spec in csv.DictReader(args.manifest.open(), delimiter="\t"):
    target = spec["target"]
    gff = Path(spec["edta_gff"])
    fa = Path(spec["fasta"])
    source = R / mp[target] / "monomer_intervals.tsv"
    lengths = {
        r[0]: int(r[1])
        for r in csv.reader(Path(str(fa) + ".fai").open(), delimiter="\t")
    }

    def seqid(raw):
        if raw in lengths:
            return raw
        candidates = [k for k in lengths if k.startswith(raw)]
        assert len(candidates) == 1, (target, raw, candidates)
        return candidates[0]

    lin = {
        r["record_id"]: tuple(
            r[k] for k in ["class", "order", "family", "sublineage", "clade"]
        )
        for r in csv.DictReader(Path(spec["lineage_table"]).open(), delimiter="\t")
    }
    by = defaultdict(list)
    for row in csv.reader(gff.open(), delimiter="\t"):
        if not row or row[0].startswith("#") or len(row) < 9:
            continue
        s, e = int(row[3]) - 1, int(row[4])
        attr = dict(x.split("=", 1) for x in row[8].split(";") if "=" in x)
        token = attr.get("classification", row[2])
        cls = prs.normalize_class(token)
        if cls is None:
            continue
        score = float(row[5]) if row[5] not in {".", ""} else 0
        chrom = seqid(row[0])
        assert 0 <= s < e <= lengths[chrom], (target, chrom, s, e)
        by[chrom].append(
            (
                s,
                e,
                score,
                e - s,
                cls,
                prs.subfamily_of(cls, token, attr.get("Name", ""), lin),
            )
        )
    for row in csv.DictReader(source.open(), delimiter="\t"):
        s, e = int(row["start0"]), int(row["end"])
        yes = row["published_homology"] == "True"
        assert row["chromosome"] in lengths and 0 <= s < e <= lengths[row["chromosome"]]
        by[row["chromosome"]].append(
            (
                s,
                e,
                1e6,
                e - s,
                "Satellite",
                (
                    "Prunus_170 homologues (163-170 bp)"
                    if yes
                    else "other tandem repeats (>=50 bp)"
                ),
            )
        )
        if yes:
            for binidx in range(s // 1000000, (e - 1) // 1000000 + 1):
                spatial.append(
                    (
                        target,
                        row["chromosome"],
                        binidx,
                        max(
                            0, min(e, (binidx + 1) * 1000000) - max(s, binidx * 1000000)
                        ),
                    )
                )
    coverage = Counter()
    for chrom, iv in by.items():
        events = defaultdict(list)
        for i, (s, e, score, span, *_) in enumerate(iv):
            events[s].append((1, i))
            events[e].append((-1, i))
        active = set()
        heap = []
        prev = None
        for x, ev in sorted(events.items()):
            while heap and heap[0][2] not in active:
                heapq.heappop(heap)
            if prev is not None and x > prev and heap:
                row = iv[heap[0][2]]
                coverage[(row[4], row[5])] += x - prev
            for sign, i in ev:
                if sign == -1:
                    active.discard(i)
                else:
                    active.add(i)
                    heapq.heappush(heap, (-iv[i][2], -iv[i][3], i))
            prev = x
    total = sum(lengths.values())
    assert sum(coverage.values()) <= total
    for (cls, sub), bp in coverage.items():
        records.append(
            dict(target=target, repeat_class=cls, subfamily=sub, bp=bp, genome_bp=total)
        )
    summary.append(
        dict(
            target=target,
            genome_bp=total,
            exclusive_repeat_bp=sum(coverage.values()),
            family_bp=sum(
                v for (c, s), v in coverage.items() if s.startswith("Prunus_170")
            ),
            other_tandem_bp=sum(
                v for (c, s), v in coverage.items() if s.startswith("other tandem")
            ),
        )
    )
    print(summary[-1], flush=True)
    for filename, rs in [
        ("exclusive_composition.tsv", records),
        ("exclusive_summary.tsv", summary),
    ]:
        with (O / filename).open("w") as f:
            wr = csv.DictWriter(f, fieldnames=list(rs[0]), delimiter="\t")
            wr.writeheader()
            wr.writerows(rs)
sp = Counter()
for t, c, b, n in spatial:
    sp[t, c, b] += n
with (O / "family_spatial_1Mb.tsv").open("w") as f:
    wr = csv.writer(f, delimiter="\t")
    wr.writerow(["target", "chromosome", "start0", "end", "homology_monomer_bp"])
    for (t, c, b), n in sorted(sp.items()):
        wr.writerow([t, c, b * 1000000, (b + 1) * 1000000, n])
(O / "COVERAGE_DONE.json").write_text(
    json.dumps(
        {
            "targets": len(summary),
            "coordinate_convention": "TRASH 1-based inclusive converted to 0-based half-open",
            "priority": "TRASH monomers >=50 bp, then EDTA maximum score/longest span; input order breaks remaining ties",
            "genome_denominator": "full FASTA index length; EDTA truncated seqids uniquely mapped to complete FASTA names",
            "family_definition": "163-170 bp monomer with published-reference homology; conservative operational census, not all possible family variants",
        },
        indent=2,
    )
)
