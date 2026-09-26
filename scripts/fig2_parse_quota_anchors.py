#!/usr/bin/env python3
"""Normalize Quota Align QA output into GTR-centered anchor tables.

The input must be the actual output of Quota Align after the LAST and Quota
Align Merge stages.  JCVI candidate tables are deliberately rejected unless
``--allow-candidate-input`` is explicitly supplied for a debugging audit.
"""

from __future__ import print_function

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def read_bed(path):
    records = {}
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            records[fields[3]] = {
                "seqid": fields[0],
                "start": fields[1],
                "end": fields[2],
            }
    return records


def read_bed_in_quota_order(path):
    rows = []
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#") or line.startswith("track"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                raise ValueError("BED has fewer than four fields at {}:{}".format(path, line_number))
            rows.append({
                "seqid": fields[0], "start": fields[1], "end": fields[2], "gene": fields[3]
            })
    rows.sort(key=lambda row: (row["seqid"], int(row["start"]), row["gene"]))
    return rows


def read_tsv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, fields, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def parse_qa(path):
    """Yield Quota Alignment clusters in their native five-column format."""
    block = -1
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, 1):
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                if stripped.startswith("###"):
                    block += 1
                continue
            fields = stripped.split("\t")
            if len(fields) != 5:
                raise ValueError("Expected five-column QA row at {}:{}".format(path, line_number))
            if block < 0:
                raise ValueError("QA anchor precedes the first cluster header at {}:{}".format(path, line_number))
            yield block, fields[0], int(fields[1]), fields[2], int(fields[3]), fields[4]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quota-output", required=True)
    parser.add_argument("--catalog", required=True, help="Directory from fig2_build_fractionation_catalogs.py")
    parser.add_argument("--query-bed", required=True, help="Tandem-reduced query BED used by Quota Align")
    parser.add_argument("--target-bed", required=True, help="Tandem-reduced GTR BED used by Quota Align")
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()
    quota = Path(args.quota_output)
    catalog = Path(args.catalog)
    outdir = Path(args.outdir)
    if not quota.is_file() or quota.stat().st_size == 0:
        raise SystemExit("ERROR: missing Quota Align output: {}".format(quota))
    if "alltaxa_fractionation_states" in str(quota) or "jcvi" in quota.name.lower():
        raise SystemExit("ERROR: candidate JCVI fractionation table cannot be formal Quota Align input")
    if outdir.exists() and any(outdir.iterdir()):
        raise SystemExit("ERROR: output directory is not empty: {}".format(outdir))
    outdir.mkdir(parents=True, exist_ok=True)
    query_order = read_bed_in_quota_order(args.query_bed)
    target_order = read_bed_in_quota_order(args.target_bed)
    gtr_order = read_tsv(catalog / "gtr_gene_order.tsv")
    gtr_by_gene = {row["gtr_gene"]: row for row in gtr_order}
    track_by_gene = {}
    for row in read_tsv(catalog / "query_gene_track_map.tsv"):
        track_by_gene[row["query_gene"]] = row

    rows = []
    block_counts = defaultdict(int)
    seen = set()
    for block, left_seqid, left_index, right_seqid, right_index, score in parse_qa(quota):
        if left_seqid.startswith("GTR__") and not right_seqid.startswith("GTR__"):
            rseqid, rindex = left_seqid, left_index
            qseqid, qindex = right_seqid, right_index
        elif right_seqid.startswith("GTR__") and not left_seqid.startswith("GTR__"):
            qseqid, qindex = left_seqid, left_index
            rseqid, rindex = right_seqid, right_index
        else:
            raise SystemExit("ERROR: cannot orient GTR/query axes in Quota Align output")
        if not 0 <= qindex < len(query_order) or not 0 <= rindex < len(target_order):
            raise SystemExit("ERROR: Quota Align coordinate outside tandem-reduced BED bounds")
        query_meta = query_order[qindex]
        target_meta = target_order[rindex]
        if query_meta["seqid"] != qseqid or target_meta["seqid"] != rseqid:
            raise SystemExit("ERROR: Quota Align chromosome/coordinate does not match its BED order")
        qgene = query_meta["gene"]
        rgene = target_meta["gene"]
        if qgene not in track_by_gene or rgene not in gtr_by_gene:
            raise SystemExit("ERROR: Quota Align coordinate maps to a gene outside the frozen catalog")
        pair = (qgene, rgene)
        if pair in seen:
            continue
        seen.add(pair)
        track = track_by_gene[qgene]
        gtr_meta = gtr_by_gene[rgene]
        block_counts[(track["track_id"], block)] += 1
        rows.append({
            "block_id": "quota_block_{:08d}".format(block),
            "track_id": track["track_id"],
            "species": track["species"],
            "copy_state": track["copy_state"],
            "query_gene": qgene,
            "query_seqid": query_meta["seqid"],
            "query_start": query_meta["start"],
            "query_end": query_meta["end"],
            "gtr_gene": rgene,
            "gtr_seqid": gtr_meta["gtr_seqid"],
            "gtr_global_index": gtr_meta["gtr_global_index"],
            "gtr_start": gtr_meta["gtr_start"],
            "gtr_end": gtr_meta["gtr_end"],
            "quota_score": score,
        })
    if not rows:
        raise SystemExit("ERROR: no query8-to-GTR pairs parsed from {}".format(quota))
    rows.sort(key=lambda row: (int(row["gtr_global_index"]), row["track_id"], row["query_gene"]))
    fields = list(rows[0].keys())
    write_tsv(outdir / "quota_filtered_anchors.tsv", fields, rows)

    presence = defaultdict(set)
    for row in rows:
        presence[(row["track_id"], row["gtr_gene"])].add(row["query_gene"])
    presence_rows = []
    for gtr in gtr_order:
        for track_id in ("PZI_A", "PZI_B", "PSP_A", "PSP_B", "PPA_A", "PPA_B", "PCAM", "PPE"):
            genes = sorted(presence.get((track_id, gtr["gtr_gene"]), set()))
            presence_rows.append({
                "track_id": track_id,
                "gtr_global_index": gtr["gtr_global_index"],
                "gtr_seqid": gtr["gtr_seqid"],
                "gtr_gene": gtr["gtr_gene"],
                "anchored": 1 if genes else 0,
                "n_query_genes": len(genes),
                "query_genes": ",".join(genes),
            })
    write_tsv(outdir / "gtr_track_presence.tsv", list(presence_rows[0].keys()), presence_rows)
    write_tsv(outdir / "quota_block_summary.tsv",
              ["track_id", "block_id", "n_pairs"],
              [{"track_id": track, "block_id": block, "n_pairs": count} for (track, block), count in sorted(block_counts.items())])
    (outdir / "method_status.txt").write_text(
        "This table was parsed from LAST -> Quota Align Merge -> Quota Align output.\n"
        "Anchored=1 means a quota-filtered syntenic anchor was observed; it is not a verified gene-loss state.\n",
        encoding="ascii",
    )
    (outdir / "COMPLETE").write_text("Formal Quota Align anchor normalization completed.\n", encoding="ascii")
    print("[OK] Parsed {} quota-filtered anchors across {} GTR loci".format(len(rows), len(gtr_order)))


if __name__ == "__main__":
    main()
