#!/usr/bin/env python3
"""Convert Quota Alignment raw hits to gene-order DAGChainer input.

The upstream ``blast_to_raw.py`` step applies the same tandem-duplicate
filter used by CoGe SynMap and writes coordinates against its sorted,
tandem-reduced BED files.  This adapter restores stable gene identifiers while
leaving the gene-order coordinates and capped -log10(E) scores unchanged.
"""

from __future__ import print_function

import argparse
import csv
from pathlib import Path


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
                "seqid": fields[0],
                "start": int(fields[1]),
                "end": int(fields[2]),
                "gene": fields[3],
            })
    rows.sort(key=lambda row: (row["seqid"], row["start"], row["gene"]))
    if not rows:
        raise ValueError("Empty BED: {}".format(path))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw", required=True, help="blast_to_raw.py .raw output")
    parser.add_argument("--query-bed", required=True, help="query8.nolocaldups.bed")
    parser.add_argument("--target-bed", required=True, help="GTR.nolocaldups.bed")
    parser.add_argument("--out", required=True)
    parser.add_argument("--stats", required=True)
    args = parser.parse_args()

    query = read_bed_in_quota_order(args.query_bed)
    target = read_bed_in_quota_order(args.target_bed)
    written = 0
    chromosome_pairs = set()
    seen = set()
    with Path(args.raw).open("r", encoding="utf-8", errors="replace") as source, \
            Path(args.out).open("w", encoding="utf-8", newline="\n") as output:
        for line_number, line in enumerate(source, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 5:
                raise ValueError("Raw row must have five columns at {}:{}".format(args.raw, line_number))
            qseqid, qindex_text, tseqid, tindex_text, score_text = fields
            qindex = int(qindex_text)
            tindex = int(tindex_text)
            score = float(score_text)
            if not 0 <= qindex < len(query) or not 0 <= tindex < len(target):
                raise ValueError("Gene-order index outside BED bounds at raw line {}".format(line_number))
            qrow = query[qindex]
            trow = target[tindex]
            if qrow["seqid"] != qseqid or trow["seqid"] != tseqid:
                raise ValueError(
                    "Raw/BED chromosome mismatch at line {}: {}:{} versus {}:{}".format(
                        line_number, qseqid, qrow["seqid"], tseqid, trow["seqid"]
                    )
                )
            key = (qrow["gene"], trow["gene"])
            if key in seen:
                continue
            seen.add(key)
            evalue = 10.0 ** (-min(max(score, 0.0), 300.0))
            output.write("\t".join([
                qseqid, qrow["gene"], str(qindex), str(qindex),
                tseqid, trow["gene"], str(tindex), str(tindex),
                "{:.6e}".format(evalue),
            ]) + "\n")
            chromosome_pairs.add((qseqid, tseqid))
            written += 1

    if written == 0:
        raise SystemExit("ERROR: no DAGChainer rows written")
    with Path(args.stats).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerow(["query_genes_after_tandem_filter", len(query)])
        writer.writerow(["target_genes_after_tandem_filter", len(target)])
        writer.writerow(["dag_anchor_rows", written])
        writer.writerow(["query_target_chromosome_pairs", len(chromosome_pairs)])
    print("[OK] Wrote {} gene-order DAGChainer anchors: {}".format(written, args.out))


if __name__ == "__main__":
    main()
