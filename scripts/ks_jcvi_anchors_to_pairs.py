#!/usr/bin/env python3
"""Convert screened JCVI anchors into standardized syntenic-pair tables."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pair-plan", required=True)
    parser.add_argument("--pair-root", required=True)
    parser.add_argument("--out-pairs", required=True)
    parser.add_argument("--out-blocks", required=True)
    parser.add_argument("--out-status", required=True)
    return parser.parse_args()


def read_tsv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_bed(path):
    records = {}
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                continue
            records[fields[3]] = {
                "seqid": fields[0],
                "start0": start,
                "end": end,
                "strand": fields[5] if len(fields) > 5 else ".",
            }
    return records


def parse_anchors(path):
    blocks = defaultdict(list)
    block_index = 1
    block_has_pairs = False
    seen = set()
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                if block_has_pairs:
                    block_index += 1
                    block_has_pairs = False
                continue
            if stripped.startswith("###"):
                if block_has_pairs:
                    block_index += 1
                    block_has_pairs = False
                continue
            if stripped.startswith("#"):
                continue
            fields = stripped.split()
            if len(fields) < 2:
                continue
            pair = (fields[0], fields[1])
            if pair in seen:
                continue
            seen.add(pair)
            blocks[block_index].append(pair)
            block_has_pairs = True
    return blocks


def write_tsv(path, rows, fields):
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    args = parse_args()
    pair_root = Path(args.pair_root)
    pair_rows = []
    block_rows = []
    status_rows = []

    for plan in read_tsv(args.pair_plan):
        query = plan["query_label"]
        ref = plan["ref_label"]
        comparison = plan.get("comparison") or "{}_vs_{}".format(query, ref)
        comparison_type = plan.get("comparison_type", "")
        work = pair_root / comparison
        anchors = work / "{}.{}.anchors.new".format(query, ref)
        qbed_path = work / "{}.bed".format(query)
        rbed_path = work / "{}.bed".format(ref)
        if not anchors.is_file() or anchors.stat().st_size == 0:
            status_rows.append({
                "comparison": comparison,
                "status": "missing_screened_anchors",
                "n_blocks": 0,
                "n_pairs": 0,
                "anchors": anchors,
            })
            continue
        qbed = read_bed(qbed_path)
        rbed = read_bed(rbed_path)
        blocks = parse_anchors(anchors)
        n_pairs = 0
        for source_block, pairs in sorted(blocks.items()):
            if not pairs:
                continue
            block_id = "{}_block{:06d}".format(comparison, source_block)
            qseqids = set()
            rseqids = set()
            for qgene, rgene in pairs:
                qrec = qbed.get(qgene, {})
                rrec = rbed.get(rgene, {})
                if qrec.get("seqid"):
                    qseqids.add(qrec["seqid"])
                if rrec.get("seqid"):
                    rseqids.add(rrec["seqid"])
                pair_rows.append({
                    "comparison": comparison,
                    "comparison_type": comparison_type,
                    "block_id": block_id,
                    "query_label": query,
                    "ref_label": ref,
                    "query_gene": qgene,
                    "ref_gene": rgene,
                    "query_seqid": qrec.get("seqid", ""),
                    "query_start0": qrec.get("start0", ""),
                    "query_end": qrec.get("end", ""),
                    "query_strand": qrec.get("strand", ""),
                    "ref_seqid": rrec.get("seqid", ""),
                    "ref_start0": rrec.get("start0", ""),
                    "ref_end": rrec.get("end", ""),
                    "ref_strand": rrec.get("strand", ""),
                })
                n_pairs += 1
            block_rows.append({
                "comparison": comparison,
                "comparison_type": comparison_type,
                "block_id": block_id,
                "query_label": query,
                "ref_label": ref,
                "n_pairs": len(pairs),
                "query_seqids": ",".join(sorted(qseqids)),
                "ref_seqids": ",".join(sorted(rseqids)),
            })
        status_rows.append({
            "comparison": comparison,
            "status": "ok" if n_pairs else "no_pairs_parsed",
            "n_blocks": sum(1 for pairs in blocks.values() if pairs),
            "n_pairs": n_pairs,
            "anchors": anchors,
        })

    write_tsv(
        args.out_pairs,
        pair_rows,
        [
            "comparison", "comparison_type", "block_id", "query_label", "ref_label",
            "query_gene", "ref_gene", "query_seqid", "query_start0", "query_end",
            "query_strand", "ref_seqid", "ref_start0", "ref_end", "ref_strand",
        ],
    )
    write_tsv(
        args.out_blocks,
        block_rows,
        [
            "comparison", "comparison_type", "block_id", "query_label", "ref_label",
            "n_pairs", "query_seqids", "ref_seqids",
        ],
    )
    write_tsv(
        args.out_status,
        status_rows,
        ["comparison", "status", "n_blocks", "n_pairs", "anchors"],
    )
    if not pair_rows or any(row["status"] != "ok" for row in status_rows):
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
