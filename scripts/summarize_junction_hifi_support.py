#!/usr/bin/env python3
"""Summarize primary HiFi read support across candidate PZI junctions.

The BAM must be aligned to the *same copy1 reference and coordinate system*
used to make the candidate windows.  A read that spans the candidate boundary
with the requested aligned sequence on both sides is recorded as support for
manual review; this script intentionally does not call a rearrangement final.
"""

from __future__ import print_function

import argparse
import csv
from pathlib import Path

import pysam


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--windows", required=True,
                        help="junction_100kb_review_windows.bed from stage 4")
    parser.add_argument("--bam", required=True,
                        help="Coordinate-sorted primary HiFi alignment BAM against copy1")
    parser.add_argument("--bai", default="",
                        help="Optional BAM index when it is not BAM.bai")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--min-flank-bp", type=int, default=1000,
                        help="Minimum aligned reference bases required on each side of the boundary")
    parser.add_argument("--min-alignment-bp", type=int, default=5000)
    parser.add_argument("--min-spanning-reads", type=int, default=3,
                        help="Project review threshold; not a final structural call")
    return parser.parse_args()


def read_windows(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"copy1_bed_chr", "copy1_fasta_chr", "start", "end", "junction_id", "boundary_mid_bp"}
    if not rows:
        return []
    missing = required.difference(rows[0])
    if missing:
        raise SystemExit("Window table is missing columns: {}".format(",".join(sorted(missing))))
    return rows


def write_tsv(path, fields, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_identity(record):
    aligned = record.query_alignment_length or 0
    if not aligned or not record.has_tag("NM"):
        return "NA"
    identity = 1.0 - float(record.get_tag("NM")) / float(aligned)
    return "{:.6f}".format(max(0.0, identity))


def main():
    args = parse_args()
    windows = read_windows(args.windows)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    index_filename = args.bai if args.bai else None
    bam = pysam.AlignmentFile(args.bam, "rb", index_filename=index_filename)
    references = set(bam.references)
    summary = []
    read_rows = []
    for window in windows:
        fasta_chrom = window["copy1_fasta_chr"]
        bed_chrom = window["copy1_bed_chr"]
        if fasta_chrom not in references:
            raise SystemExit(
                "BAM has no sequence named {}. It is not aligned to the required copy1 coordinate system.".format(fasta_chrom))
        start = int(float(window["start"]))
        end = int(float(window["end"]))
        boundary = int(float(window["boundary_mid_bp"]))
        total_primary = 0
        mapq_pass = 0
        length_pass = 0
        spanning = 0
        identities = []
        for record in bam.fetch(fasta_chrom, start, end):
            if record.is_unmapped or record.is_secondary or record.is_supplementary:
                continue
            if record.is_duplicate or record.is_qcfail:
                continue
            total_primary += 1
            if record.mapping_quality < args.min_mapq:
                continue
            mapq_pass += 1
            aligned = record.query_alignment_length or 0
            if aligned < args.min_alignment_bp:
                continue
            length_pass += 1
            ref_start = record.reference_start
            ref_end = record.reference_end
            if ref_start is None or ref_end is None:
                continue
            if ref_start > boundary - args.min_flank_bp or ref_end < boundary + args.min_flank_bp:
                continue
            spanning += 1
            identity = read_identity(record)
            if identity != "NA":
                identities.append(float(identity))
            read_rows.append({
                "junction_id": window["junction_id"], "copy1_bed_chr": bed_chrom, "copy1_fasta_chr": fasta_chrom,
                "boundary_mid_bp": boundary, "read_id": record.query_name,
                "mapping_quality": record.mapping_quality, "reference_start": ref_start,
                "reference_end": ref_end, "query_alignment_length": aligned,
                "strand": "-" if record.is_reverse else "+", "identity_from_NM": identity,
            })
        status = "read_supported_for_manual_review" if spanning >= args.min_spanning_reads else "read_not_adequately_supported"
        summary.append({
            "junction_id": window["junction_id"], "copy1_bed_chr": bed_chrom, "copy1_fasta_chr": fasta_chrom,
            "window_start": start, "window_end": end, "boundary_mid_bp": boundary,
            "primary_reads_overlapping_window": total_primary,
            "mapq_pass_reads": mapq_pass, "alignment_length_pass_reads": length_pass,
            "boundary_spanning_reads": spanning,
            "mean_identity_from_NM": "{:.6f}".format(sum(identities) / len(identities)) if identities else "NA",
            "status": status,
            "claim_boundary": "Read support requires Hi-C/manual review and does not alone confirm a rejoined junction.",
        })
    bam.close()

    write_tsv(outdir / "junction_hifi_support.tsv", [
        "junction_id", "copy1_bed_chr", "copy1_fasta_chr", "window_start", "window_end", "boundary_mid_bp",
        "primary_reads_overlapping_window", "mapq_pass_reads", "alignment_length_pass_reads",
        "boundary_spanning_reads", "mean_identity_from_NM", "status", "claim_boundary",
    ], summary)
    write_tsv(outdir / "junction_hifi_spanning_reads.tsv", [
        "junction_id", "copy1_bed_chr", "copy1_fasta_chr", "boundary_mid_bp", "read_id", "mapping_quality",
        "reference_start", "reference_end", "query_alignment_length", "strand", "identity_from_NM",
    ], read_rows)
    write_tsv(outdir / "parameters.tsv", ["parameter", "value"], [
        {"parameter": "min_mapq", "value": args.min_mapq},
        {"parameter": "min_flank_bp", "value": args.min_flank_bp},
        {"parameter": "min_alignment_bp", "value": args.min_alignment_bp},
        {"parameter": "min_spanning_reads", "value": args.min_spanning_reads},
        {"parameter": "interpretation", "value": "project review threshold; candidate remains unconfirmed without Hi-C/manual inspection"},
    ])
    print("[OK] Reviewed {} candidate windows; {} have >= {} spanning reads.".format(
        len(summary), sum(row["status"] == "read_supported_for_manual_review" for row in summary), args.min_spanning_reads))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
