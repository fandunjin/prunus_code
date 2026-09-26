#!/usr/bin/env python3
import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


def parse_gff3(path):
    counts = Counter()
    genes_by_seq = Counter()
    lengths = defaultdict(int)
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            seqid, source, feature, start, end, score, strand, phase, attrs = fields
            counts[feature] += 1
            if feature == "gene":
                genes_by_seq[seqid] += 1
            try:
                lengths[feature] += int(end) - int(start) + 1
            except ValueError:
                pass
    return counts, genes_by_seq, lengths


def parse_args():
    parser = argparse.ArgumentParser(description="Summarize feature counts in a GFF3 file.")
    parser.add_argument("--gff3", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--out-dir", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    gff3 = Path(args.gff3)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    counts, genes_by_seq, lengths = parse_gff3(gff3)

    with (out_dir / "{}.gff3_feature_counts.tsv".format(args.prefix)).open("w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["feature", "count", "total_bp"])
        for feature in sorted(counts):
            writer.writerow([feature, counts[feature], lengths.get(feature, 0)])

    with (out_dir / "{}.genes_by_sequence.tsv".format(args.prefix)).open("w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["seq_id", "gene_count"])
        for seq_id in sorted(genes_by_seq):
            writer.writerow([seq_id, genes_by_seq[seq_id]])


if __name__ == "__main__":
    main()
